#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <utility>
#include <vector>

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLockFile>
#include <QSaveFile>
#include <QStringList>
#include <QTemporaryDir>

#include <rbl_application_state.h>
#include <rbl_error.h>

#include <rml_file_utils.h>
#include <rml_structural_product_profile.h>

#include "rsolver.h"
#include "rsolveracoustic.h"
#include "rsolverfluidparticle.h"
#include "rsolverelectrostatics.h"
#include "rsolvermagnetostatics.h"
#include "rsolverfluid.h"
#include "rsolverfluidheat.h"
#include "rsolverheat.h"
#include "rsolvermesh.h"
#include "rsolverradiativeheat.h"
#include "rsolverstress.h"
#include "rstructural_cable_model_adapter.h"
#include "rstructural_piezo_model_adapter.h"
#include "rparallel_execution_context.h"
#include "rsolverwave.h"

namespace
{

bool isHeatStressSequence(const RProblemTaskItem &task)
{
    if (task.getProblemType() != R_PROBLEM_NONE || task.getNChildren() != 2)
    {
        return false;
    }

    const RProblemTaskItem &heat = task.getChild(0);
    const RProblemTaskItem &stress = task.getChild(1);
    return heat.getProblemType() == R_PROBLEM_HEAT && heat.getNChildren() == 0 &&
           stress.getProblemType() == R_PROBLEM_STRESS && stress.getNChildren() == 0;
}

QString getCanonicalOutputLinkPath(const QString &modelFileName)
{
    const QString linkFileName = RFileUtils::getFileNameWithOutTimeStep(modelFileName);
    const QFileInfo linkInfo(linkFileName);
    if (linkInfo.isSymLink())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transient Heat->Stress output cannot use a symbolic-link target '%s'.",
                     linkFileName.toUtf8().constData());
    }
    QString directory = QDir(linkInfo.absolutePath()).canonicalPath();
    if (directory.isEmpty())
    {
        directory = QDir::cleanPath(linkInfo.absolutePath());
    }
    return QDir(directory).filePath(linkInfo.fileName());
}

QString getCanonicalFilePath(const QString &fileName)
{
    const QFileInfo info(fileName);
    if (info.isSymLink())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transactional transient output cannot replace a symbolic-link file '%s'.",
                     fileName.toUtf8().constData());
    }
    QString directory = QDir(info.absolutePath()).canonicalPath();
    if (directory.isEmpty())
    {
        directory = QDir::cleanPath(info.absolutePath());
    }
    return QDir(directory).filePath(info.fileName());
}

QString getTransientExactTargetLockPath(const QString &canonicalTarget)
{
    return canonicalTarget+QStringLiteral(".transient.lock");
}

class TransientOutputLock
{
private:
    std::unique_ptr<QLockFile> lockFile;

public:
    explicit TransientOutputLock(const QString &modelFileName)
    {
        if (modelFileName.isEmpty())
        {
            return;
        }
        const QString canonicalLinkPath = getCanonicalOutputLinkPath(modelFileName);
        const QString lockPath = getTransientExactTargetLockPath(canonicalLinkPath);
        if (QFileInfo(lockPath).isSymLink())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transient Heat->Stress output lock path is a symbolic link: '%s'.",
                         lockPath.toUtf8().constData());
        }
        this->lockFile.reset(new QLockFile(lockPath));
        // A slow solve must not lose its lock just because it exceeds a
        // timeout. A leftover lock after process termination is handled by
        // failing closed for manual inspection.
        this->lockFile->setStaleLockTime(0);
        if (!this->lockFile->tryLock(0))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transient Heat->Stress output is locked by another process or a stale lock: '%s'.",
                         lockPath.toUtf8().constData());
        }
    }
};

bool isOutputStepScheduled(const RTimeSolver &timeSolver, uint outputStep)
{
    return outputStep == timeSolver.getNTimeSteps() ||
           (timeSolver.getOutputFrequency() > 0 &&
            outputStep % timeSolver.getOutputFrequency() == 0);
}

bool areOutputPathsEqual(const QString &first, const QString &second)
{
#ifdef Q_OS_WIN
    return first.compare(second,Qt::CaseInsensitive) == 0;
#else
    return first == second;
#endif
}

Qt::CaseSensitivity outputPathCaseSensitivity()
{
#ifdef Q_OS_WIN
    return Qt::CaseInsensitive;
#else
    return Qt::CaseSensitive;
#endif
}

bool isTransientModelOutputPath(const QString &canonicalTarget,
                                const QString &modelFileName,
                                const RTimeSolver &timeSolver)
{
    if (modelFileName.isEmpty())
    {
        return false;
    }
    const QString modelLink = getCanonicalOutputLinkPath(modelFileName);
    if (areOutputPathsEqual(canonicalTarget,modelLink))
    {
        return true;
    }

    const QFileInfo linkInfo(modelLink);
    const QFileInfo targetInfo(canonicalTarget);
    if (!areOutputPathsEqual(targetInfo.absolutePath(),linkInfo.absolutePath()))
    {
        return false;
    }
    const QString prefix = linkInfo.completeBaseName()+QStringLiteral("-");
    const QString suffix = linkInfo.suffix().isEmpty()
        ? QString() : QStringLiteral(".")+linkInfo.suffix();
    const QString targetName = targetInfo.fileName();
    if (!targetName.startsWith(prefix,outputPathCaseSensitivity()) ||
        !targetName.endsWith(suffix,outputPathCaseSensitivity()) ||
        targetName.size() < prefix.size()+5+suffix.size())
    {
        return false;
    }
    const QString stepText = targetName.mid(prefix.size(),
        targetName.size()-prefix.size()-suffix.size());
    if (stepText.size() < 5)
    {
        return false;
    }
    for (const QChar digit : stepText)
    {
        if (!digit.isDigit())
        {
            return false;
        }
    }
    bool parsed = false;
    const uint outputStep = stepText.toUInt(&parsed);
    if (!parsed || outputStep == 0 || outputStep > timeSolver.getNTimeSteps() ||
        QStringLiteral("%1").arg(outputStep,5,10,QChar('0')) != stepText)
    {
        return false;
    }
    return isOutputStepScheduled(timeSolver,outputStep);
}

class TransientMonitoringOutputRunLock
{
private:
    std::vector<std::unique_ptr<QLockFile>> locks;

public:
    TransientMonitoringOutputRunLock(const RMonitoringPointManager &monitoringPoints,
                                     const QString &modelFileName,
                                     const RTimeSolver &timeSolver)
    {
        std::vector<QString> targetFileNames;
        for (const RMonitoringPoint &point : monitoringPoints)
        {
            if (point.getOutputFileName().isEmpty())
            {
                continue;
            }
            const QString canonicalFileName = getCanonicalFilePath(point.getOutputFileName());
            if (isTransientModelOutputPath(canonicalFileName,modelFileName,timeSolver))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Transient monitoring output cannot target a model result link or scheduled time-step record '%s'.",
                             canonicalFileName.toUtf8().constData());
            }
            if (std::none_of(targetFileNames.cbegin(),targetFileNames.cend(),
                             [&canonicalFileName](const QString &existingFileName)
                             {
                                 return areOutputPathsEqual(existingFileName,canonicalFileName);
                             }))
            {
                targetFileNames.push_back(canonicalFileName);
            }
        }

        std::sort(targetFileNames.begin(),targetFileNames.end(),
                  [](const QString &first,const QString &second)
                  {
                      return first.compare(second,outputPathCaseSensitivity()) < 0;
                  });

        this->locks.reserve(targetFileNames.size()*2);
        for (const QString &targetFileName : targetFileNames)
        {
            const QStringList lockPaths = {
                getTransientExactTargetLockPath(targetFileName),
                // Honor lock files left by a process running the preceding
                // monitor-transaction implementation.
                targetFileName+QStringLiteral(".transient-monitor.lock")
            };
            for (const QString &lockPath : lockPaths)
            {
                if (QFileInfo(lockPath).isSymLink())
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Transient monitoring output lock path is a symbolic link: '%s'.",
                                 lockPath.toUtf8().constData());
                }
                std::unique_ptr<QLockFile> lock(new QLockFile(lockPath));
                lock->setStaleLockTime(0);
                if (!lock->tryLock(0))
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Transient monitoring output is locked by another process or a stale lock: '%s'.",
                                 lockPath.toUtf8().constData());
                }
                this->locks.push_back(std::move(lock));
            }
        }
    }
};

void rejectUnpublishedOutputRecords(const QString &modelFileName,
                                    const RTimeSolver &timeSolver,
                                    bool restart)
{
    if (modelFileName.isEmpty() || timeSolver.getNTimeSteps() == 0)
    {
        return;
    }

    const QString linkFileName = RFileUtils::getFileNameWithOutTimeStep(modelFileName);
    const QFileInfo linkInfo(linkFileName);
    const QString prefix = linkInfo.completeBaseName()+QStringLiteral("-");
    const QString suffix = QStringLiteral(".")+linkInfo.suffix();
    const QDir outputDirectory(linkInfo.absolutePath());

    const uint acceptedOutputStep = restart ? timeSolver.getCurrentTimeStep()+1 : 0;
    const uint nextTrialOutputStep = acceptedOutputStep+1;
    const bool nextTrialWillPublish = nextTrialOutputStep <= timeSolver.getNTimeSteps() &&
                                      isOutputStepScheduled(timeSolver,nextTrialOutputStep);
    const QStringList entries = outputDirectory.entryList(
        QDir::Files | QDir::Hidden | QDir::System | QDir::NoDotAndDotDot,
        QDir::NoSort);
    for (const QString &entry : entries)
    {
        if (!entry.startsWith(prefix) || !entry.endsWith(suffix) ||
            entry.size() < prefix.size()+5+suffix.size())
        {
            continue;
        }
        const QString stepText = entry.mid(prefix.size(),
                                           entry.size()-prefix.size()-suffix.size());
        if (stepText.size() < 5)
        {
            continue;
        }
        bool allDigits = true;
        for (const QChar digit : stepText)
        {
            if (!digit.isDigit())
            {
                allDigits = false;
                break;
            }
        }
        if (!allDigits)
        {
            continue;
        }
        bool parsed = false;
        const uint outputStep = stepText.toUInt(&parsed);
        if (!parsed || outputStep <= acceptedOutputStep ||
            outputStep > timeSolver.getNTimeSteps())
        {
            continue;
        }
        if (nextTrialWillPublish && outputStep == nextTrialOutputStep)
        {
            continue;
        }
        const QString stalePath = outputDirectory.filePath(entry);
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transient Heat->Stress run found a pre-existing output record '%s' beyond accepted history; remove or archive it before continuing.",
                     stalePath.toUtf8().constData());
    }
}

class OutputFileSnapshot
{
private:
    enum class OriginalState { Absent, RegularFile, Other };

    QString path;
    QString backupPath;
    OriginalState originalState = OriginalState::Absent;

    bool matchesBackup() const noexcept
    {
        const QFileInfo originalInfo(this->path);
        const QFileInfo backupInfo(this->backupPath);
        if (originalInfo.isSymLink() || !originalInfo.isFile() || !backupInfo.isFile() ||
            originalInfo.size() != backupInfo.size())
        {
            return false;
        }
        QFile original(this->path);
        QFile backup(this->backupPath);
        if (!original.open(QIODevice::ReadOnly) || !backup.open(QIODevice::ReadOnly))
        {
            return false;
        }
        while (!original.atEnd() && !backup.atEnd())
        {
            const QByteArray originalBlock = original.read(1024*1024);
            const QByteArray backupBlock = backup.read(1024*1024);
            if (original.error() != QFileDevice::NoError || backup.error() != QFileDevice::NoError ||
                originalBlock != backupBlock)
            {
                return false;
            }
        }
        return original.atEnd() && backup.atEnd();
    }

    bool restoreFile() const noexcept
    {
        if (this->originalState == OriginalState::Other)
        {
            return true;
        }

        if (this->originalState == OriginalState::Absent)
        {
            const QFileInfo current(this->path);
            return !current.exists() || !current.isFile() || QFile::remove(this->path);
        }
        if (this->matchesBackup())
        {
            return true;
        }

        QFile source(this->backupPath);
        if (!source.open(QIODevice::ReadOnly))
        {
            return false;
        }
        QSaveFile destination(this->path);
        if (!destination.open(QIODevice::WriteOnly))
        {
            return false;
        }
        while (!source.atEnd())
        {
            const QByteArray block = source.read(1024*1024);
            if (source.error() != QFileDevice::NoError ||
                (!block.isEmpty() && destination.write(block) != block.size()))
            {
                destination.cancelWriting();
                return false;
            }
        }
        return destination.commit();
    }

public:
    void capture(const QString &fileName, const QString &backupFileName)
    {
        this->path = fileName;
        this->backupPath = backupFileName;
        const QFileInfo info(fileName);
        if (info.isSymLink())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transactional transient output cannot replace a symbolic-link file '%s'.",
                         fileName.toUtf8().constData());
        }
        if (!info.exists())
        {
            this->originalState = OriginalState::Absent;
            return;
        }
        if (!info.isFile())
        {
            // Directories and other non-file paths cannot be replaced by the
            // model writer. Preserve them as-is and let publication report it.
            this->originalState = OriginalState::Other;
            return;
        }
        this->originalState = OriginalState::RegularFile;
        if (!QFile::copy(fileName,backupFileName))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Could not preserve existing transient output file '%s'.",
                         fileName.toUtf8().constData());
        }
    }

    bool restore() const noexcept
    {
        if (!this->restoreFile())
        {
            return false;
        }
        return true;
    }

    const QString &getPath() const noexcept
    {
        return this->path;
    }
};

class DeferredOutputTransaction
{
private:
    std::unique_ptr<QTemporaryDir> backupDirectory;
    std::unique_ptr<QLockFile> recordLock;
    OutputFileSnapshot record;
    OutputFileSnapshot link;
    bool prepared = false;
    bool committed = false;

public:
    void prepare(const QString &modelFileName, const RTimeSolver &timeSolver)
    {
        if (modelFileName.isEmpty() || !timeSolver.getEnabled())
        {
            return;
        }
        const uint outputStep = timeSolver.getCurrentTimeStep()+1;
        const bool shouldPublish = outputStep == timeSolver.getNTimeSteps() ||
            (timeSolver.getOutputFrequency() > 0 &&
             outputStep % timeSolver.getOutputFrequency() == 0);
        if (!shouldPublish)
        {
            return;
        }

        this->backupDirectory.reset(new QTemporaryDir(
            QDir::tempPath()+QStringLiteral("/range-fea-step-XXXXXX")));
        if (!this->backupDirectory->isValid())
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Could not create a temporary directory for transient output rollback.");
        }
        const QString recordFileName = RFileUtils::getFileNameWithTimeStep(modelFileName,outputStep);
        const QString linkFileName = RFileUtils::getFileNameWithOutTimeStep(modelFileName);
        const QString canonicalRecordFileName = getCanonicalFilePath(recordFileName);
        const QString recordLockPath = getTransientExactTargetLockPath(canonicalRecordFileName);
        if (QFileInfo(recordLockPath).isSymLink())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transient model result record lock path is a symbolic link: '%s'.",
                         recordLockPath.toUtf8().constData());
        }
        this->recordLock.reset(new QLockFile(recordLockPath));
        this->recordLock->setStaleLockTime(0);
        if (!this->recordLock->tryLock(0))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transient model result record is locked by another run or a stale lock: '%s'.",
                         canonicalRecordFileName.toUtf8().constData());
        }
        this->record.capture(recordFileName,this->backupDirectory->filePath(QStringLiteral("record")));
        this->link.capture(linkFileName,this->backupDirectory->filePath(QStringLiteral("link")));
        this->prepared = true;
    }

    void commit() noexcept
    {
        this->committed = true;
    }

    QString rollback() noexcept
    {
        if (!this->prepared || this->committed)
        {
            return QString();
        }
        const bool recordRestored = this->record.restore();
        const bool linkRestored = this->link.restore();
        if (!recordRestored || !linkRestored)
        {
            QStringList failedPaths;
            if (!recordRestored)
            {
                failedPaths.append(this->record.getPath());
            }
            if (!linkRestored)
            {
                failedPaths.append(this->link.getPath());
            }
            this->backupDirectory->setAutoRemove(false);
            this->committed = true;
            return QStringLiteral("Could not restore %1. Rollback backups were retained at '%2'.")
                .arg(failedPaths.join(QStringLiteral(", ")),
                     this->backupDirectory->path());
        }
        this->committed = true;
        return QString();
    }
};

class TransientMonitoringOutputTransaction
{
private:
    struct Target
    {
        QString fileName;
        QString stagingFileName;
        OutputFileSnapshot snapshot;
    };

    RMonitoringPointManager *manager = nullptr;
    std::vector<QString> originalFileNames;
    std::vector<int> pointTargetIndices;
    std::vector<Target> targets;
    std::unique_ptr<QTemporaryDir> temporaryDirectory;
    bool pathsRestored = false;
    bool prepared = false;
    bool committed = false;

    static void copyContents(QFile &source, QSaveFile &destination,
                             const QString &sourceName, const QString &destinationName)
    {
        while (!source.atEnd())
        {
            const QByteArray block = source.read(1024*1024);
            if (source.error() != QFileDevice::NoError ||
                (!block.isEmpty() && destination.write(block) != block.size()))
            {
                throw RError(RError::Type::WriteFile,R_ERROR_REF,
                             "Could not copy transactional monitoring output '%s' to '%s'.",
                             sourceName.toUtf8().constData(),destinationName.toUtf8().constData());
            }
        }
    }

    void restoreFileNames() noexcept
    {
        if (!this->manager || this->pathsRestored)
        {
            return;
        }
        const size_t n = std::min(this->manager->size(),this->originalFileNames.size());
        for (size_t i=0;i<n;i++)
        {
            (*this->manager)[static_cast<uint>(i)].setOutputFileName(this->originalFileNames[i]);
        }
        this->pathsRestored = true;
    }

public:
    void prepare(RMonitoringPointManager &monitoringPoints,
                 const QString &modelFileName,
                 const RTimeSolver &timeSolver)
    {
        this->manager = &monitoringPoints;
        this->originalFileNames.reserve(monitoringPoints.size());
        this->pointTargetIndices.resize(monitoringPoints.size(),-1);
        for (const RMonitoringPoint &point : monitoringPoints)
        {
            this->originalFileNames.push_back(point.getOutputFileName());
        }

        for (uint i=0;i<monitoringPoints.size();i++)
        {
            const QString &requestedFileName = monitoringPoints[i].getOutputFileName();
            if (requestedFileName.isEmpty())
            {
                continue;
            }
            const QString canonicalFileName = getCanonicalFilePath(requestedFileName);
            size_t targetIndex = this->targets.size();
            for (size_t j=0;j<this->targets.size();j++)
            {
                if (areOutputPathsEqual(this->targets[j].fileName,canonicalFileName))
                {
                    targetIndex = j;
                    break;
                }
            }
            if (targetIndex == this->targets.size())
            {
                Target target;
                target.fileName = canonicalFileName;
                this->targets.push_back(std::move(target));
            }
            this->pointTargetIndices[i] = static_cast<int>(targetIndex);
        }

        if (this->targets.empty())
        {
            this->pathsRestored = true;
            return;
        }
        for (const Target &target : this->targets)
        {
            if (isTransientModelOutputPath(target.fileName,modelFileName,timeSolver))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Transient monitoring output cannot target a model result link or scheduled time-step record '%s'.",
                             target.fileName.toUtf8().constData());
            }
        }

        this->temporaryDirectory.reset(new QTemporaryDir(
            QDir::tempPath()+QStringLiteral("/range-fea-monitor-XXXXXX")));
        if (!this->temporaryDirectory->isValid())
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Could not create a temporary directory for transient monitoring output.");
        }
        for (size_t i=0;i<this->targets.size();i++)
        {
            Target &target = this->targets[i];
            const QString backupFileName = this->temporaryDirectory->filePath(
                QStringLiteral("monitor-%1-backup").arg(i));
            target.snapshot.capture(target.fileName,backupFileName);
            target.stagingFileName = this->temporaryDirectory->filePath(
                QStringLiteral("monitor-%1-stage").arg(i));
        }

        this->prepared = true;
        for (uint i=0;i<monitoringPoints.size();i++)
        {
            if (this->pointTargetIndices[i] < 0)
            {
                continue;
            }
            const size_t targetIndex = static_cast<size_t>(this->pointTargetIndices[i]);
            if (targetIndex >= this->targets.size())
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Transient monitoring output target mapping is invalid.");
            }
            monitoringPoints[i].setOutputFileName(this->targets[targetIndex].stagingFileName);
        }
    }

    void restoreOutputFileNames() noexcept
    {
        this->restoreFileNames();
    }

    void publish()
    {
        this->restoreFileNames();
        for (Target &target : this->targets)
        {
            QFile staging(target.stagingFileName);
            const QFileInfo stagingInfo(target.stagingFileName);
            if (!stagingInfo.exists() || stagingInfo.size() == 0)
            {
                continue;
            }
            if (!staging.open(QIODevice::ReadOnly))
            {
                throw RError(RError::Type::OpenFile,R_ERROR_REF,
                             "Could not read staged monitoring output '%s'.",
                             target.stagingFileName.toUtf8().constData());
            }

            const QFileInfo targetInfo(target.fileName);
            if (targetInfo.isSymLink() || (targetInfo.exists() && !targetInfo.isFile()))
            {
                throw RError(RError::Type::WriteFile,R_ERROR_REF,
                             "Transient monitoring output target is not a regular file: '%s'.",
                             target.fileName.toUtf8().constData());
            }
            QSaveFile output(target.fileName);
            if (!output.open(QIODevice::WriteOnly))
            {
                throw RError(RError::Type::OpenFile,R_ERROR_REF,
                             "Could not open transient monitoring output '%s': %s",
                             target.fileName.toUtf8().constData(),
                             output.errorString().toUtf8().constData());
            }
            if (targetInfo.exists())
            {
                QFile original(target.fileName);
                if (!original.open(QIODevice::ReadOnly))
                {
                    output.cancelWriting();
                    throw RError(RError::Type::OpenFile,R_ERROR_REF,
                                 "Could not read existing transient monitoring output '%s'.",
                                 target.fileName.toUtf8().constData());
                }
                copyContents(original,output,target.fileName,target.fileName);
            }
            copyContents(staging,output,target.stagingFileName,target.fileName);
            if (!output.commit())
            {
                throw RError(RError::Type::WriteFile,R_ERROR_REF,
                             "Could not atomically publish transient monitoring output '%s': %s",
                             target.fileName.toUtf8().constData(),
                             output.errorString().toUtf8().constData());
            }
        }
    }

    void commit() noexcept
    {
        this->committed = true;
    }

    QString rollback() noexcept
    {
        this->restoreFileNames();
        if (!this->prepared || this->committed)
        {
            return QString();
        }
        QStringList failedPaths;
        for (auto it=this->targets.rbegin();it!=this->targets.rend();++it)
        {
            if (!it->snapshot.restore())
            {
                failedPaths.append(it->snapshot.getPath());
            }
        }
        if (!failedPaths.isEmpty())
        {
            this->temporaryDirectory->setAutoRemove(false);
            this->committed = true;
            return QStringLiteral("Could not restore monitoring output %1. Rollback backups were retained at '%2'.")
                .arg(failedPaths.join(QStringLiteral(", ")),
                     this->temporaryDirectory->path());
        }
        this->committed = true;
        return QString();
    }
};

}

void RSolver::_init()
{
    RProblemTypeMask problemTypeMask = this->pModel->getProblemTaskTree().getProblemTypeMask();
    if (problemTypeMask == R_PROBLEM_NONE)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "No solver problem task is configured.");
    }

    if (this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled)
    {
        // This gate deliberately precedes solver construction and the result
        // cleanup below for every execution policy.
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
    }

    if (this->executionPolicy == ExecutionPolicy::StructuralProduct)
    {
        const RProductProfileValidation profileValidation =
            RStructuralProductProfile::validate(*this->pModel);
        if (!profileValidation.isAccepted())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         profileValidation.getDiagnostic());
        }
    }

    RRVector initialElementTemperature(this->pModel->getNElements(),
        RVariable::getInitValue(R_VARIABLE_TEMPERATURE));
    const uint temperaturePosition = this->pModel->findVariable(R_VARIABLE_TEMPERATURE);
    if (temperaturePosition != RConstants::eod)
    {
        const RVariable &temperature = this->pModel->getVariable(temperaturePosition);
        if (temperature.getNVectors() != 1)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Persisted temperature must contain exactly one absolute-kelvin component.");
        }
        if (temperature.getApplyType() == R_VARIABLE_APPLY_NODE &&
            temperature.getNValues() == this->pModel->getNNodes())
        {
            RRVector nodalTemperature(this->pModel->getNNodes());
            for (uint node=0;node<this->pModel->getNNodes();node++)
            {
                const double value = temperature.getValue(0,node);
                if (!std::isfinite(value) || value < 0.0)
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Persisted nodal temperature is not finite absolute kelvin at node %u.",node);
                nodalTemperature[node] = value;
            }
            this->pModel->convertNodeToElementVector(nodalTemperature,initialElementTemperature);
        }
        else if (temperature.getApplyType() == R_VARIABLE_APPLY_ELEMENT &&
                 temperature.getNValues() == this->pModel->getNElements())
        {
            for (uint element=0;element<this->pModel->getNElements();element++)
            {
                const double value = temperature.getValue(0,element);
                if (!std::isfinite(value) || value < 0.0)
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Persisted element temperature is not finite absolute kelvin at element %u.",element);
                initialElementTemperature[element] = value;
            }
        }
        else
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Persisted temperature has an unsupported scope or sample count.");
        }
    }
    this->sharedData.addData("element-temperature",initialElementTemperature);

    std::vector<RProblemType> problemTypes = RProblem::getTypes(problemTypeMask);
    for (uint i=0;i<problemTypes.size();i++)
    {
        QString problemConvergenceFileName(RFileUtils::getFileNameWithSuffix(this->convergenceFileName,RProblem::getId(problemTypes[i])));
        if (problemTypes[i] == R_PROBLEM_ACOUSTICS)
        {
            this->solvers[problemTypes[i]] = new RSolverAcoustic(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_FLUID_PARTICLE)
        {
            this->solvers[problemTypes[i]] = new RSolverFluidParticle(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_ELECTROSTATICS)
        {
            this->solvers[problemTypes[i]] = new RSolverElectrostatics(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_MAGNETOSTATICS)
        {
            this->solvers[problemTypes[i]] = new RSolverMagnetostatics(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_FLUID)
        {
            this->solvers[problemTypes[i]] = new RSolverFluid(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_FLUID_HEAT)
        {
            this->solvers[problemTypes[i]] = new RSolverFluidHeat(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_RADIATIVE_HEAT)
        {
            this->solvers[problemTypes[i]] = new RSolverRadiativeHeat(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_HEAT)
        {
            this->solvers[problemTypes[i]] = new RSolverHeat(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_STRESS)
        {
            this->solvers[problemTypes[i]] = new RSolverStress(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData,false);
        }
        else if (problemTypes[i] == R_PROBLEM_STRESS_MODAL)
        {
            this->solvers[problemTypes[i]] = new RSolverStress(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData,true);
        }
        else if (problemTypes[i] == R_PROBLEM_POTENTIAL)
        {
            // NOT IMPLEMENTED
        }
        else if (problemTypes[i] == R_PROBLEM_WAVE)
        {
            // NOT WORKING
            this->solvers[problemTypes[i]] = new RSolverWave(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
        else if (problemTypes[i] == R_PROBLEM_MESH)
        {
            this->solvers[problemTypes[i]] = new RSolverMesh(this->pModel,this->modelFileName,problemConvergenceFileName,this->sharedData);
        }
    }

    // if not restarting clear results data.
    if (!this->pModel->getProblemSetup().getRestart())
    {
        std::vector<RVariable> initialKinematics;
        if (this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled())
        {
            // Existing nodal kinematics are an intentional initial state for
            // a new dynamic run (e.g. a selected modal shape). Other results
            // remain stale and must still be cleared.
            for (RVariableType type : {R_VARIABLE_DISPLACEMENT,R_VARIABLE_ROTATION,
                                       R_VARIABLE_VELOCITY,R_VARIABLE_ANGULAR_VELOCITY})
            {
                const uint position = this->pModel->findVariable(type);
                if (position != RConstants::eod)
                {
                    const RVariable &variable = this->pModel->getVariable(position);
                    if (variable.getApplyType() == R_VARIABLE_APPLY_NODE &&
                        variable.getNValues() == this->pModel->getNNodes() &&
                        variable.getNVectors() == 3)
                    {
                        initialKinematics.push_back(variable);
                    }
                }
            }
        }
        this->pModel->removeAllVariables();
        for (const RVariable &variable : initialKinematics)
        {
            this->pModel->addVariable(variable);
        }
    }
}

RSolver::RSolver(RModel &model,
                 const QString &modelFileName,
                 const QString &convergenceFileName,
                 ExecutionPolicy executionPolicy)
    : pModel(&model)
    , modelFileName(modelFileName)
    , convergenceFileName(convergenceFileName)
    , executionPolicy(executionPolicy)
    , transientHeatStressRetryPending(false)
{
    this->_init();
}

RSolver::~RSolver()
{
    foreach (RSolverGeneric *solver, this->solvers)
    {
        delete solver;
    }
}

void RSolver::setParallelExecutionReportObserver(ParallelExecutionReportObserver observer)
{
    this->parallelExecutionReportObserver = std::move(observer);
}

void RSolver::setStructuralSupportResponseObserver(StructuralSupportResponseObserver observer)
{
    if ((this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
         this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled) && observer)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural support-response observers are unsupported in static cable/piezo mode.");
    }
    this->structuralSupportResponseObserver = std::move(observer);
    auto it = this->solvers.find(R_PROBLEM_STRESS);
    if (it != this->solvers.end())
    {
        static_cast<RSolverStress*>(it.value())->setSupportResponseObserver(
            this->structuralSupportResponseObserver);
    }
}

std::vector<RStructuralSupportSample> RSolver::getStructuralSupportSamples() const
{
    if (!this->structuralSupportResponseObserver)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural support samples require an observer registered before run().");
    }
    auto it = this->solvers.constFind(R_PROBLEM_STRESS);
    return it == this->solvers.cend() ? std::vector<RStructuralSupportSample>()
                                     : static_cast<const RSolverStress*>(it.value())->getAcceptedSupportSamples();
}

void RSolver::setStructuralControlCallback(const QString &controllerId,
                                           const QString &controllerRevision,
                                          StructuralControlCallback callback)
{
    const bool cableDynamics =
        this->pModel->getProblemSetup().getStructuralCableSetup().enabled &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    const bool unsupportedSpecializedMode =
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled ||
        (this->pModel->getProblemSetup().getStructuralCableSetup().enabled &&
         !cableDynamics);
    if (unsupportedSpecializedMode && callback)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control callbacks are supported only for transient cable NodalForce v1.");
    }
    this->structuralControllerId = controllerId;
    this->structuralControllerRevision = controllerRevision;
    this->structuralControlCallback = std::move(callback);
    auto it = this->solvers.find(R_PROBLEM_STRESS);
    if (it != this->solvers.end())
        static_cast<RSolverStress*>(it.value())->setControlCallback(
            controllerId,controllerRevision,this->structuralControlCallback);
}

void RSolver::setStructuralControlStepObserver(StructuralControlStepObserver observer)
{
    const bool cableDynamics =
        this->pModel->getProblemSetup().getStructuralCableSetup().enabled &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    const bool unsupportedSpecializedMode =
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled ||
        (this->pModel->getProblemSetup().getStructuralCableSetup().enabled &&
         !cableDynamics);
    if (unsupportedSpecializedMode && observer)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control observers are supported only for transient cable NodalForce v1.");
    }
    this->structuralControlStepObserver = std::move(observer);
    auto it = this->solvers.find(R_PROBLEM_STRESS);
    if (it != this->solvers.end())
        static_cast<RSolverStress*>(it.value())->setControlStepObserver(this->structuralControlStepObserver);
}

std::vector<RStructuralControlChannelSample> RSolver::getStructuralControlSamples() const
{
    if (!this->structuralControlCallback)
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control samples require a registered controller.");
    auto it = this->solvers.constFind(R_PROBLEM_STRESS);
    return it == this->solvers.cend() ? std::vector<RStructuralControlChannelSample>()
                                     : static_cast<const RSolverStress*>(it.value())->getAcceptedControlSamples();
}

void RSolver::setResultFrameObserver(ResultFrameObserver observer)
{
    for (auto it=this->solvers.begin();it!=this->solvers.end();++it)
    {
        it.value()->setResultFrameObserver(observer);
    }
}

void RSolver::run()
{
    const RProblemTaskItem &taskTree = this->pModel->getProblemTaskTree();
    const bool transactionalHeatStress = isHeatStressSequence(taskTree);
    const bool transientHeatStress = transactionalHeatStress &&
        this->pModel->getTimeSolver().getEnabled();
    const bool resumingFailedTransientStep = transientHeatStress &&
        this->transientHeatStressRetryPending;
    if (this->transientHeatStressRetryPending && !transientHeatStress)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "A failed transient Heat->Stress step must be retried with the same task shape and time solver enabled.");
    }
    if (transactionalHeatStress && taskTree.getNIterations() != 1)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transactional Heat->Stress task publication supports one root iteration.");
    }
    const RProblemTypeMask taskMask = taskTree.getProblemTypeMask();
    if (this->pModel->getTimeSolver().getEnabled() &&
        (taskMask & R_PROBLEM_HEAT) && (taskMask & R_PROBLEM_STRESS) &&
        !transactionalHeatStress)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transient Heat->Stress execution requires one root iteration with exactly Heat then Stress leaves.");
    }
    std::unique_ptr<TransientOutputLock> transientOutputLock;
    std::unique_ptr<TransientMonitoringOutputRunLock> transientMonitoringOutputRunLock;
    if (transientHeatStress)
    {
        transientOutputLock.reset(new TransientOutputLock(this->modelFileName));
        transientMonitoringOutputRunLock.reset(new TransientMonitoringOutputRunLock(
            this->pModel->getMonitoringPointManager(),this->modelFileName,
            this->pModel->getTimeSolver()));
    }

    if (this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled)
    {
        // Recheck because the referenced RModel remains mutable after _init().
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
        const bool cableDynamics =
            this->pModel->getProblemSetup().getStructuralCableSetup().enabled &&
            this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
        const bool cableControlEnabled = cableDynamics &&
            this->pModel->getProblemSetup().getStructuralControlSetup().enabled;
        if (this->structuralSupportResponseObserver ||
            ((this->structuralControlCallback || this->structuralControlStepObserver) &&
             !cableControlEnabled))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural support response and control are supported only for transient cable NodalForce v1.");
        }
    }
    RParallelExecutionContext parallelContext(this->pModel->getProblemSetup().getParallelExecutionConfig());
    const RParallelExecutionReport &parallel = parallelContext.getReport();
    if (this->parallelExecutionReportObserver)
    {
        this->parallelExecutionReportObserver(parallel);
    }
    RLogger::info("Parallel execution: requested=%u applied=%u actual=%u backend=%s deterministic=%s\n",
                  parallel.requestedThreads,parallel.appliedThreads,parallel.actualThreads,
                  parallel.backend.toUtf8().constData(),parallel.deterministicReduction ? "yes" : "no");
    if (!parallel.diagnostic.isEmpty())
    {
        RLogger::warning("Parallel execution: %s\n",parallel.diagnostic.toUtf8().constData());
    }
    RTimeSolver &timeSolver = this->pModel->getTimeSolver();
    const RStructuralDynamicsSetup &dynamics = this->pModel->getProblemSetup().getStructuralDynamicsSetup();
    const RStructuralControlSetup &control = this->pModel->getProblemSetup().getStructuralControlSetup();
    control.validate();
    if (control.enabled && (!dynamics.getEnabled() || !this->structuralControlCallback ||
                            this->structuralControllerId != control.controllerId ||
                            this->structuralControllerRevision.isEmpty()))
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Enabled structural control requires dynamics, matching controller ID and nonempty revision fingerprint.");
    if (!control.enabled && this->structuralControlCallback)
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Registered structural controller requires enabled model control setup.");
    if (dynamics.getEnabled())
    {
        dynamics.validate();
        RStructuralDynamicsRestartData &checkpoint =
            this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
        if (this->pModel->getProblemSetup().getRestart())
        {
            checkpoint.validate();
            if (checkpoint.controlAvailable != control.enabled)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Structural restart control enablement differs from its accepted checkpoint.");
            if (!checkpoint.available || checkpoint.meshNodeCount != this->pModel->getNNodes() ||
                checkpoint.meshElementCount != this->pModel->getNElements() ||
                checkpoint.meshFingerprint != RStructuralDynamicsRestartData::computeMeshFingerprint(*this->pModel) ||
                checkpoint.completedSteps != timeSolver.getCurrentTimeStep()+1 ||
                std::fabs(checkpoint.time-timeSolver.getCurrentTime()) >
                    1.0e-10*std::max(1.0,std::fabs(checkpoint.time)) ||
                std::fabs(checkpoint.fixedTimeStep-dynamics.getTimeStep()) >
                    1.0e-12*std::max(checkpoint.fixedTimeStep,dynamics.getTimeStep()))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Structural restart checkpoint does not match the mesh, accepted time step or fixed dt.");
            }
        }
        else
        {
            checkpoint.clear();
            checkpoint.meshNodeCount = this->pModel->getNNodes();
            checkpoint.meshElementCount = this->pModel->getNElements();
            checkpoint.meshFingerprint = RStructuralDynamicsRestartData::computeMeshFingerprint(*this->pModel);
        }
        const RProblemTypeMask problemTypes = this->pModel->getProblemTaskTree().getProblemTypeMask();
        if (problemTypes != R_PROBLEM_STRESS)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Second-order structural dynamics currently requires a stress-only task tree.");
        }
        timeSolver.setEnabled(true);
        timeSolver.setInputNTimeSteps(dynamics.getNSteps());
        timeSolver.setInputTimeStepSize(dynamics.getTimeStep());
        timeSolver.setOutputFrequency(dynamics.getOutputFrequency());
    }
    if (!resumingFailedTransientStep)
    {
        timeSolver.harmonizeTimesWithInput(this->pModel->getProblemSetup().getRestart());
    }
    if (transientHeatStress && !resumingFailedTransientStep)
    {
        // There is no safe point at which to delete old future records if the
        // run may stop before reaching them. Reject stale in-schedule records
        // before the first trial instead of letting them masquerade as output.
        rejectUnpublishedOutputRecords(this->modelFileName,timeSolver,
                                       this->pModel->getProblemSetup().getRestart());
    }

    // A harmonic acoustic analysis is a frequency domain analysis - the solver
    // sweeps frequencies internally and must not be driven by the time loop.
    bool harmonicAcoustics = (this->pModel->getProblemTaskTree().getProblemTypeMask() & R_PROBLEM_ACOUSTICS) &&
                             (this->pModel->getProblemSetup().getAcousticSetup().getAnalysisType() == R_ACOUSTIC_ANALYSIS_HARMONIC);

    if (harmonicAcoustics)
    {
        if (timeSolver.getEnabled())
        {
            RLogger::warning("Harmonic acoustic analysis is solved in the frequency domain - time stepping is skipped.\n");
        }
        // The time solver setup is left untouched so that switching back to a
        // transient analysis does not lose it.
        this->runSingle();
    }
    else if (timeSolver.getEnabled() && RProblem::getTimeSolverEnabled(this->pModel->getProblemTaskTree().getProblemTypeMask()))
    {
        if (this->pModel->getProblemSetup().getRestart() && !resumingFailedTransientStep)
        {
            timeSolver.setNextTimeStep();
        }
        // The legacy pre-pass rewrites prior records and removes future ones
        // before a trial step is accepted. Structural dynamics publishes its
        // checkpoint/results transactionally; preserve existing files on a
        // failed step and let RModel::write publish the new latest link only
        // after its target record succeeds.
        if (!dynamics.getEnabled() && !transientHeatStress)
        {
            RSolverGeneric::updateOldRecords(timeSolver,this->modelFileName);
        }

        do {
            RLogger::info("time step: %9u of %-9u | time = % 12e [sec] | dt = % 12e [sec]\n",
                          timeSolver.getCurrentTimeStep()+1,
                          timeSolver.getNTimeSteps(),
                          timeSolver.getCurrentTime(),
                          timeSolver.getCurrentTimeStepSize());
            const double acceptedComputedTime = timeSolver.getComputedTime();
            timeSolver.setComputedTime(timeSolver.getCurrentTime());

            RLogger::indent();
            try
            {
                this->runSingle();
            }
            catch (...)
            {
                if (transientHeatStress)
                {
                    // Keep the failed index selected so the same trial step can
                    // be retried, while computedTime continues to identify the
                    // last accepted coupled step.
                    timeSolver.setComputedTime(acceptedComputedTime);
                    this->transientHeatStressRetryPending = true;
                }
                RLogger::unindent();
                throw;
            }
            RLogger::unindent();
            if (transientHeatStress)
            {
                this->transientHeatStressRetryPending = false;
            }
            if (control.enabled)
            {
                const auto *stress = static_cast<const RSolverStress*>(this->solvers.value(R_PROBLEM_STRESS));
                const QString failure = stress->getPendingControlObserverFailure();
                if (!failure.isEmpty())
                    throw RError(RError::Type::Application,R_ERROR_REF,
                                 "Control audit observer failed after committed time step: %s",
                                 failure.toUtf8().constData());
            }

            if (RApplicationState::getInstance().getStateType() == R_APPLICATION_STATE_STOP)
            {
                break;
            }
        } while (timeSolver.setNextTimeStep() != RConstants::eod);
    }
    else
    {
        timeSolver.setCurrentTimeStep(0);
        timeSolver.setComputedTime(0.0);
        timeSolver.setTimes(RTimeSolver::findTimesVector(1, 0.0, 0.0));
        this->runSingle();
    }
}

void RSolver::runSingle()
{
    const RProblemTaskItem &taskTree = this->pModel->getProblemTaskTree();
    if (!isHeatStressSequence(taskTree))
    {
        this->runProblemTask(taskTree,0);
        return;
    }

    RSolverGeneric *heatSolver = this->solvers.value(R_PROBLEM_HEAT,nullptr);
    RSolverGeneric *stressSolver = this->solvers.value(R_PROBLEM_STRESS,nullptr);
    if (!heatSolver || !stressSolver)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Transactional Heat->Stress task is missing a registered sub-solver.");
    }
    if (taskTree.getNIterations() != 1)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transactional Heat->Stress task publication supports one root iteration.");
    }

    const RModel acceptedModel(*this->pModel);
    const RSolverSharedData acceptedSharedData(this->sharedData);
    const QMap<RProblemType,uint> acceptedExecutionCount(this->solversExecutionCount);
    DeferredOutputTransaction outputTransaction;
    TransientMonitoringOutputTransaction monitoringOutputTransaction;
    const bool transient = this->pModel->getTimeSolver().getEnabled();
    heatSolver->setResultWritesDeferred(true);
    stressSolver->setResultWritesDeferred(true);
    heatSolver->setMonitoringOutputFailuresAreFatal(transient);
    stressSolver->setMonitoringOutputFailuresAreFatal(transient);
    try
    {
        outputTransaction.prepare(this->modelFileName,this->pModel->getTimeSolver());
        if (transient)
        {
            monitoringOutputTransaction.prepare(this->pModel->getMonitoringPointManager(),
                                                this->modelFileName,
                                                this->pModel->getTimeSolver());
        }
        if (!this->runProblemTask(taskTree,0))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Sequential Heat->Stress task did not converge; transaction rolled back.");
        }

        // Both sub-solvers have populated the same in-memory model. Publish it
        // exactly once, after the enclosing task reports convergence.
        monitoringOutputTransaction.restoreOutputFileNames();
        stressSolver->publishDeferredResults();
        monitoringOutputTransaction.publish();
        outputTransaction.commit();
        monitoringOutputTransaction.commit();
        heatSolver->setResultWritesDeferred(false);
        stressSolver->setResultWritesDeferred(false);
        heatSolver->setMonitoringOutputFailuresAreFatal(false);
        stressSolver->setMonitoringOutputFailuresAreFatal(false);
        stressSolver->publishResultFrameObserver();
    }
    catch (...)
    {
        const std::exception_ptr originalFailure = std::current_exception();
        heatSolver->discardPendingResultFrameObserver();
        stressSolver->discardPendingResultFrameObserver();
        const QString outputRollbackFailure = outputTransaction.rollback();
        const QString monitoringRollbackFailure = monitoringOutputTransaction.rollback();
        heatSolver->setResultWritesDeferred(false);
        stressSolver->setResultWritesDeferred(false);
        heatSolver->setMonitoringOutputFailuresAreFatal(false);
        stressSolver->setMonitoringOutputFailuresAreFatal(false);
        bool modelRollbackFailed = false;
        try
        {
            *this->pModel = acceptedModel;
            this->sharedData = acceptedSharedData;
            this->solversExecutionCount = acceptedExecutionCount;
        }
        catch (...)
        {
            modelRollbackFailed = true;
        }
        if (!outputRollbackFailure.isEmpty() || !monitoringRollbackFailure.isEmpty() || modelRollbackFailed)
        {
            const QString modelDiagnostic = modelRollbackFailed
                ? QStringLiteral(" Model/shared-state restoration also failed.") : QString();
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Heat->Stress transaction failed and rollback was incomplete. %s%s",
                         (outputRollbackFailure+monitoringRollbackFailure).toUtf8().constData(),
                         modelDiagnostic.toUtf8().constData());
        }
        std::rethrow_exception(originalFailure);
    }
}

bool RSolver::runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration)
{
    bool converged = false;

    if (problemTaskItem.getProblemType() == R_PROBLEM_NONE)
    {
        for (uint i=0;i<problemTaskItem.getNIterations();i++)
        {
            uint nConverged = 0;
            RLogger::info("Problem task iteration: %u of %u\n",i+1,problemTaskItem.getNIterations());
            RLogger::indent();
            for (uint j=0;j<problemTaskItem.getNChildren();j++)
            {
                if (this->runProblemTask(problemTaskItem.getChild(j),i))
                {
                    nConverged++;
                }
            }
            if (nConverged == problemTaskItem.getNChildren())
            {
                converged = true;
            }
            RLogger::unindent();
            if (converged)
            {
                RLogger::info("All sub-tasks have converged.\n");
                break;
            }

            if (RApplicationState::getInstance().getStateType() == R_APPLICATION_STATE_STOP)
            {
                break;
            }
        }
    }
    else
    {
        if (this->solvers.contains(problemTaskItem.getProblemType()))
        {
            RProblemType problemType = problemTaskItem.getProblemType();
            if (!this->solversExecutionCount.contains(problemType))
            {
                this->solversExecutionCount[problemType] = 0;
            }
            uint executionCount = this->solversExecutionCount[problemType];

            bool firstExecution = (executionCount == 0);

            RLogger::info("Solving problem task: %s\n",RProblem::getName(problemTaskItem.getProblemType()).toUtf8().constData());
            RLogger::indent();

            this->solvers[problemType]->run(firstExecution,taskIteration);
            converged = this->solvers[problemTaskItem.getProblemType()]->hasConverged();
            if (this->solvers[problemTaskItem.getProblemType()]->getMeshChanged())
            {
                RProblemTypeMask problemTypeMask = this->pModel->getProblemTaskTree().getProblemTypeMask();
                std::vector<RProblemType> problemTypes = RProblem::getTypes(problemTypeMask);
                for (uint i=0;i<problemTypes.size();i++)
                {
                    if (!this->solvers.contains(problemTypes[i]))
                    {
                        throw RError(RError::Type::Application,R_ERROR_REF,
                                     "No solver is registered for enabled problem type '%s'.",
                                     RProblem::getId(problemTypes[i]).toUtf8().constData());
                    }
                    this->solvers[problemTypes[i]]->setMeshChanged(true);
                }
            }

            RLogger::unindent();

            this->solversExecutionCount[problemType]++;
        }
        else
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "No solver is registered for enabled problem type '%s'.",
                         RProblem::getId(problemTaskItem.getProblemType()).toUtf8().constData());
        }
    }

    return converged;
}
