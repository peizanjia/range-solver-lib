#include <QFile>
#include <QFileInfo>
#include <cstdio>
#include <exception>
#include <memory>
#include <omp.h>

#include <rbl_application_state.h>

#include <rml_file_utils.h>

#include "rsolvergeneric.h"


const double RSolverGeneric::sigma = 5.670367e-8; // Stefan-Boltzman constant [W/(m*m*K*K*K*K)]
const double RSolverGeneric::e0 = 8.854187817e-12; // Vacuum permittivity [F/m]

RSolverGeneric::RSolverGeneric(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : meshChanged(true)
    , problemType(R_PROBLEM_NONE)
    , pModel(pModel)
    , modelFileName(modelFileName)
    , convergenceFileName(convergenceFileName)
    , pSharedData(&sharedData)
    , firstRun(false)
    , resultWritesDeferred(false)
    , monitoringOutputFailuresAreFatal(false)
    , taskIteration(0)
    , computableElements(this->pModel->getNElements(),false)
{
    this->elementTemperature = this->pSharedData->findData("element-temperature");
    this->elementTemperature.resize(this->pModel->getNElements(),RVariable::getInitValue(R_VARIABLE_TEMPERATURE));
    this->localRotations.resize(this->pModel->getNNodes());
}

RSolverGeneric::~RSolverGeneric()
{

}

void RSolverGeneric::beginStepTransaction() {}
void RSolverGeneric::rollbackStepTransaction() {}
void RSolverGeneric::commitStepTransaction() {}
void RSolverGeneric::validateRunConfiguration() const {}

void RSolverGeneric::run(bool firstExecution, uint taskIteration)
{
    this->validateRunConfiguration();
    this->firstRun = firstExecution && !this->pModel->getProblemSetup().getRestart();
    this->taskIteration = taskIteration;

    this->elementTemperature.resize(this->pModel->getNElements(),RVariable::getInitValue(R_VARIABLE_TEMPERATURE));

    if (this->taskIteration == 0)
    {
        this->findComputableElements(this->problemType);
        this->findIncludableElements();
        this->findInwardElements();
        this->updateLocalRotations();
    }

    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        RModalSetup &modalSetup = this->pModel->getProblemSetup().getModalSetup();

        this->updateScales();
        this->scales.downscale(*this->pModel);

        this->initialize();
        this->recoverSharedData();
        this->recover();
        this->prepare();
        this->solve();

        this->scales.upscale(*this->pModel);

        uint nModes = modalSetup.getNModesToExtract();
        if (modalSetup.getMethod() == R_MODAL_DOMINANT_MODE)
        {
            nModes = 1;
        }

        // Only the modes which were actually produced can be stored. The eigen
        // value solver can return fewer than were asked for.
        uint nComputedModes = this->getNComputedModes();
        if (nComputedModes < nModes)
        {
            RLogger::warning("Only %u of the %u requested modes were computed - the remaining modes are not stored.\n",nComputedModes,nModes);
            nModes = nComputedModes;
        }

        for (uint i=0;i<nModes;i++)
        {
            uint mode = nModes - (i + 1);
            RLogger::info("Extracting mode %d\n",mode+1);
            modalSetup.setMode(mode);

            RLogger::indent();

            this->scales.downscale(*this->pModel);

            this->process();
            this->store();
            this->storeSharedData();

            this->scales.upscale(*this->pModel);

            this->writeResults();
            this->statistics();

            RLogger::unindent();

            if (RApplicationState::getInstance().getStateType() == R_APPLICATION_STATE_STOP)
            {
                break;
            }
        }
    }
    else if (this->problemType == R_PROBLEM_ACOUSTICS &&
             this->pModel->getProblemSetup().getAcousticSetup().getAnalysisType() == R_ACOUSTIC_ANALYSIS_HARMONIC)
    {
        // Frequency-domain acoustics - sweep the requested frequencies, each one
        // producing its own result record (see RModel::write).
        RAcousticSetup &acousticSetup = this->pModel->getProblemSetup().getAcousticSetup();

        this->applyDisplacement();

        if (firstExecution)
        {
            this->updateScales();
        }

        uint nFrequencies = std::max(acousticSetup.getNFrequencies(),uint(1));

        for (uint i=0;i<nFrequencies;i++)
        {
            acousticSetup.setFrequencyIndex(i);

            RLogger::info("Frequency: %9u of %-9u | f = % 12e [Hz]\n",
                          i+1,
                          nFrequencies,
                          acousticSetup.getFrequency());
            RLogger::indent();

            this->scales.downscale(*this->pModel);

            this->initialize();
            this->recoverSharedData();
            this->recover();
            this->prepare();
            this->solve();
            this->process();
            this->store();
            this->storeSharedData();

            this->scales.upscale(*this->pModel);

            this->writeResults();
            this->statistics();

            RLogger::unindent();

            if (RApplicationState::getInstance().getStateType() == R_APPLICATION_STATE_STOP)
            {
                break;
            }
        }

        this->removeDisplacement();
    }
    else
    {
        const bool structuralDynamics = this->problemType == R_PROBLEM_STRESS &&
            this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
        // The time driver has already selected the trial time. Preserve the
        // previously accepted physical fields and shared data before any
        // solve/store mutation. Copying is deliberately conservative here;
        // a result/setup journal should replace this O(model size) per-step
        // copy after the rollback contract is established.
        std::unique_ptr<RModel> acceptedModel;
        std::unique_ptr<RSolverSharedData> acceptedSharedData;
        if (structuralDynamics)
        {
            acceptedModel.reset(new RModel(*this->pModel));
            acceptedSharedData.reset(new RSolverSharedData(*this->pSharedData));
        }
        bool modelScaled = false;
        bool transactionStarted = false;
        try
        {
            if (structuralDynamics)
            {
                this->beginStepTransaction();
                transactionStarted = true;
            }
            if (this->problemType != R_PROBLEM_STRESS && this->problemType != R_PROBLEM_STRESS_MODAL && this->problemType != R_PROBLEM_MESH)
            {
                this->applyDisplacement();
            }

            if (firstExecution || this->problemType == R_PROBLEM_MESH)
            {
                this->updateScales();
            }
            if (this->problemType != R_PROBLEM_MESH)
            {
                this->scales.downscale(*this->pModel);
                modelScaled = true;
            }

            this->initialize();
            this->recoverSharedData();
            this->recover();
            this->prepare();
            this->solve();
            this->process();
            this->store();
            this->storeSharedData();

            if (modelScaled)
            {
                this->scales.upscale(*this->pModel);
                modelScaled = false;
            }

            if (this->problemType != R_PROBLEM_STRESS && this->problemType != R_PROBLEM_STRESS_MODAL)
            {
                this->removeDisplacement();
            }

            this->statistics();
            this->writeResults();
            if (transactionStarted)
            {
                this->commitStepTransaction();
            }
        }
        catch (...)
        {
            if (structuralDynamics)
            {
                const std::exception_ptr originalFailure = std::current_exception();
                try
                {
                    // A late write failure is as much a rejected step as an
                    // assembly/solve failure. Restore all in-memory state,
                    // then return to the last accepted checkpoint time.
                    *this->pModel = *acceptedModel;
                    *this->pSharedData = *acceptedSharedData;
                    RTimeSolver &time = this->pModel->getTimeSolver();
                    const auto &accepted = this->pModel->getProblemSetup()
                        .getStructuralDynamicsRestartData();
                    if (accepted.available && accepted.completedSteps > 0)
                    {
                        time.setCurrentTimeStep(accepted.completedSteps-1);
                        time.setComputedTime(accepted.time);
                    }
                    else
                    {
                        time.setCurrentTimeStep(0);
                        time.setComputedTime(time.getPreviousTime());
                    }
                    if (transactionStarted)
                    {
                        this->rollbackStepTransaction();
                    }
                }
                catch (const RError &rollbackError)
                {
                    RLogger::warning("Structural step rollback also failed: %s\n",
                                     rollbackError.getMessage().toUtf8().constData());
                }
                catch (const std::exception &rollbackError)
                {
                    RLogger::warning("Structural step rollback also failed: %s\n",rollbackError.what());
                }
                catch (...)
                {
                    RLogger::warning("Structural step rollback also failed with an unknown error.\n");
                }
                std::rethrow_exception(originalFailure);
            }
            else if (modelScaled)
            {
                // Preserve the existing exception-safe scale restoration for
                // all non-transactional solver paths.
                this->scales.upscale(*this->pModel);
            }
            throw;
        }
    }

    this->meshChanged = (this->problemType == R_PROBLEM_MESH);
}

uint RSolverGeneric::getNComputedModes() const
{
    return 0;
}

bool RSolverGeneric::getMeshChanged() const
{
    return this->meshChanged;
}

void RSolverGeneric::setMeshChanged(bool meshChanged)
{
    this->meshChanged = meshChanged;
}

void RSolverGeneric::updateOldRecords(const RTimeSolver &rTimeSolver, const QString &modelFileName)
{
    if (rTimeSolver.getEnabled())
    {
        // Update previous records.
        for (unsigned int i=0;i<=rTimeSolver.getCurrentTimeStep();i++)
        {
            QString fileName = RFileUtils::getFileNameWithTimeStep(modelFileName,i);

            if (RFileUtils::fileExists(fileName))
            {
                RLogger::info("Updating model file \'%s\'\n",fileName.toUtf8().constData());

                RModel model;
                model.read(fileName);
                model.getTimeSolver().setTimes(rTimeSolver.getTimes());
                model.getTimeSolver().setInputNTimeSteps(rTimeSolver.getInputNTimeSteps());
                model.getTimeSolver().setInputStartTime(rTimeSolver.getInputStartTime());
                model.getTimeSolver().setInputTimeStepSize(rTimeSolver.getInputTimeStepSize());
                model.write(fileName);
            }
        }
        // Delete records which will be computed again.
        for (unsigned int i=rTimeSolver.getCurrentTimeStep();i<rTimeSolver.getNTimeSteps();i++)
        {
            QString fileName = RFileUtils::getFileNameWithTimeStep(modelFileName,i+1);

            if (RFileUtils::fileExists(fileName))
            {
                RLogger::info("Removing model file \'%s\'\n",fileName.toUtf8().constData());
                try
                {
                    RFileUtils::remove(fileName);
                }
                catch (RError &error)
                {
                    RLogger::warning("Failed to remove model file \'%s\'. ERROR: %s\n",fileName.toUtf8().constData(),error.getMessage().toUtf8().constData());
                }
            }
        }
    }
}

RRVector RSolverGeneric::findElementSizes() const
{
    uint ne = this->pModel->getNElements();
    RRVector elementSizes(ne,0.0);
    #pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(ne);i++)
    {
        elementSizes[uint(i)] = this->pModel->getElement(uint(i)).findSize(this->pModel->getNodes());
    }
    return elementSizes;
}

double RSolverGeneric::findElementScale(bool onlyVolumes) const
{
    RRVector elementSizes = this->findElementSizes();

    if (onlyVolumes)
    {
        uint ne = this->pModel->getNElements();

        std::vector<bool> volumeElements;
        volumeElements.resize(ne,false);

        for (uint i=0;i<this->pModel->getNLines();i++)
        {
            const RLine &rLine = this->pModel->getLine(i);
            if (rLine.getCrossArea() > RConstants::eps)
            {
                for (uint j=0;j<rLine.size();j++)
                {
                    volumeElements[rLine.get(j)] = true;
                }
            }
        }
        for (uint i=0;i<this->pModel->getNSurfaces();i++)
        {
            const RSurface &rSurface = this->pModel->getSurface(i);
            if (rSurface.getThickness() > RConstants::eps)
            {
                for (uint j=0;j<rSurface.size();j++)
                {
                    volumeElements[rSurface.get(j)] = true;
                }
            }
        }
        for (uint i=0;i<this->pModel->getNVolumes();i++)
        {
            const RVolume &rVolume = this->pModel->getVolume(i);
            for (uint j=0;j<rVolume.size();j++)
            {
                volumeElements[rVolume.get(j)] = true;
            }
        }

        uint nve = 0;
        for (uint i=0;i<ne;i++)
        {
            if (volumeElements[i])
            {
                elementSizes[nve] = elementSizes[i];
                nve ++;
            }
        }
        elementSizes.resize(nve);
    }
    double elementSize = RStatistics::findMedianValue(elementSizes);

    if (elementSize < RConstants::eps)
    {
        elementSize = 1.0;
    }
    return 1.0 / elementSize;
}

double RSolverGeneric::findMeshScale() const
{
    return this->pModel->findNodeScale();
}

void RSolverGeneric::updateLocalRotations()
{
    std::vector<RR3Vector> nodeNormals(this->pModel->getNNodes());
    RBVector nodeNormalsBook(this->pModel->getNNodes(),false);
    RR3Vector normal;
    RRMatrix R;

    // Calculate directions for surface elements
    for (unsigned int i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);
        bool findLocalDirections = false;
        bool useExplicitDirection = false;
        RR3Vector explicitDirection;

        for (unsigned int j=0;j<rSurface.getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = rSurface.getBoundaryCondition(j);
            if (RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType)
            {
                if (bc.getHasLocalDirection())
                {
                    findLocalDirections = true;
                    // An explicitly entered direction replaces the one which
                    // would otherwise be averaged from the element normals.
                    useExplicitDirection = bc.getExplicitLocalDirection();
                    explicitDirection = bc.getLocalDirection();
                    break;
                }
            }
        }
        if (useExplicitDirection && explicitDirection.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of surface entity \'%s\' has zero length.",
                         rSurface.getName().toUtf8().constData());
        }
        for (unsigned int j=0;j<rSurface.size();j++)
        {
            if (findLocalDirections)
            {
                const RElement &rElement = this->pModel->getElement(rSurface.get(j));

                if (useExplicitDirection)
                {
                    normal = explicitDirection;
                }
                else if (!rElement.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]))
                {
                    throw RError(RError::Type::Application,R_ERROR_REF,"Could not calculate element normal for element# = %u.",rSurface.get(j));
                }
                for (unsigned int k=0;k<rElement.size();k++)
                {
                    uint nodeId = rElement.getNodeId(k);
                    nodeNormals[nodeId][0] += normal[0];
                    nodeNormals[nodeId][1] += normal[1];
                    nodeNormals[nodeId][2] += normal[2];
                    nodeNormalsBook[nodeId] = true;
                }
            }
        }
    }

    // Apply directions for point elements
    for (unsigned int i=0;i<this->pModel->getNPoints();i++)
    {
        RR3Vector direction;
        bool hasLocalDirections = false;
        const RPoint &rPoint = this->pModel->getPoint(i);
        for (unsigned int j=0;j<rPoint.getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = rPoint.getBoundaryCondition(j);
            if (RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType)
            {
                if (bc.getHasLocalDirection())
                {
                    hasLocalDirections = true;
                    direction = bc.getLocalDirection();
                    break;
                }
            }
        }
        if (hasLocalDirections)
        {
            for (uint j=0;j<rPoint.size();j++)
            {
                const RElement &rElement = this->pModel->getElement(rPoint.get(j));
                for (unsigned int k=0;k<rElement.size();k++)
                {
                    uint nodeId = rElement.getNodeId(k);
                    nodeNormals[nodeId][0] = direction[0];
                    nodeNormals[nodeId][1] = direction[1];
                    nodeNormals[nodeId][2] = direction[2];
                    nodeNormalsBook[nodeId] = true;
                }
            }
        }
    }

    // Convert local directions to rotation matricies
    this->localRotations.resize(this->pModel->getNNodes());
    for (uint i=0;i<nodeNormals.size();i++)
    {
        if (nodeNormalsBook[i])
        {
            nodeNormals[i].normalize();
            nodeNormals[i].findRotationMatrix(R);
            this->localRotations[i].activate(R);
        }
        else
        {
            this->localRotations[i].deactivate();
        }
    }

    // Apply local rotation for line elements
    for (unsigned int i=0;i<this->pModel->getNLines();i++)
    {
        RR3Vector direction;
        bool findLocalDirections = false;
        bool useExplicitDirection = false;
        const RLine &rLine = this->pModel->getLine(i);
        for (unsigned int j=0;j<rLine.getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = rLine.getBoundaryCondition(j);
            if (RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType)
            {
                if (bc.getHasLocalDirection())
                {
                    findLocalDirections = true;
                    // An explicitly entered direction replaces the frame which
                    // would otherwise follow the element direction.
                    useExplicitDirection = bc.getExplicitLocalDirection();
                    direction = bc.getLocalDirection();
                    break;
                }
            }
        }
        if (findLocalDirections)
        {
            if (useExplicitDirection)
            {
                RR3Vector localDirection(direction);
                if (localDirection.normalize() == 0.0)
                {
                    throw RError(RError::Type::Application,R_ERROR_REF,
                                 "Local direction of line entity \'%s\' has zero length.",
                                 rLine.getName().toUtf8().constData());
                }
                localDirection.findRotationMatrix(R);

                for (uint j=0;j<rLine.size();j++)
                {
                    const RElement &rElement = this->pModel->getElement(rLine.get(j));
                    for (uint k=0;k<rElement.size();k++)
                    {
                        this->localRotations[rElement.getNodeId(k)].activate(R);
                    }
                }
                continue;
            }

            for (uint j=0;j<rLine.size();j++)
            {
                const RElement &rElement = this->pModel->getElement(rLine.get(j));

                RR3Vector d1,d2,d3;
                RSegment s(this->pModel->getNode(rElement.getNodeId(0)),
                           this->pModel->getNode(rElement.getNodeId(1)));
                s.findPerpendicularVectors(d1,d2,d3);

                R.resize(3,3,0.0);
                R[0][0]=d1[0]; R[0][1]=d1[1]; R[0][2]=d1[2];
                R[1][0]=d2[0]; R[1][1]=d2[1]; R[1][2]=d2[2];
                R[2][0]=d3[0]; R[2][1]=d3[1]; R[2][2]=d3[2];
                R.invert();

                for (uint k=0;k<rElement.size();k++)
                {
                    this->localRotations[rElement.getNodeId(k)].activate(R);
                }
            }
        }
    }
}

void RSolverGeneric::clearSharedData()
{
    this->pSharedData->clearData();
}

void RSolverGeneric::storeSharedData()
{
    this->pSharedData->addData("element-temperature",this->elementTemperature);
}

void RSolverGeneric::recoverSharedData()
{
    if (this->pSharedData->hasData("element-temperature",this->pModel->getNElements()))
    {
        this->elementTemperature = this->pSharedData->findData("element-temperature");
    }
}

void RSolverGeneric::writeResults()
{
    if (this->resultWritesDeferred || this->modelFileName.isEmpty())
    {
        return;
    }

    // A harmonic acoustic analysis writes one record per frequency and is not
    // subject to the time solver output frequency.
    bool harmonicAcoustics = (this->problemType == R_PROBLEM_ACOUSTICS) &&
                             (this->pModel->getProblemSetup().getAcousticSetup().getAnalysisType() == R_ACOUSTIC_ANALYSIS_HARMONIC);

    bool canWrite = true;
    if (!harmonicAcoustics && this->pModel->getTimeSolver().getEnabled())
    {
        canWrite = false;
        uint outputFrequency = this->pModel->getTimeSolver().getOutputFrequency();
        if ((this->pModel->getTimeSolver().getCurrentTimeStep()+1) == this->pModel->getTimeSolver().getNTimeSteps() ||
            (outputFrequency > 0 && (this->pModel->getTimeSolver().getCurrentTimeStep()+1) % outputFrequency == 0))
        {
            canWrite = true;
        }
    }

    if (canWrite)
    {
        this->modelFileName = this->pModel->write(this->modelFileName);
    }
}

void RSolverGeneric::setResultWritesDeferred(bool deferred)
{
    this->resultWritesDeferred = deferred;
}

void RSolverGeneric::setMonitoringOutputFailuresAreFatal(bool fatal)
{
    this->monitoringOutputFailuresAreFatal = fatal;
}

void RSolverGeneric::publishDeferredResults()
{
    const bool wasDeferred = this->resultWritesDeferred;
    this->resultWritesDeferred = false;
    try
    {
        this->writeResults();
    }
    catch (...)
    {
        this->resultWritesDeferred = wasDeferred;
        throw;
    }
    this->resultWritesDeferred = wasDeferred;
}

void RSolverGeneric::applyDisplacement()
{
    uint variablePosition = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
    if (variablePosition == RConstants::eod)
    {
        return;
    }

    const RVariable &rVariable = this->pModel->getVariable(variablePosition);

    R_ERROR_ASSERT(rVariable.getNValues() == this->pModel->getNNodes());
    R_ERROR_ASSERT(rVariable.getNVectors() == 3);

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        RRVector u = rVariable.getValueVector(i);
        this->pModel->getNode(i).move(RR3Vector(u[0],u[1],u[2]));
    }
}

void RSolverGeneric::removeDisplacement()
{
    uint variablePosition = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
    if (variablePosition == RConstants::eod)
    {
        return;
    }

    const RVariable &rVariable = this->pModel->getVariable(variablePosition);

    R_ERROR_ASSERT(rVariable.getNValues() == this->pModel->getNNodes());
    R_ERROR_ASSERT(rVariable.getNVectors() == 3);

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        RRVector u = rVariable.getValueVector(i);
        u *= -1.0;
        this->pModel->getNode(i).move(RR3Vector(u[0],u[1],u[2]));
    }
}

void RSolverGeneric::generateNodeBook(RProblemType problemType)
{
    if (problemType == R_PROBLEM_FLUID)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Unsupported problem type \'%s\' for node book generation.",
                     RProblem::getName(problemType).toUtf8().constData());
    }

    this->nodeBook.resize(this->pModel->getNNodes());
    RBVector disabledPositions(this->nodeBook.size(),false);

    for (unsigned int i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        bool hasExplicitBc = false;
        for (unsigned int j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType)
            {
                if (bc.getExplicit())
                {
                    hasExplicitBc = true;
                    break;
                }
            }
        }
        for (unsigned int j=0;j<pElementGroup->size();j++)
        {
            if (hasExplicitBc)
            {
                const RElement &element = this->pModel->getElement(pElementGroup->get(j));
                for (unsigned int k=0;k<element.size();k++)
                {
                    disabledPositions[element.getNodeId(k)] = true;
                }
            }
        }
    }
    RBVector computableNodes(this->pModel->getNNodes(),false);
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (this->computableElements[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            for (uint j=0;j<rElement.size();j++)
            {
                computableNodes[rElement.getNodeId(j)] = true;
            }
        }
    }
    for (uint i=0;i<computableNodes.size();i++)
    {
        if (!computableNodes[i])
        {
            disabledPositions[i] = true;
        }
    }

    RSolverGeneric::rebuildNodeBook(this->nodeBook,disabledPositions);
}

void RSolverGeneric::rebuildNodeBook(RBook &nodeBook, const RBVector &disabledPositions)
{
    R_ERROR_ASSERT(nodeBook.size() == disabledPositions.size());

    uint nextPosition = 0;
    for (uint i=0;i<nodeBook.size();i++)
    {
        if (disabledPositions[i])
        {
            nodeBook.setValue(i,RConstants::eod);
        }
        else
        {
            nodeBook.setValue(i,nextPosition);
            nextPosition++;
        }
    }
}

void RSolverGeneric::generateMaterialVecor(RMaterialProperty::Type materialPropertyType, RRVector &materialPropertyValues) const
{
    unsigned int ne = this->pModel->getNElements();

    materialPropertyValues.resize(ne,0.0);

    for (unsigned int i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        const RMaterial &material = pElementGroup->getMaterial();

        // Get material property
        unsigned int materialPropertyPosition = material.findPosition(materialPropertyType);
        if (materialPropertyPosition == material.size())
        {
//            RLogger::warning("Could not find material property \'%s\' in element group \'%s\'\n",
//                             RMaterialProperty::getName(materialPropertyType).toUtf8().constData(),
//                             pElementGroup->getName().toUtf8().constData());
            continue;
        }
        // Apply values to element vectors.
        for (unsigned int j=0;j<pElementGroup->size();j++)
        {
            unsigned int elementID = pElementGroup->get(j);
            double temperature = this->elementTemperature[elementID];
            RMaterialEvaluationContext context;
            context.setField(R_VARIABLE_TEMPERATURE,temperature);
            materialPropertyValues[elementID] = material.evaluateProperty(materialPropertyType,context).getValue();
        }
    }
}

void RSolverGeneric::generateVariableVector(RVariableType variableType,
                                            RRVector &variableValues,
                                            RBVector &setValues,
                                            bool boundaryConditions,
                                            bool initialConditions,
                                            bool environmentConditions,
                                            bool onlyExplicitBcs) const
{
    unsigned int ne = this->pModel->getNElements();

    variableValues.resize(ne,0.0);
    setValues.resize(ne,false);

    for (unsigned int i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        if (environmentConditions)
        {
            // Apply environment conditions.
            for (unsigned int j=0;j<pElementGroup->getNEnvironmentConditions();j++)
            {
                const REnvironmentCondition &ec = pElementGroup->getEnvironmentCondition(j);
                unsigned int ecComponentPosition = ec.findComponentPosition(variableType);
                if (ecComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = ec.getComponent(ecComponentPosition);
                if (!conditionComponent.getEnabled())
                {
                    continue;
                }
                double value = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());

                for (unsigned int k=0;k<pElementGroup->size();k++)
                {
                    variableValues[pElementGroup->get(k)] = value;
                    setValues[pElementGroup->get(k)] = true;
                }
            }
        }

        if (initialConditions)
        {
            // Apply initial conditions.
            for (unsigned int j=0;j<pElementGroup->getNInitialConditions();j++)
            {
                const RInitialCondition &ic = pElementGroup->getInitialCondition(j);
                unsigned int icComponentPosition = ic.findComponentPosition(variableType);
                if (icComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = ic.getComponent(icComponentPosition);
                if (!conditionComponent.getEnabled())
                {
                    continue;
                }
                double value = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());

                for (unsigned int k=0;k<pElementGroup->size();k++)
                {
                    variableValues[pElementGroup->get(k)] = value;
                    setValues[pElementGroup->get(k)] = true;
                }
            }
        }

        if (boundaryConditions)
        {
            // Apply boundary conditions.
            for (unsigned int j=0;j<pElementGroup->getNBoundaryConditions();j++)
            {
                const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
                if (onlyExplicitBcs && !bc.getExplicit())
                {
                    continue;
                }
                unsigned int bcComponentPosition = bc.findComponentPosition(variableType);
                if (bcComponentPosition == RConstants::eod)
                {
                    continue;
                }
                const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
                // A component of an optional condition can be switched off - it
                // then prescribes nothing at all.
                if (!conditionComponent.getEnabled())
                {
                    continue;
                }
                double value = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());

                for (unsigned int k=0;k<pElementGroup->size();k++)
                {
                    variableValues[pElementGroup->get(k)] = value;
                    setValues[pElementGroup->get(k)] = true;
                }
            }
        }
    }
}

void RSolverGeneric::generateHeatVector(RRVector &heatValues, RBVector &setValues) const
{
    unsigned int ne = this->pModel->getNElements();

    heatValues.resize(ne);
    heatValues.fill(0.0);
    setValues.resize(ne);
    setValues.fill(false);

    for (unsigned int i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        for (unsigned int j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            unsigned int bcComponentPosition = bc.findComponentPosition(R_VARIABLE_HEAT);
            if (bcComponentPosition == RConstants::eod)
            {
                continue;
            }
            const RConditionComponent &conditionComponent = bc.getComponent(bcComponentPosition);
            if (!conditionComponent.getEnabled())
            {
                continue;
            }
            double totalHeat = conditionComponent.get(this->pModel->getTimeSolver().getCurrentTime());

            double measure = this->findElementGroupMeasure(*pElementGroup);
            if (measure < RConstants::eps)
            {
                RLogger::warning("Heat boundary condition on entity \'%s\' is ignored - the entity has no computable elements.\n",
                                 pElementGroup->getName().toUtf8().constData());
                continue;
            }

            for (unsigned int k=0;k<pElementGroup->size();k++)
            {
                unsigned int elementID = pElementGroup->get(k);
                if (!this->computableElements[elementID])
                {
                    continue;
                }
                heatValues[elementID] = totalHeat / measure;
                setValues[elementID] = true;
            }
        }
    }
}

double RSolverGeneric::findElementGroupMeasure(const RElementGroup &rElementGroup) const
{
    double measure = 0.0;

    for (unsigned int i=0;i<rElementGroup.size();i++)
    {
        unsigned int elementID = rElementGroup.get(i);
        if (!this->computableElements[elementID])
        {
            continue;
        }

        const RElement &rElement = this->pModel->getElement(elementID);
        double elementMeasure = 0.0;

        switch (RElementGroup::getGroupType(rElement.getType()))
        {
            case R_ENTITY_GROUP_POINT:
                // A point has no measure to spread the heat over - count it instead.
                elementMeasure = 1.0;
                break;
            case R_ENTITY_GROUP_LINE:
                rElement.findLength(this->pModel->getNodes(),elementMeasure);
                break;
            case R_ENTITY_GROUP_SURFACE:
                rElement.findArea(this->pModel->getNodes(),elementMeasure);
                break;
            case R_ENTITY_GROUP_VOLUME:
                rElement.findVolume(this->pModel->getNodes(),elementMeasure);
                break;
            default:
                break;
        }

        measure += elementMeasure;
    }

    return measure;
}

void RSolverGeneric::recoverVariable(RVariableType variableType,
                                     RVariableApplyType applyType,
                                     uint expectedSize,
                                     uint vectorPosition,
                                     RRVector &variableValues,
                                     double defaultValue) const
{
    variableValues.resize(expectedSize);
    variableValues.fill(defaultValue);

    uint variablePosition = this->pModel->findVariable(variableType);

    if (variablePosition != RConstants::eod)
    {
        const RVariable &rVariable = this->pModel->getVariable(variablePosition);
        if (rVariable.getApplyType() == applyType)
        {
            variableValues = rVariable.getValues(vectorPosition);
        }
    }

    R_ERROR_ASSERT(variableValues.size() == expectedSize);
}

void RSolverGeneric::syncShared(const QString &keyName, RRVector &values)
{
    if (this->taskIteration == 0)
    {
        this->pSharedData->addData(keyName, values);
    }
    else
    {
        if (this->pSharedData->hasData(keyName,uint(values.size())))
        {
            values = this->pSharedData->findData(keyName);
        }
    }
}

void RSolverGeneric::findComputableElements(RProblemType problemType)
{
    this->computableElements.resize(this->pModel->getNElements());
    this->computableElements.fill(false);

    QList<RMaterialProperty::Type> propList;

    bool needsMaterial = (RMaterialProperty::getTypes(problemType).size() > 0);

    if (needsMaterial)
    {
        for (uint i=0;i<RMaterialProperty::nTypes;i++)
        {
            if (problemType & RMaterialProperty::getProblemTypeMask(RMaterialProperty::Type(i)))
            {
                propList.push_back(RMaterialProperty::Type(i));
            }
        }
    }
    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        const RMaterial &rMaterial = pElementGroup->getMaterial();
        if (!needsMaterial || rMaterial.hasProperties(propList))
        {
            for (uint j=0;j<pElementGroup->size();j++)
            {
                this->computableElements[pElementGroup->get(j)] = true;
            }
            continue;
        }
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &rBoundaryCondition = pElementGroup->getBoundaryCondition(j);
            if (problemType & RBoundaryCondition::getProblemTypeMask(rBoundaryCondition.getType()))
            {
                for (uint k=0;k<pElementGroup->size();k++)
                {
                    this->computableElements[pElementGroup->get(k)] = true;
                }
            }
        }
    }
}

void RSolverGeneric::findIncludableElements()
{
    this->includableElements = this->computableElements;

    RBVector computableNodes(this->pModel->getNNodes(),false);

    for (uint i=0;i<this->computableElements.size();i++)
    {
        if (this->computableElements[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            for (uint j=0;j<rElement.size();j++)
            {
                computableNodes[rElement.getNodeId(j)] = true;
            }
        }
    }

    for (uint i=0;i<this->includableElements.size();i++)
    {
        if (!this->includableElements[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            uint nComputable = 0;
            for (uint j=0;j<rElement.size();j++)
            {
                if (computableNodes[rElement.getNodeId(j)])
                {
                    nComputable++;
                }
            }
            if (nComputable == rElement.size())
            {
                this->includableElements[i] = true;
            }
        }
    }

}

void RSolverGeneric::findInwardElements()
{
    this->inwardElements.resize(this->pModel->getNElements());
    this->inwardElements.fill(false);

    // Build node → computable volume element index to avoid O(N²) search below.
    std::vector<std::vector<uint>> nodeToVolElems(this->pModel->getNNodes());
    for (uint l=0;l<this->pModel->getNElements();l++)
    {
        if (!this->computableElements[l])
        {
            continue;
        }
        const RElement &rEl = this->pModel->getElement(l);
        if (!R_ELEMENT_TYPE_IS_VOLUME(rEl.getType()))
        {
            continue;
        }
        for (uint m=0;m<rEl.size();m++)
        {
            nodeToVolElems[rEl.getNodeId(m)].push_back(l);
        }
    }

    // Assuming that surface normals are synchronized.

    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);

        bool isInward = false;

        for (uint j=0;j<rSurface.size();j++)
        {
            uint elementID = rSurface.get(j);
            uint volumeElementID = RConstants::eod;
            if (!this->includableElements[elementID])
            {
                continue;
            }
            const RElement &rElement = this->pModel->getElement(elementID);
            // Find related volume element via node index (O(degree) instead of O(N_elements)).
            const std::vector<uint> &candidates = nodeToVolElems[rElement.getNodeId(0)];
            for (uint l : candidates)
            {
                const RElement &rVolumeElement = this->pModel->getElement(l);
                uint nNodesFound = 0;
                for (uint m=0;m<rElement.size();m++)
                {
                    if (rVolumeElement.hasNodeId(rElement.getNodeId(m)))
                    {
                        nNodesFound++;
                    }
                }
                if (nNodesFound == rElement.size())
                {
                    volumeElementID = l;
                    break;
                }
            }
            if (volumeElementID != RConstants::eod)
            {
                const RElement &rVolumeElement = this->pModel->getElement(volumeElementID);

                std::vector<RElement> edgeElements = rVolumeElement.generateEdgeElements();

                RR3Vector elementNormal;
                rElement.findNormal(this->pModel->getNodes(),elementNormal[0],elementNormal[1],elementNormal[2]);

                for (uint l=0;l<edgeElements.size();l++)
                {
                    uint nNodesFound = 0;
                    for (uint m=0;m<rElement.size();m++)
                    {
                        if (edgeElements[l].hasNodeId(rElement.getNodeId(m)))
                        {
                            nNodesFound++;
                        }
                    }
                    if (nNodesFound == rElement.size())
                    {
                        RR3Vector edgeElementNormal;
                        edgeElements[l].findNormal(this->pModel->getNodes(),edgeElementNormal[0],edgeElementNormal[1],edgeElementNormal[2]);

                        if (std::fabs(RR3Vector::angle(elementNormal,edgeElementNormal)) < RConstants::pi/2.0)
                        {
                            isInward = true;
                        }
                    }
                }

                break;
            }
        }
        for (uint j=0;j<rSurface.size();j++)
        {
            this->inwardElements[rSurface.get(j)] = isInward;
        }
    }
}

void RSolverGeneric::processMonitoringPoints() const
{
    for (uint i=0;i<this->pModel->getMonitoringPointManager().size();i++)
    {
        RMonitoringPoint &rMonitorinPoint = this->pModel->getMonitoringPointManager()[i];

        unsigned int variableID = this->pModel->findVariable(rMonitorinPoint.getVariableType());
        if (variableID == RConstants::eod)
        {
            RLogger::warning("Monitoring point variable (%s) was not found among the results\n",
                             RVariable::getName(rMonitorinPoint.getVariableType()).toUtf8().constData());
            continue;
        }
        RVariable &rVariable = this->pModel->getVariable(variableID);

        RNode iNode(rMonitorinPoint.getPosition());

        unsigned int elementID = RConstants::eod;
        for (unsigned int j=0;j<this->pModel->getNElements();j++)
        {
            if (this->pModel->getElement(j).isInside(this->pModel->getNodes(),iNode))
            {
                elementID = j;
                break;
            }
        }
        if (elementID == RConstants::eod)
        {
            RLogger::warning("Monitoring point [%g %g %g] is outside of the model\n",
                             rMonitorinPoint.getPosition()[0],
                             rMonitorinPoint.getPosition()[1],
                             rMonitorinPoint.getPosition()[2]);
            continue;
        }

        RElement &rElement = this->pModel->getElement(elementID);

        RValueVector valueVector;
        valueVector.resize(rVariable.getNVectors());

        if (rVariable.getApplyType() == R_VARIABLE_APPLY_ELEMENT)
        {
            RRVector values(rVariable.getValueVector(elementID));
            for (unsigned int j=0;j<values.size();j++)
            {
                valueVector[j] = values[j];
            }
        }
        else if (rVariable.getApplyType() == R_VARIABLE_APPLY_NODE)
        {
            std::vector<RRVector> values(rElement.size());
            for (unsigned int j=0;j<rElement.size();j++)
            {
                values[j] = rVariable.getValueVector(rElement.getNodeId(j));
            }

            for (unsigned int j=0;j<rVariable.getNVectors();j++)
            {
                RRVector nodeValues(rElement.size());
                for (unsigned int k=0;k<rElement.size();k++)
                {
                    nodeValues[k] = values[k][j];
                }
                valueVector[j] = rElement.interpolate(this->pModel->getNodes(),iNode,nodeValues);
            }
        }

        try
        {
            const QString outputFileName = rMonitorinPoint.getOutputFileName();
            const QFileInfo beforeInfo(outputFileName);
            const qint64 previousSize = beforeInfo.exists() ? beforeInfo.size() : 0;
            rMonitorinPoint.writeValueVectorToFile(this->pModel->getTimeSolver().getCurrentTimeStep()+1,valueVector);
            if (this->monitoringOutputFailuresAreFatal && !outputFileName.isEmpty())
            {
                const QFileInfo afterInfo(outputFileName);
                if (!afterInfo.isFile() || afterInfo.size() <= previousSize)
                {
                    throw RError(RError::Type::WriteFile,R_ERROR_REF,
                                 "Monitoring point output did not append a record to '%s'.",
                                 outputFileName.toUtf8().constData());
                }
            }
        }
        catch (const RError &error)
        {
            if (this->monitoringOutputFailuresAreFatal &&
                !rMonitorinPoint.getOutputFileName().isEmpty())
            {
                throw;
            }
            RLogger::error("Failed to write monitoring point. %s\n",error.getMessage().toUtf8().constData());
        }
    }
}

void RSolverGeneric::printStats(RVariableType variableType) const
{
    unsigned int variablePosition = this->pModel->findVariable(variableType);
    if (variablePosition == RConstants::eod)
    {
        return;
    }
    const RVariable &rVariable = this->pModel->getVariable(variablePosition);

    RStatistics stat(rVariable.getValues(),100,true);

    RLogger::info("%s [%s]\n",rVariable.getName().toUtf8().constData(),rVariable.getUnits().toUtf8().constData());
    RLogger::indent();
    stat.print();
    RLogger::unindent();
}
