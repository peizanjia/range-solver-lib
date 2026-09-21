#include <QFile>
#include <QFileInfo>
#include <QTextStream>

#include <cmath>

#include <rbl_logger.h>
#include <rbl_utils.h>

#include "riterationinfo.h"

void RIterationInfo::_init(const RIterationInfo *pIterationInfo)
{
    if (pIterationInfo)
    {
        this->nIterations = pIterationInfo->nIterations;
        this->convergenceValue = pIterationInfo->convergenceValue;
        this->outputFrequency = pIterationInfo->outputFrequency;
        this->iteration = pIterationInfo->iteration;
        this->error = pIterationInfo->error;
        this->trend = pIterationInfo->trend;
        this->equationScale = pIterationInfo->equationScale;
        this->outputFileName = pIterationInfo->outputFileName;
    }
}

RIterationInfo::RIterationInfo()
    : nIterations(0)
    , convergenceValue(RConstants::eps)
    , outputFrequency(1)
    , iteration(0)
    , error(0.0)
    , trend(0.0)
    , equationScale(1.0)
{
    this->_init();
}

RIterationInfo::RIterationInfo(const RIterationInfo &iterationInfo)
{
    this->_init(&iterationInfo);
}

RIterationInfo::~RIterationInfo()
{
}

RIterationInfo &RIterationInfo::operator =(const RIterationInfo &iterationInfo)
{
    this->_init(&iterationInfo);
    return (*this);
}

void RIterationInfo::setNIterations(unsigned int nIterations)
{
    this->nIterations = nIterations;
}

void RIterationInfo::setOutputFrequency(unsigned int outputFrequency)
{
    this->outputFrequency = outputFrequency;
}

void RIterationInfo::setConvergenceValue(double convergenceValue)
{
    this->convergenceValue = convergenceValue;
}

void RIterationInfo::setIteration(unsigned int iteration)
{
    this->iteration = iteration;
}

void RIterationInfo::setError(double error)
{
    this->trend = error - this->error;
    this->error = error;
}

void RIterationInfo::setEquationScale(double equationScale)
{
    this->equationScale = equationScale;
}

void RIterationInfo::setOutputFileName(const QString &outputFileName)
{
    this->outputFileName = outputFileName;
}

bool RIterationInfo::hasConverged() const
{
    if (this->hasDiverged())
    {
        return false;
    }
    if (this->error < this->convergenceValue)
    {
        return true;
    }
    return false;
}

bool RIterationInfo::hasDiverged() const
{
    return (!std::isfinite(this->error) || !std::isfinite(this->trend));
}

void RIterationInfo::printHeader(const QString &title) const
{
    if (this->outputFrequency == 0)
    {
        return;
    }
    RLogger::info("-------------------------------------------------------------------------------\n");
    RLogger::info("%s\n",title.toUtf8().constData());
    RLogger::info("-------------------------------------------------------------------------------\n");
    RLogger::info(">    Iteration     |   Residual   |    Trend     |                             \n");
    RLogger::info("-------------------------------------------------------------------------------\n");
}

void RIterationInfo::printFooter() const
{
    if (this->outputFrequency == 0)
    {
        return;
    }
    RLogger::info("-------------------------------------------------------------------------------\n");
}

void RIterationInfo::printIteration() const
{
    if (this->outputFrequency == 0)
    {
        return;
    }
    if (!this->hasConverged() && this->iteration % this->outputFrequency != 0)
    {
        return;
    }
    RLogger::info("> %6u of %-6u |% 12e |% 12e |\n",this->iteration+1,this->nIterations, this->error, this->trend);
    this->writeToFile();
}

void RIterationInfo::writeToFile(const QString &fileName, uint iteration, const std::vector<RIterationInfoValue> &iterationValues)
{
    if (fileName.isEmpty())
    {
        return;
    }
    bool firstWrite = !QFileInfo::exists(fileName);

    QFile file(fileName);

    if (!file.open(QIODevice::Append | QIODevice::Text))
    {
        RLogger::warning("Failed to open the file \'%s\' for writing.\n",fileName.toUtf8().constData());
        return;
    }

    QTextStream out(&file);
    out.setRealNumberNotation(QTextStream::ScientificNotation);
    out.setRealNumberPrecision(7);

    if (firstWrite)
    {
        out << "Iteration number";
        for (uint i=0;i<iterationValues.size();i++)
        {
            out << "|" << iterationValues[i].getName();
        }
        out << "\n";
    }

    out << iteration + 1;
    for (uint i=0;i<iterationValues.size();i++)
    {
        out << "|" << iterationValues[i].getValue();
    }
    out << "\n";

    file.close();
}

void RIterationInfo::writeToFile() const
{
    std::vector<RIterationInfoValue> cvgValues;
    cvgValues.push_back(RIterationInfoValue("Solver residual",this->error));
    cvgValues.push_back(RIterationInfoValue("Solver convergence",this->trend));
    RIterationInfo::writeToFile(this->outputFileName,this->iteration,cvgValues);
}
