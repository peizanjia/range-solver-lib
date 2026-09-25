#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

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
    if ((this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
         this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled) && callback)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control callbacks are unsupported in static cable/piezo mode.");
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
    if ((this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
         this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled) && observer)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control observers are unsupported in static cable/piezo mode.");
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

void RSolver::run()
{
    const RProblemTaskItem &taskTree = this->pModel->getProblemTaskTree();
    if (isHeatStressSequence(taskTree) &&
        (taskTree.getNIterations() != 1 || this->pModel->getTimeSolver().getEnabled()))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transactional Heat->Stress task publication supports one steady root iteration only.");
    }

    if (this->pModel->getProblemSetup().getStructuralCableSetup().enabled ||
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled)
    {
        // Recheck because the referenced RModel remains mutable after _init().
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
        if (this->structuralSupportResponseObserver || this->structuralControlCallback ||
            this->structuralControlStepObserver)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural support/control observers are unsupported in static cable/piezo mode.");
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
    timeSolver.harmonizeTimesWithInput(this->pModel->getProblemSetup().getRestart());

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
        if (this->pModel->getProblemSetup().getRestart())
        {
            timeSolver.setNextTimeStep();
        }
        // The legacy pre-pass rewrites prior records and removes future ones
        // before a trial step is accepted. Structural dynamics publishes its
        // checkpoint/results transactionally; preserve existing files on a
        // failed step and let RModel::write publish the new latest link only
        // after its target record succeeds.
        if (!dynamics.getEnabled())
        {
            RSolverGeneric::updateOldRecords(timeSolver,this->modelFileName);
        }

        do {
            RLogger::info("time step: %9u of %-9u | time = % 12e [sec] | dt = % 12e [sec]\n",
                          timeSolver.getCurrentTimeStep()+1,
                          timeSolver.getNTimeSteps(),
                          timeSolver.getCurrentTime(),
                          timeSolver.getCurrentTimeStepSize());
            timeSolver.setComputedTime(timeSolver.getCurrentTime());

            RLogger::indent();
            this->runSingle();
            RLogger::unindent();
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
    if (taskTree.getNIterations() != 1 || this->pModel->getTimeSolver().getEnabled())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Transactional Heat->Stress task publication supports one steady root iteration only.");
    }

    const RModel acceptedModel(*this->pModel);
    const RSolverSharedData acceptedSharedData(this->sharedData);
    const QMap<RProblemType,uint> acceptedExecutionCount(this->solversExecutionCount);
    heatSolver->setResultWritesDeferred(true);
    stressSolver->setResultWritesDeferred(true);
    try
    {
        if (!this->runProblemTask(taskTree,0))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Sequential Heat->Stress task did not converge; transaction rolled back.");
        }

        // Both sub-solvers have populated the same in-memory model. Publish it
        // exactly once, after the enclosing task reports convergence.
        stressSolver->publishDeferredResults();
        heatSolver->setResultWritesDeferred(false);
        stressSolver->setResultWritesDeferred(false);
    }
    catch (...)
    {
        const std::exception_ptr originalFailure = std::current_exception();
        heatSolver->setResultWritesDeferred(false);
        stressSolver->setResultWritesDeferred(false);
        try
        {
            *this->pModel = acceptedModel;
            this->sharedData = acceptedSharedData;
            this->solversExecutionCount = acceptedExecutionCount;
        }
        catch (...)
        {
            RLogger::warning("Heat->Stress task rollback failed while restoring its snapshot.\n");
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
