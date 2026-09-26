#include <cmath>
#include <array>
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include <rbl_error.h>

#include <rml_line.h>
#include <rml_material.h>
#include <rml_model.h>
#include <rml_structural_control_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_interaction.h>
#include <rml_variable.h>

#include "rstructural_cable_model_adapter.h"
#include "rstructural_cable_dynamics_restart_adapter.h"
#include "rstructural_cable_kernel.h"
#include "rstructural_static_truss_model.h"

namespace
{

struct CableModelData
{
    RStructuralStaticTrussModelData common;
    std::vector<RStructuralCableEquilibriumMember> members;
    std::vector<double> nodalLumpedMass;
};

[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}

constexpr double referenceTemperature = 293.15;

void validateDynamicsSetup(const RStructuralDynamicsSetup &setup)
{
    setup.validate();
    if (!setup.getEnabled())
    {
        fail("Cable transient adapter requires enabled structural dynamics.");
    }
    if (setup.getMethod() != RStructuralTimeIntegrationMethod::Newmark ||
        std::fabs(setup.getNewmarkGamma()-0.5) > 1.0e-12 ||
        std::fabs(setup.getNewmarkBeta()-0.25) > 1.0e-12)
    {
        fail("Cable transient mode supports Newmark average acceleration only (gamma=0.5, beta=0.25).");
    }
    if (setup.getRayleighAlpha() != 0.0 || setup.getRayleighBeta() != 0.0 ||
        !setup.getModalDampingRatios().empty())
    {
        fail("Cable transient mode does not support Rayleigh or modal damping.");
    }
}

void validateCableControlCatalog(const RModel &model,bool dynamics)
{
    const RStructuralControlSetup &control =
        model.getProblemSetup().getStructuralControlSetup();
    if (!control.enabled) return;
    if (!dynamics)
        fail("Static cable mode does not support Stage-6 structural control.");
    control.validate();

    const RStructuralInteractionSetup &interaction =
        model.getProblemSetup().getStructuralInteractionSetup();
    const RStructuralValidationReport report = interaction.validate(model);
    if (report.hasErrors())
    {
        for (const RStructuralDiagnostic &diagnostic : report.getDiagnostics())
        {
            if (diagnostic.severity == RStructuralDiagnosticSeverity::Error)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable control catalog is invalid (%s): %s",
                             diagnostic.objectId.toUtf8().constData(),
                             diagnostic.message.toUtf8().constData());
            }
        }
        fail("Cable control catalog is invalid.");
    }
    if (interaction.getActuators().empty())
        fail("Enabled cable control requires at least one point-node global NodalForce actuator.");

    for (const RStructuralActuator &actuator : interaction.getActuators())
    {
        if (actuator.type != RStructuralActuatorType::NodalForce ||
            actuator.target.scope != RStructuralTargetScope::Node ||
            actuator.frame != RStructuralFrame::Global)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control actuator '%s' requires point-node global NodalForce v1.",
                         actuator.id.toUtf8().constData());
        }
        if (actuator.variable != R_VARIABLE_FORCE &&
            actuator.variable != R_VARIABLE_FORCE_X &&
            actuator.variable != R_VARIABLE_FORCE_Y &&
            actuator.variable != R_VARIABLE_FORCE_Z)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control actuator '%s' command must be a physical NodalForce component.",
                         actuator.id.toUtf8().constData());
        }
    }
    for (const RStructuralSensor &sensor : interaction.getSensors())
    {
        if (sensor.target.scope != RStructuralTargetScope::Node ||
            sensor.aggregation != RStructuralSensorAggregation::Point ||
            sensor.frame != RStructuralFrame::Global ||
            (sensor.quantity != RStructuralSensorQuantity::Displacement &&
             sensor.quantity != RStructuralSensorQuantity::Velocity &&
             sensor.quantity != RStructuralSensorQuantity::Acceleration))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control sensor '%s' requires point-node global displacement, velocity or acceleration v1.",
                         sensor.id.toUtf8().constData());
        }
    }
}

void clearControlCatalogForTrussExtraction(RModel &model)
{
    model.getProblemSetup().getStructuralControlSetup() = RStructuralControlSetup();
    RStructuralInteractionSetup &interaction =
        model.getProblemSetup().getStructuralInteractionSetup();
    interaction.getActuators().clear();
    interaction.getSensors().clear();
}

CableModelData buildCableModelData(const RModel &model, bool dynamics)
{
    const RProblemSetup &problemSetup = model.getProblemSetup();
    const RStructuralCableSetup &setup = problemSetup.getStructuralCableSetup();
    if (!setup.enabled)
    {
        fail("Cable model adapter requires an enabled structural cable setup.");
    }
    if (problemSetup.getStructuralPiezoSetup().enabled)
    {
        fail("Structural cable and piezo modes cannot be enabled together.");
    }
    setup.validate(model);

    const RStructuralDynamicsSetup &dynamicsSetup =
        problemSetup.getStructuralDynamicsSetup();
    if (dynamics)
    {
        validateDynamicsSetup(dynamicsSetup);
    }
    else if (dynamicsSetup.getEnabled())
    {
        fail("Static cable mode does not support second-order dynamics.");
    }
    validateCableControlCatalog(model,dynamics);

    CableModelData data;
    if (dynamics &&
        (problemSetup.getRestart() ||
         problemSetup.getStructuralDynamicsRestartData().available ||
         problemSetup.getStructuralControlSetup().enabled))
    {
        // The shared static truss extractor intentionally rejects restart and
        // Stage-6 interaction data. Build physical FE data from a local clean
        // copy while retaining every accepted value on the live model.
        RModel physicalModel(model);
        physicalModel.getProblemSetup().setRestart(false);
        physicalModel.getProblemSetup().getStructuralDynamicsRestartData().clear();
        if (problemSetup.getStructuralControlSetup().enabled)
            clearControlCatalogForTrussExtraction(physicalModel);
        data.common = RStructuralStaticTrussModel::build(physicalModel,true,false);
    }
    else
    {
        data.common = RStructuralStaticTrussModel::build(model,dynamics,false);
    }

    std::vector<const RStructuralCableMember*> setupByElement(model.getNElements(),nullptr);
    for (const RStructuralCableMember &member : setup.members)
    {
        if (member.elementId >= setupByElement.size())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable setup references invalid element %u.",member.elementId);
        }
        if (setupByElement[member.elementId])
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable setup configures element %u more than once.",member.elementId);
        }
        setupByElement[member.elementId] = &member;
    }

    data.members.reserve(data.common.members.size());
    for (const RStructuralStaticTrussMember &member : data.common.members)
    {
        const RStructuralCableMember *cableMember = setupByElement[member.elementId];
        if (!cableMember)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable TRUSS1 element %u has no explicit cable setup record.",member.elementId);
        }
        RStructuralCableEquilibriumMember equilibriumMember;
        equilibriumMember.node1 = member.node1;
        equilibriumMember.node2 = member.node2;
        // The per-cable setup referenceLength is constitutive L0. The shared
        // member referenceLength remains the imported geometric mesh length.
        equilibriumMember.referenceLength = cableMember->referenceLength;
        equilibriumMember.initialTension = cableMember->initialTension;
        equilibriumMember.axialRigidity = member.area*member.youngModulus;
        if (!std::isfinite(equilibriumMember.axialRigidity) || equilibriumMember.axialRigidity <= 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable TRUSS1 element %u has invalid axial rigidity EA.",member.elementId);
        }
        data.members.push_back(equilibriumMember);
    }

    if (dynamics)
    {
        for (std::size_t dof=0;dof<data.common.prescribedMask.size();dof++)
        {
            if (data.common.prescribedMask[dof] &&
                std::fabs(data.common.prescribedDisplacements[dof]) > 1.0e-12)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient mode supports fixed zero-valued displacement supports only.");
            }
        }

        std::vector<double> densityByElement(model.getNElements(),0.0);
        for (uint lineId=0;lineId<model.getNLines();lineId++)
        {
            const RLine &line = model.getLine(lineId);
            if (line.empty())
            {
                continue;
            }
            const RMaterial &material = line.getMaterial();
            const uint densityPosition = material.findPosition(RMaterialProperty::Density);
            if (densityPosition >= material.size())
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient line group %u has no mass density.",lineId);
            }
            const RMaterialProperty &densityProperty = material.get(densityPosition);
            if (densityProperty.size() == 0)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient line group %u has an empty mass-density table.",lineId);
            }
            const double density = densityProperty.get(referenceTemperature);
            if (!std::isfinite(density) || density <= 0.0)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient line group %u requires finite positive density at 293.15 K.",lineId);
            }
            for (uint member=0;member<line.size();member++)
            {
                const uint elementId = line.get(member);
                if (elementId >= densityByElement.size())
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Cable transient line group %u references invalid element %u.",
                                 lineId,elementId);
                }
                densityByElement[elementId] = density;
            }
        }

        data.nodalLumpedMass.assign(data.common.nodes.size(),0.0);
        for (const RStructuralStaticTrussMember &member : data.common.members)
        {
            const double density = densityByElement[member.elementId];
            const double memberMass = density*member.area*member.referenceLength;
            if (!std::isfinite(memberMass) || memberMass <= 0.0)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient element %u has invalid lumped mass.",member.elementId);
            }
            const double halfMass = 0.5*memberMass;
            data.nodalLumpedMass[member.node1] += halfMass;
            data.nodalLumpedMass[member.node2] += halfMass;
            if (!std::isfinite(data.nodalLumpedMass[member.node1]) ||
                !std::isfinite(data.nodalLumpedMass[member.node2]))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable transient nodal lumped mass is non-finite.");
            }
        }
    }
    return data;
}

void storeNodalKinematics(RModel &model,RVariableType type,
                          const std::vector<double> &values)
{
    const std::size_t expected = std::size_t(3)*model.getNNodes();
    if (values.size() != expected)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Converged cable transient kinematics do not match the current node count.");
    }
    for (double value : values)
    {
        if (!std::isfinite(value))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Converged cable transient kinematics contain a non-finite value.");
        }
    }
    uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        position = model.addVariable(type);
    }
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(3,model.getNNodes());
    for (uint node=0;node<model.getNNodes();node++)
    {
        variable.setValue(0,node,values[3*node]);
        variable.setValue(1,node,values[3*node+1]);
        variable.setValue(2,node,values[3*node+2]);
    }
}

void copyResultVariable(RModel &model,const RModel &source,RVariableType type)
{
    const uint sourcePosition = source.findVariable(type);
    if (sourcePosition == RConstants::eod)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Cable transient mechanical publisher omitted a result variable.");
    }
    uint targetPosition = model.findVariable(type);
    if (targetPosition == RConstants::eod)
    {
        targetPosition = model.addVariable(type);
    }
    model.getVariable(targetPosition) = source.getVariable(sourcePosition);
}

void validateTransientResult(const RModel &model,
                             const RStructuralCableDynamicsResult &result,
                             std::size_t memberCount)
{
    const std::size_t expectedDofs = std::size_t(3)*model.getNNodes();
    if (result.displacements.size() != expectedDofs ||
        result.velocities.size() != expectedDofs ||
        result.accelerations.size() != expectedDofs ||
        result.internalForces.size() != expectedDofs ||
        result.elementTensions.size() != memberCount)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Converged cable transient response dimensions do not match the current model.");
    }
    const std::vector<double> *vectors[] =
    {
        &result.displacements,&result.velocities,&result.accelerations,&result.internalForces
    };
    for (const std::vector<double> *values : vectors)
    {
        for (double value : *values)
        {
            if (!std::isfinite(value))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Converged cable transient response contains a non-finite nodal value.");
            }
        }
    }
    for (double tension : result.elementTensions)
    {
        if (!std::isfinite(tension) || tension < 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Converged cable transient member tension must be finite and nonnegative.");
        }
    }
}

RStructuralCableDynamicsOptions cableDynamicsOptions(const RModel &model)
{
    const RStructuralCableSetup &cable =
        model.getProblemSetup().getStructuralCableSetup();
    const RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    RStructuralCableDynamicsOptions options;
    options.dt = dynamics.getTimeStep();
    options.beta = dynamics.getNewmarkBeta();
    options.gamma = dynamics.getNewmarkGamma();
    options.maxIterations = cable.maxIterations;
    options.maxLineSearchIterations = cable.maxLineSearchIterations;
    options.absForceTol = cable.absoluteForceTolerance;
    options.relForceTol = cable.relativeForceTolerance;
    return options;
}

QString cableDynamicsSystemFingerprint(const RModel &model,
                                       const CableModelData &data,
                                       const RStructuralCableDynamicsOptions &options)
{
    std::vector<double> elementAreas;
    elementAreas.reserve(data.common.members.size());
    for (const RStructuralStaticTrussMember &member : data.common.members)
    {
        elementAreas.push_back(member.area);
    }
    return RStructuralCableDynamicsRestartAdapter::computeSystemFingerprint(
        model,data.common.nodes,data.members,elementAreas,data.nodalLumpedMass,
        data.common.prescribedMask,data.common.prescribedDisplacements,
        data.common.deadNodalLoads,options);
}

std::vector<double> totalCableExternalForces(
    const CableModelData &data,const std::vector<double> &additionalNodalForces)
{
    std::vector<double> forces = data.common.deadNodalLoads;
    if (!additionalNodalForces.empty() &&
        additionalNodalForces.size() != forces.size())
    {
        fail("Cable transient additional nodal-force vector does not match the xyz DOF count.");
    }
    for (std::size_t dof=0;dof<additionalNodalForces.size();dof++)
    {
        const double controlForce = additionalNodalForces[dof];
        if (!std::isfinite(controlForce))
            fail("Cable transient additional nodal force is non-finite.");
        forces[dof] += controlForce;
        if (!std::isfinite(forces[dof]))
            fail("Cable transient total nodal force is non-finite.");
    }
    return forces;
}

void refreshCableAcceptedAcceleration(
    const CableModelData &data,
    const std::vector<double> &displacements,
    const std::vector<double> &externalForces,
    std::vector<double> &accelerations)
{
    const std::size_t nDofs = std::size_t(3)*data.common.nodes.size();
    if (displacements.size() != nDofs || accelerations.size() != nDofs ||
        externalForces.size() != nDofs || data.common.prescribedMask.size() != nDofs)
    {
        fail("Cable transient accepted acceleration refresh received inconsistent state dimensions.");
    }

    std::vector<double> internalForces(nDofs,0.0);
    for (const RStructuralCableEquilibriumMember &member : data.members)
    {
        RStructuralCableVector3 node1{},node2{};
        for (std::size_t axis=0;axis<3;axis++)
        {
            node1[axis] = data.common.nodes[member.node1][axis]+
                          displacements[std::size_t(member.node1)*3+axis];
            node2[axis] = data.common.nodes[member.node2][axis]+
                          displacements[std::size_t(member.node2)*3+axis];
        }
        const RStructuralCableKernelResult response = RStructuralCableKernel::evaluate(
            node1,node2,member.referenceLength,member.initialTension,
            member.axialRigidity);
        const std::size_t dofs[6] = {
            std::size_t(member.node1)*3,std::size_t(member.node1)*3+1,
            std::size_t(member.node1)*3+2,std::size_t(member.node2)*3,
            std::size_t(member.node2)*3+1,std::size_t(member.node2)*3+2
        };
        for (std::size_t local=0;local<6;local++)
        {
            internalForces[dofs[local]] += response.internalForce[local];
            if (!std::isfinite(internalForces[dofs[local]]))
                fail("Cable transient accepted internal force is non-finite during acceleration refresh.");
        }
    }

    for (std::size_t dof=0;dof<nDofs;dof++)
    {
        if (data.common.prescribedMask[dof])
        {
            accelerations[dof] = 0.0;
            continue;
        }
        const double mass = data.nodalLumpedMass[dof/3];
        if (!std::isfinite(mass) || mass <= 0.0)
            fail("Cable transient free DOF needs positive lumped mass for acceleration refresh.");
        accelerations[dof] = (externalForces[dof]-internalForces[dof])/mass;
        if (!std::isfinite(accelerations[dof]))
            fail("Cable transient right-limit acceleration is non-finite after a control-force change.");
    }
}

}

void RStructuralCableModelAdapter::validate(const RModel &model)
{
    if (!model.getProblemSetup().getStructuralCableSetup().enabled)
    {
        return;
    }
    const bool dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    const CableModelData data = buildCableModelData(model,dynamics);
    if (dynamics)
    {
        const RStructuralCableDynamicsOptions options = cableDynamicsOptions(model);
        const QString fingerprint = cableDynamicsSystemFingerprint(model,data,options);
        RStructuralCableDynamicsRestartAdapter::validateAcceptedState(
            model,fingerprint,data.common.deadNodalLoads);
    }
}

QString RStructuralCableModelAdapter::computeDynamicsSystemFingerprint(const RModel &model)
{
    const CableModelData data = buildCableModelData(model,true);
    return cableDynamicsSystemFingerprint(model,data,cableDynamicsOptions(model));
}

std::vector<bool> RStructuralCableModelAdapter::prescribedDofMask(const RModel &model)
{
    return buildCableModelData(model,true).common.prescribedMask;
}

RStructuralCableEquilibriumResult RStructuralCableModelAdapter::solve(const RModel &model)
{
    const CableModelData data = buildCableModelData(model,false);
    const RStructuralCableSetup &setup = model.getProblemSetup().getStructuralCableSetup();
    RStructuralCableEquilibriumOptions options;
    options.maxIterations = setup.maxIterations;
    options.maxLineSearchIterations = setup.maxLineSearchIterations;
    options.absForceTol = setup.absoluteForceTolerance;
    options.relForceTol = setup.relativeForceTolerance;
    return RStructuralCableEquilibriumSolver::solve(
        data.common.nodes,data.members,data.common.prescribedMask,
        data.common.prescribedDisplacements,data.common.deadNodalLoads,options);
}

RStructuralCableDynamicsResult RStructuralCableModelAdapter::solveDynamics(
    const RModel &model,
    const std::vector<double> &previousDisplacements,
    const std::vector<double> &previousVelocities,
    const std::vector<double> &previousAccelerations,
    const std::vector<double> &additionalNodalForces,
    bool refreshAcceptedAcceleration)
{
    const CableModelData data = buildCableModelData(model,true);
    const RStructuralCableDynamicsOptions options = cableDynamicsOptions(model);
    const QString fingerprint = cableDynamicsSystemFingerprint(model,data,options);
    std::vector<double> acceptedDisplacements = previousDisplacements;
    std::vector<double> acceptedVelocities = previousVelocities;
    std::vector<double> acceptedAccelerations = previousAccelerations;
    RStructuralCableDynamicsRestartAdapter::restoreAcceptedState(
        model,fingerprint,data.common.deadNodalLoads,
        acceptedDisplacements,acceptedVelocities,acceptedAccelerations);

    const std::vector<double> externalForces =
        totalCableExternalForces(data,additionalNodalForces);
    if (refreshAcceptedAcceleration)
    {
        refreshCableAcceptedAcceleration(data,acceptedDisplacements,externalForces,
                                         acceptedAccelerations);
    }

    return RStructuralCableDynamicsSolver::step(
        data.common.nodes,data.members,data.nodalLumpedMass,
        acceptedDisplacements,acceptedVelocities,acceptedAccelerations,
        data.common.prescribedMask,data.common.prescribedDisplacements,
        externalForces,options);
}

void RStructuralCableModelAdapter::initializeFreshStartAcceleration(
    const RModel &model,
    const std::vector<double> &displacements,
    std::vector<double> &accelerations)
{
    const RProblemSetup &problemSetup = model.getProblemSetup();
    if (model.getTimeSolver().getCurrentTimeStep() != 0 ||
        problemSetup.getRestart() ||
        problemSetup.getStructuralDynamicsRestartData().available)
    {
        return;
    }

    const CableModelData data = buildCableModelData(model,true);
    const std::vector<double> externalForces = totalCableExternalForces(data,{});
    refreshCableAcceptedAcceleration(data,displacements,externalForces,accelerations);
}

void RStructuralCableModelAdapter::clearResults(RModel &model)
{
    RStructuralStaticTrussModel::clearResults(model);
}

void RStructuralCableModelAdapter::publish(
    RModel &model,
    const RStructuralCableEquilibriumResult &result)
{
    const CableModelData data = buildCableModelData(model,false);
    if (result.elementTensions.size() != data.common.members.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Converged cable tension count does not match the current model.");
    }
    for (double tension : result.elementTensions)
    {
        if (!std::isfinite(tension) || tension < 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Converged cable member tension must be finite and nonnegative.");
        }
    }
    RStructuralStaticTrussModel::publishMechanical(
        model,result.displacements,result.internalForces,result.elementTensions);
}

void RStructuralCableModelAdapter::publishDynamics(
    RModel &model,
    const RStructuralCableDynamicsResult &result)
{
    const CableModelData data = buildCableModelData(model,true);
    validateTransientResult(model,result,data.common.members.size());
    const RStructuralCableDynamicsOptions options = cableDynamicsOptions(model);
    const QString fingerprint = cableDynamicsSystemFingerprint(model,data,options);
    RStructuralDynamicsRestartData checkpoint =
        RStructuralCableDynamicsRestartAdapter::stageAcceptedState(
            model,result,fingerprint,data.common.deadNodalLoads);

    // The static publisher's common data builder rejects persisted restart
    // records. Keep the live accepted checkpoint untouched while publishing
    // outputs through a clean, temporary copy in that case.
    const RProblemSetup &problemSetup = model.getProblemSetup();
    if (problemSetup.getRestart() ||
        problemSetup.getStructuralDynamicsRestartData().available ||
        problemSetup.getStructuralControlSetup().enabled)
    {
        RModel publicationModel(model);
        publicationModel.getProblemSetup().setRestart(false);
        publicationModel.getProblemSetup().getStructuralDynamicsRestartData().clear();
        if (problemSetup.getStructuralControlSetup().enabled)
            clearControlCatalogForTrussExtraction(publicationModel);
        RStructuralStaticTrussModel::publishMechanical(
            publicationModel,result.displacements,result.internalForces,
            result.elementTensions,true);
        const RVariableType publishedTypes[] = {
            R_VARIABLE_DISPLACEMENT,
            R_VARIABLE_FORCE,
            R_VARIABLE_STRESS_NORMAL,
            R_VARIABLE_STRESS_VON_MISES,
            R_VARIABLE_STRESS_X
        };
        for (RVariableType type : publishedTypes)
        {
            copyResultVariable(model,publicationModel,type);
        }
    }
    else
    {
        RStructuralStaticTrussModel::publishMechanical(
            model,result.displacements,result.internalForces,result.elementTensions,true);
    }
    storeNodalKinematics(model,R_VARIABLE_VELOCITY,result.velocities);
    storeNodalKinematics(model,R_VARIABLE_ACCELERATION,result.accelerations);
    RStructuralCableDynamicsRestartAdapter::commitAcceptedState(
        model,std::move(checkpoint));
}
