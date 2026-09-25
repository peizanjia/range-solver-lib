#include <cmath>
#include <array>
#include <algorithm>
#include <limits>
#include <vector>

#include <rbl_error.h>

#include <rml_line.h>
#include <rml_material.h>
#include <rml_model.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_structural_cable_setup.h>
#include <rml_variable.h>

#include "rstructural_cable_model_adapter.h"
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

    CableModelData data;
    data.common = RStructuralStaticTrussModel::build(model,dynamics,false);

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

}

void RStructuralCableModelAdapter::validate(const RModel &model)
{
    if (!model.getProblemSetup().getStructuralCableSetup().enabled)
    {
        return;
    }
    const bool dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    (void)buildCableModelData(model,dynamics);
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
    const std::vector<double> &previousAccelerations)
{
    const CableModelData data = buildCableModelData(model,true);
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

    return RStructuralCableDynamicsSolver::step(
        data.common.nodes,data.members,data.nodalLumpedMass,
        previousDisplacements,previousVelocities,previousAccelerations,
        data.common.prescribedMask,data.common.prescribedDisplacements,
        data.common.deadNodalLoads,options);
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

    // The common publisher validates the complete mechanical response and the
    // current model before its first write. Kinematic vectors were validated
    // above, so the remaining writes cannot expose a partial failed candidate.
    RStructuralStaticTrussModel::publishMechanical(
        model,result.displacements,result.internalForces,result.elementTensions,true);
    storeNodalKinematics(model,R_VARIABLE_VELOCITY,result.velocities);
    storeNodalKinematics(model,R_VARIABLE_ACCELERATION,result.accelerations);
}
