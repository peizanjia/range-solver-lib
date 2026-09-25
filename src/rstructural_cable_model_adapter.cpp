#include <cmath>
#include <vector>

#include <rbl_error.h>

#include <rml_model.h>
#include <rml_structural_cable_setup.h>

#include "rstructural_cable_model_adapter.h"
#include "rstructural_static_truss_model.h"

namespace
{

struct CableModelData
{
    RStructuralStaticTrussModelData common;
    std::vector<RStructuralCableEquilibriumMember> members;
};

[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}

CableModelData buildCableModelData(const RModel &model)
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

    CableModelData data;
    data.common = RStructuralStaticTrussModel::build(model);

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
    return data;
}

}

void RStructuralCableModelAdapter::validate(const RModel &model)
{
    if (!model.getProblemSetup().getStructuralCableSetup().enabled)
    {
        return;
    }
    (void)buildCableModelData(model);
}

RStructuralCableEquilibriumResult RStructuralCableModelAdapter::solve(const RModel &model)
{
    const CableModelData data = buildCableModelData(model);
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

void RStructuralCableModelAdapter::clearResults(RModel &model)
{
    RStructuralStaticTrussModel::clearResults(model);
}

void RStructuralCableModelAdapter::publish(
    RModel &model,
    const RStructuralCableEquilibriumResult &result)
{
    const CableModelData data = buildCableModelData(model);
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
