#include "rstructural_piezo_model_adapter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <rbl_error.h>
#include <rbl_rvector.h>

#include <rml_matrix_solver_conf.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_piezo_setup.h>
#include <rml_variable.h>

#include "rmatrixpreconditioner.h"
#include "rmatrixsolver.h"
#include "rstructural_piezo_kernel.h"
#include "rstructural_static_truss_model.h"

namespace
{
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();

[[noreturn]] void fail(const QString &message)
{
    const QByteArray utf8 = message.toUtf8();
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",utf8.constData());
}

double stableNorm(const std::vector<double> &values)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            fail("Piezo static solve encountered a non-finite vector value.");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        fail("Piezo static solve vector norm is non-finite.");
    return norm;
}

struct PiezoModelData
{
    RStructuralStaticTrussModelData truss;
    //! Port pointers in the same sorted-by-element order as truss.members.
    std::vector<const RStructuralPiezoPort*> portsByMember;
};

PiezoModelData buildPiezoModelData(const RModel &model)
{
    const RProblemSetup &problemSetup = model.getProblemSetup();
    const RStructuralPiezoSetup &setup = problemSetup.getStructuralPiezoSetup();
    if (!setup.enabled)
        fail("Piezo static adapter requires an enabled structural piezo setup.");
    if (problemSetup.getStructuralCableSetup().enabled)
        fail("Structural cable and piezo configurations cannot be enabled together.");

    setup.validate(model);
    PiezoModelData data;
    data.truss = RStructuralStaticTrussModel::build(model);
    if (setup.ports.size() != data.truss.members.size())
        fail("Piezo setup must configure exactly one terminal port for every TRUSS1 member.");

    std::map<uint,const RStructuralPiezoPort*> portByElement;
    for (const RStructuralPiezoPort &port : setup.ports)
    {
        if (!std::isfinite(port.coupling) || !std::isfinite(port.permittivity) ||
            port.permittivity <= 0.0 || !std::isfinite(port.voltage))
        {
            fail(QString("Piezo port '%1' has invalid coupling, permittivity or voltage.").arg(port.id));
        }
        if (!portByElement.emplace(port.elementId,&port).second)
            fail("Piezo setup has more than one port for a TRUSS1 member.");
    }

    data.portsByMember.reserve(data.truss.members.size());
    for (const RStructuralStaticTrussMember &member : data.truss.members)
    {
        const auto found = portByElement.find(member.elementId);
        if (found == portByElement.end())
        {
            fail(QString("TRUSS1 element %1 has no explicit piezo port.").arg(member.elementId));
        }
        data.portsByMember.push_back(found->second);
    }
    return data;
}

struct FreeSystem
{
    RSparseMatrix tangent;
    std::vector<double> rhs;
    std::vector<std::size_t> globalToFree;
    std::vector<std::size_t> freeToGlobal;
    double freeExternalLoadNorm = 0.0;
};

FreeSystem assembleFreeSystem(const PiezoModelData &data)
{
    const std::size_t dofCount = data.truss.prescribedMask.size();
    if (dofCount > std::numeric_limits<uint>::max())
        fail("Piezo model has too many DOFs for the sparse solver index type.");

    FreeSystem system;
    system.globalToFree.assign(dofCount,noFreeDof);
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!data.truss.prescribedMask[dof])
        {
            system.globalToFree[dof] = system.freeToGlobal.size();
            system.freeToGlobal.push_back(dof);
        }
    }
    if (system.freeToGlobal.size() > std::numeric_limits<uint>::max())
        fail("Piezo model has too many free DOFs for the sparse solver index type.");

    system.tangent.setNRows(uint(system.freeToGlobal.size()));
    system.rhs.assign(system.freeToGlobal.size(),0.0);
    for (std::size_t free=0;free<system.freeToGlobal.size();free++)
    {
        const double load = data.truss.deadNodalLoads[system.freeToGlobal[free]];
        if (!std::isfinite(load))
            fail("Piezo model contains a non-finite free nodal load.");
        system.rhs[free] = load;
        system.freeExternalLoadNorm = std::hypot(system.freeExternalLoadNorm,load);
    }
    if (!std::isfinite(system.freeExternalLoadNorm))
        fail("Piezo model free-load norm is non-finite.");

    for (std::size_t memberPosition=0;memberPosition<data.truss.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.truss.members[memberPosition];
        const RStructuralPiezoPort &port = *data.portsByMember[memberPosition];
        const auto &x1 = data.truss.nodes[member.node1];
        const auto &x2 = data.truss.nodes[member.node2];
        const double inverseLength = 1.0/member.referenceLength;
        const std::array<double,3> direction = {{
            (x2[0]-x1[0])*inverseLength,
            (x2[1]-x1[1])*inverseLength,
            (x2[2]-x1[2])*inverseLength
        }};
        const std::array<double,6> b = {{
            -direction[0],-direction[1],-direction[2],
             direction[0], direction[1], direction[2]
        }};
        const std::array<std::size_t,6> globalDofs = {{
            std::size_t(member.node1)*3,std::size_t(member.node1)*3+1,std::size_t(member.node1)*3+2,
            std::size_t(member.node2)*3,std::size_t(member.node2)*3+1,std::size_t(member.node2)*3+2
        }};
        const double axialStiffness = member.youngModulus*member.area*inverseLength;
        const double couplingScale = port.coupling*member.area*inverseLength;
        if (!std::isfinite(axialStiffness) || axialStiffness <= 0.0 ||
            !std::isfinite(couplingScale))
        {
            fail(QString("Piezo element %1 has non-finite axial or coupling stiffness.").arg(member.elementId));
        }

        for (std::size_t row=0;row<6;row++)
        {
            const std::size_t freeRow = system.globalToFree[globalDofs[row]];
            if (freeRow == noFreeDof)
                continue;

            // Internal force is K*u + b*(e*A/L)*V, so move the prescribed-
            // voltage term to the external side with a negative sign.
            system.rhs[freeRow] -= b[row]*couplingScale*port.voltage;
            if (!std::isfinite(system.rhs[freeRow]))
                fail("Piezo equivalent voltage load is non-finite.");

            for (std::size_t column=0;column<6;column++)
            {
                const double value = axialStiffness*b[row]*b[column];
                if (!std::isfinite(value))
                    fail("Piezo element tangent contains a non-finite entry.");
                if (value == 0.0)
                    continue;
                const std::size_t freeColumn = system.globalToFree[globalDofs[column]];
                if (freeColumn == noFreeDof)
                {
                    system.rhs[freeRow] -= value*data.truss.prescribedDisplacements[globalDofs[column]];
                }
                else
                {
                    const double updated = system.tangent.findValue(uint(freeRow),uint(freeColumn))+value;
                    if (!std::isfinite(updated))
                        fail("Piezo sparse tangent assembly is non-finite.");
                    system.tangent.addValue(uint(freeRow),uint(freeColumn),value);
                }
                if (!std::isfinite(system.rhs[freeRow]))
                    fail("Piezo prescribed-displacement load shift is non-finite.");
            }
        }
    }
    return system;
}

RStructuralPiezoModelResult recoverResult(const PiezoModelData &data,
                                          const std::vector<double> &displacements)
{
    RStructuralPiezoModelResult result;
    result.displacements = displacements;
    result.internalForces.assign(displacements.size(),0.0);
    result.memberAxialForces.resize(data.truss.members.size(),0.0);
    result.memberRightCharges.resize(data.truss.members.size(),0.0);

    for (std::size_t memberPosition=0;memberPosition<data.truss.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.truss.members[memberPosition];
        const RStructuralPiezoPort &port = *data.portsByMember[memberPosition];
        const auto &x1 = data.truss.nodes[member.node1];
        const auto &x2 = data.truss.nodes[member.node2];
        const double inverseLength = 1.0/member.referenceLength;
        const std::array<double,3> direction = {{
            (x2[0]-x1[0])*inverseLength,
            (x2[1]-x1[1])*inverseLength,
            (x2[2]-x1[2])*inverseLength
        }};
        const std::array<double,6> b = {{
            -direction[0],-direction[1],-direction[2],
             direction[0], direction[1], direction[2]
        }};
        const std::array<std::size_t,6> globalDofs = {{
            std::size_t(member.node1)*3,std::size_t(member.node1)*3+1,std::size_t(member.node1)*3+2,
            std::size_t(member.node2)*3,std::size_t(member.node2)*3+1,std::size_t(member.node2)*3+2
        }};
        double elongation = 0.0;
        for (std::size_t local=0;local<6;local++)
            elongation += b[local]*displacements[globalDofs[local]];
        if (!std::isfinite(elongation))
            fail(QString("Piezo element %1 has non-finite axial elongation.").arg(member.elementId));

        const RStructuralPiezoKernelResult response = RStructuralPiezoKernel::evaluate(
            elongation,port.voltage,member.youngModulus,port.coupling,
            port.permittivity,member.area,member.referenceLength);
        result.memberAxialForces[memberPosition] = response.mechanicalForce;
        result.memberRightCharges[memberPosition] = -response.dielectricDisplacementFlux;
        result.enthalpy += response.enthalpy;
        if (!std::isfinite(result.enthalpy))
            fail("Total piezo electric enthalpy is non-finite.");

        for (std::size_t local=0;local<6;local++)
        {
            const std::size_t dof = globalDofs[local];
            result.internalForces[dof] += b[local]*response.mechanicalForce;
            if (!std::isfinite(result.internalForces[dof]))
                fail("Piezo internal nodal force is non-finite.");
        }
    }

    std::vector<double> freeResidual;
    freeResidual.reserve(displacements.size());
    for (std::size_t dof=0;dof<displacements.size();dof++)
    {
        if (!data.truss.prescribedMask[dof])
        {
            freeResidual.push_back(result.internalForces[dof]-data.truss.deadNodalLoads[dof]);
        }
    }
    result.residualNorm = stableNorm(freeResidual);
    return result;
}

void storeElectricCharge(RModel &model,const std::vector<double> &elementCharges)
{
    const uint position = model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod
        ? model.addVariable(R_VARIABLE_ELECTRIC_CHARGE)
        : model.findVariable(R_VARIABLE_ELECTRIC_CHARGE);
    RVariable &charge = model.getVariable(position);
    charge.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    charge.resize(1,uint(elementCharges.size()));
    double minimum = 0.0;
    double maximum = 0.0;
    for (std::size_t element=0;element<elementCharges.size();element++)
    {
        const double value = elementCharges[element];
        charge.setValue(0,uint(element),value);
        minimum = std::min(minimum,value);
        maximum = std::max(maximum,value);
    }
    charge.getVariableData().setMinMaxDisplayValue(minimum,maximum);
}
}

void RStructuralPiezoModelAdapter::validate(const RModel &model)
{
    if (!model.getProblemSetup().getStructuralPiezoSetup().enabled)
        return;
    (void)buildPiezoModelData(model);
}

RStructuralPiezoModelResult RStructuralPiezoModelAdapter::solve(const RModel &model)
{
    const PiezoModelData data = buildPiezoModelData(model);
    FreeSystem system = assembleFreeSystem(data);
    const RMatrixSolverConf &solverConf = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    const double configuredTolerance = solverConf.getSolverCvgValue();
    if (!std::isfinite(configuredTolerance) || configuredTolerance < 0.0)
        fail("Piezo CG convergence tolerance must be finite and nonnegative.");

    std::vector<double> displacement(data.truss.prescribedDisplacements);
    if (!system.freeToGlobal.empty())
    {
        for (std::size_t row=0;row<system.freeToGlobal.size();row++)
        {
            const double diagonal = system.tangent.findValue(uint(row),uint(row));
            if (!std::isfinite(diagonal) || diagonal <= 0.0)
            {
                fail("Piezo free tangent has a nonpositive diagonal; it is singular or contains an unconstrained DOF.");
            }
        }

        RRVector rhs(uint(system.rhs.size()),0.0);
        for (std::size_t row=0;row<system.rhs.size();row++) rhs[uint(row)] = system.rhs[row];
        RMatrixSolver linearSolver(solverConf);
        linearSolver.disableConvergenceLogFile();
        RRVector freeSolution;
        linearSolver.solve(system.tangent,rhs,freeSolution,
                           R_MATRIX_PRECONDITIONER_JACOBI,1);
        if (freeSolution.getNRows() != rhs.getNRows())
            fail("Piezo CG returned a solution with invalid dimensions.");

        std::vector<double> freeValues(freeSolution.begin(),freeSolution.end());
        const double effectiveRhsNorm = stableNorm(system.rhs);
        const double linearTolerance =
            std::max(1.0e-8,configuredTolerance)*
            std::max(effectiveRhsNorm,system.freeExternalLoadNorm)+1.0e-12;
        if (!std::isfinite(linearTolerance))
            fail("Piezo linear residual tolerance is non-finite.");

        RRVector product;
        RSparseMatrix::mlt(system.tangent,freeSolution,product);
        if (product.getNRows() != rhs.getNRows())
            fail("Piezo fresh linear residual has invalid dimensions.");
        std::vector<double> freshLinearResidual(system.rhs.size(),0.0);
        for (std::size_t row=0;row<system.rhs.size();row++)
        {
            if (!std::isfinite(freeValues[row]) || !std::isfinite(product[uint(row)]))
                fail("Piezo CG returned a non-finite solution or matrix product.");
            freshLinearResidual[row] = system.rhs[row]-product[uint(row)];
        }
        if (stableNorm(freshLinearResidual) > linearTolerance)
            fail("Piezo CG failed fresh b-Au residual certification.");

        for (std::size_t free=0;free<system.freeToGlobal.size();free++)
            displacement[system.freeToGlobal[free]] = freeValues[free];
    }

    RStructuralPiezoModelResult result = recoverResult(data,displacement);
    const double effectiveRhsNorm = stableNorm(system.rhs);
    const double physicalTolerance = std::max(1.0e-8,configuredTolerance)*
        std::max(effectiveRhsNorm,system.freeExternalLoadNorm)+1.0e-12;
    if (!std::isfinite(physicalTolerance) || result.residualNorm > physicalTolerance)
        fail("Piezo recovered physical free-force residual failed certification.");
    return result;
}

void RStructuralPiezoModelAdapter::clearResults(RModel &model)
{
    RStructuralStaticTrussModel::clearResults(model);
}

void RStructuralPiezoModelAdapter::publish(
    RModel &model,const RStructuralPiezoModelResult &result)
{
    const PiezoModelData data = buildPiezoModelData(model);
    const std::size_t dofCount = std::size_t(model.getNNodes())*3;
    if (result.displacements.size() != dofCount || result.internalForces.size() != dofCount ||
        result.memberAxialForces.size() != data.truss.members.size() ||
        result.memberRightCharges.size() != data.truss.members.size() ||
        !std::isfinite(result.residualNorm) || result.residualNorm < 0.0 ||
        !std::isfinite(result.enthalpy))
    {
        fail("Piezo converged response dimensions or scalar values do not match the current model.");
    }
    for (double value : result.displacements)
        if (!std::isfinite(value)) fail("Piezo displacement contains a non-finite value.");
    for (double value : result.internalForces)
        if (!std::isfinite(value)) fail("Piezo internal force contains a non-finite value.");
    for (double value : result.memberAxialForces)
        if (!std::isfinite(value)) fail("Piezo axial force contains a non-finite value.");
    for (double value : result.memberRightCharges)
        if (!std::isfinite(value)) fail("Piezo right-electrode charge contains a non-finite value.");

    std::vector<double> elementCharges(model.getNElements(),0.0);
    for (std::size_t member=0;member<data.truss.members.size();member++)
    {
        const RStructuralStaticTrussMember &trussMember = data.truss.members[member];
        if (trussMember.elementId >= elementCharges.size())
            fail("Piezo charge result references an invalid element ID.");
        const double stress = result.memberAxialForces[member]/trussMember.area;
        if (!std::isfinite(stress))
            fail(QString("Piezo stress recovery is non-finite for element %1.").arg(trussMember.elementId));
        elementCharges[trussMember.elementId] = result.memberRightCharges[member];
    }
    for (double value : elementCharges)
        if (!std::isfinite(value)) fail("Piezo element charge output contains a non-finite value.");

    // All dimensions, recovered quantities and element slots are certified
    // before the shared mechanics helper or charge variable mutates the model.
    RStructuralStaticTrussModel::publishMechanical(
        model,result.displacements,result.internalForces,result.memberAxialForces);
    storeElectricCharge(model,elementCharges);
}
