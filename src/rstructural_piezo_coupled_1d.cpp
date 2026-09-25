#include "rstructural_piezo_coupled_1d.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <rml_matrix_solver_conf.h>
#include <rml_sparse_matrix.h>
#include <rbl_rvector.h>

#include "rmatrixpreconditioner.h"
#include "rmatrixsolver.h"
#include "rstructural_piezo_kernel.h"

namespace
{
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();

[[noreturn]] void invalidInput(const std::string &message)
{
    throw std::invalid_argument("Coupled 1-D piezo solve: " + message);
}

[[noreturn]] void failedSolve(const std::string &message)
{
    throw std::runtime_error("Coupled 1-D piezo solve: " + message);
}

double stableNorm(const std::vector<double> &values,const char *description)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            failedSolve(std::string(description) + " contains a non-finite value.");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        failedSolve(std::string(description) + " norm is non-finite.");
    return norm;
}

class DisjointSets
{
    private:
        std::vector<uint> parent;
        std::vector<uint> rank;

    public:
        explicit DisjointSets(uint size)
            : parent(size),rank(size,0)
        {
            std::iota(parent.begin(),parent.end(),uint(0));
        }

        uint find(uint value)
        {
            if (parent[value] != value)
                parent[value] = find(parent[value]);
            return parent[value];
        }

        void unite(uint first,uint second)
        {
            first = find(first);
            second = find(second);
            if (first == second)
                return;
            if (rank[first] < rank[second])
                std::swap(first,second);
            parent[second] = first;
            if (rank[first] == rank[second])
                ++rank[first];
        }
};

struct ValidatedInput
{
    std::vector<double> lengths;
    std::vector<double> directions;
    std::vector<uint> componentByNode;
    uint componentCount = 0;
};

ValidatedInput validateInput(const RStructuralPiezoCoupled1DInput &input)
{
    const std::size_t nNodes = input.coordinates.size();
    if (nNodes == 0)
        invalidInput("at least one node is required.");
    if (nNodes > std::numeric_limits<uint>::max()/2)
        invalidInput("the node count exceeds the sparse solver index limit.");
    if (input.members.empty())
        invalidInput("at least one piezo rod is required.");

    const auto requireNodeSize = [nNodes](std::size_t size,const char *name)
    {
        if (size != nNodes)
            invalidInput(std::string(name) + " must contain one entry per node.");
    };
    requireNodeSize(input.mechanicalLoads.size(),"mechanicalLoads");
    requireNodeSize(input.prescribedConductorCharges.size(),"prescribedConductorCharges");
    requireNodeSize(input.displacementPrescribed.size(),"displacementPrescribed");
    requireNodeSize(input.prescribedDisplacements.size(),"prescribedDisplacements");
    requireNodeSize(input.potentialPrescribed.size(),"potentialPrescribed");
    requireNodeSize(input.prescribedPotentials.size(),"prescribedPotentials");

    if (input.members.size() > std::numeric_limits<uint>::max())
        invalidInput("the member count exceeds the solver index limit.");
    if (input.gmresRestart == 0 || input.gmresCycles == 0 ||
        !std::isfinite(input.gmresRelativeTolerance) ||
        input.gmresRelativeTolerance <= 0.0 || input.gmresRelativeTolerance >= 1.0)
    {
        invalidInput("GMRES restart, cycle count and relative tolerance must be positive and finite.");
    }
    const std::array<double,4> physicalTolerances = {{
        input.mechanicalAbsoluteTolerance,input.mechanicalRelativeTolerance,
        input.electricAbsoluteTolerance,input.electricRelativeTolerance
    }};
    for (double tolerance : physicalTolerances)
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0)
            invalidInput("mechanical and electrical residual tolerances must be finite and nonnegative.");
    }
    if (input.mechanicalAbsoluteTolerance == 0.0 && input.mechanicalRelativeTolerance == 0.0)
        invalidInput("at least one mechanical residual tolerance must be positive.");
    if (input.electricAbsoluteTolerance == 0.0 && input.electricRelativeTolerance == 0.0)
        invalidInput("at least one electrical residual tolerance must be positive.");

    for (double coordinate : input.coordinates)
        if (!std::isfinite(coordinate))
            invalidInput("node coordinates must be finite.");

    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!std::isfinite(input.mechanicalLoads[node]) ||
            !std::isfinite(input.prescribedConductorCharges[node]))
        {
            invalidInput("nodal loads and conductor charges must be finite.");
        }
        if (!std::isfinite(input.prescribedDisplacements[node]) ||
            !std::isfinite(input.prescribedPotentials[node]))
        {
            invalidInput("prescribed displacements and potentials must be finite.");
        }
        if (input.potentialPrescribed[node] && input.prescribedConductorCharges[node] != 0.0)
        {
            invalidInput("a prescribed-potential node cannot also have a prescribed conductor charge.");
        }
    }

    ValidatedInput validated;
    validated.lengths.resize(input.members.size());
    validated.directions.resize(input.members.size());
    std::vector<bool> incident(nNodes,false);
    DisjointSets sets{uint(nNodes)};
    for (std::size_t memberIndex=0;memberIndex<input.members.size();memberIndex++)
    {
        const RStructuralPiezoCoupled1DMember &member = input.members[memberIndex];
        if (member.node0 >= nNodes || member.node1 >= nNodes || member.node0 == member.node1)
            invalidInput("each member must connect two distinct valid nodes.");
        if (!std::isfinite(member.area) || member.area <= 0.0 ||
            !std::isfinite(member.youngModulus) || member.youngModulus <= 0.0 ||
            !std::isfinite(member.coupling) || !std::isfinite(member.permittivity) ||
            member.permittivity <= 0.0)
        {
            invalidInput("member area, Young modulus and permittivity must be positive and finite; coupling must be finite.");
        }
        const double signedLength = input.coordinates[member.node1]-input.coordinates[member.node0];
        if (!std::isfinite(signedLength) || signedLength == 0.0)
            invalidInput("member endpoints must have distinct finite coordinates.");
        const double length = std::abs(signedLength);
        if (!std::isfinite(length) || length <= 0.0)
            invalidInput("member length must be positive and finite.");
        validated.lengths[memberIndex] = length;
        validated.directions[memberIndex] = std::copysign(1.0,signedLength);
        incident[member.node0] = true;
        incident[member.node1] = true;
        sets.unite(member.node0,member.node1);
    }

    std::vector<bool> hasMechanicalAnchor(nNodes,false);
    std::vector<bool> hasPotentialGauge(nNodes,false);
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!incident[node])
            invalidInput("isolated nodes are unsupported.");
        const uint component = sets.find(uint(node));
        if (input.displacementPrescribed[node])
            hasMechanicalAnchor[component] = true;
        if (input.potentialPrescribed[node])
            hasPotentialGauge[component] = true;
    }
    for (std::size_t node=0;node<nNodes;node++)
    {
        const uint component = sets.find(uint(node));
        if (!hasMechanicalAnchor[component])
            invalidInput("each connected component requires at least one prescribed displacement to remove rigid translation.");
        if (!hasPotentialGauge[component])
            invalidInput("each connected component requires at least one prescribed potential to fix the electric gauge.");
    }

    validated.componentByNode.resize(nNodes);
    std::vector<uint> rootToComponent(nNodes,std::numeric_limits<uint>::max());
    for (std::size_t node=0;node<nNodes;node++)
    {
        const uint root = sets.find(uint(node));
        if (rootToComponent[root] == std::numeric_limits<uint>::max())
            rootToComponent[root] = validated.componentCount++;
        validated.componentByNode[node] = rootToComponent[root];
    }
    return validated;
}

struct FreeSystem
{
    RSparseMatrix matrix;
    std::vector<double> rhs;
    std::vector<double> diagonalScale;
    std::vector<std::size_t> globalToFree;
    std::vector<std::size_t> freeToGlobal;
};

FreeSystem assembleFreeSystem(const RStructuralPiezoCoupled1DInput &input,
                              const ValidatedInput &validated)
{
    const std::size_t nNodes = input.coordinates.size();
    const std::size_t nDofs = 2*nNodes;
    FreeSystem system;
    system.globalToFree.assign(nDofs,noFreeDof);
    for (std::size_t dof=0;dof<nDofs;dof++)
    {
        const std::size_t node = dof % nNodes;
        const bool prescribed = dof < nNodes
            ? input.displacementPrescribed[node]
            : input.potentialPrescribed[node];
        if (!prescribed)
        {
            system.globalToFree[dof] = system.freeToGlobal.size();
            system.freeToGlobal.push_back(dof);
        }
    }
    if (system.freeToGlobal.size() > std::numeric_limits<uint>::max())
        invalidInput("the free DOF count exceeds the sparse solver index limit.");

    const uint nFree = uint(system.freeToGlobal.size());
    system.matrix.setNRows(nFree);
    system.rhs.assign(nFree,0.0);
    for (std::size_t freeDof=0;freeDof<system.freeToGlobal.size();freeDof++)
    {
        const std::size_t globalDof = system.freeToGlobal[freeDof];
        const std::size_t node = globalDof % nNodes;
        system.rhs[freeDof] = globalDof < nNodes
            ? input.mechanicalLoads[node]
            : -input.prescribedConductorCharges[node];
        if (!std::isfinite(system.rhs[freeDof]))
            invalidInput("the free mechanical/electric right-hand side is non-finite.");
    }

    for (std::size_t memberIndex=0;memberIndex<input.members.size();memberIndex++)
    {
        const RStructuralPiezoCoupled1DMember &member = input.members[memberIndex];
        const double direction = validated.directions[memberIndex];
        const double length = validated.lengths[memberIndex];
        const std::array<double,2> axialShape = {{-direction,direction}};
        const std::array<double,2> potentialShape = {{-1.0,1.0}};
        const std::array<std::size_t,4> globalDofs = {{
            std::size_t(member.node0),std::size_t(member.node1),
            nNodes+std::size_t(member.node0),nNodes+std::size_t(member.node1)
        }};
        const RStructuralPiezoKernelResult tangent = RStructuralPiezoKernel::evaluate(
            0.0,0.0,member.youngModulus,member.coupling,member.permittivity,
            member.area,length);
        std::array<std::array<double,4>,4> elementTangent{};
        for (std::size_t row=0;row<4;row++)
        {
            const bool rowIsMechanical = row < 2;
            const std::size_t rowShapeIndex = rowIsMechanical ? row : row-2;
            for (std::size_t column=0;column<4;column++)
            {
                const bool columnIsMechanical = column < 2;
                const std::size_t columnShapeIndex = columnIsMechanical ? column : column-2;
                double value = 0.0;
                if (rowIsMechanical && columnIsMechanical)
                {
                    value = tangent.hessian[0][0]*
                        axialShape[rowShapeIndex]*axialShape[columnShapeIndex];
                }
                else if (rowIsMechanical && !columnIsMechanical)
                {
                    value = tangent.hessian[0][1]*
                        axialShape[rowShapeIndex]*potentialShape[columnShapeIndex];
                }
                else if (!rowIsMechanical && columnIsMechanical)
                {
                    value = tangent.hessian[1][0]*
                        potentialShape[rowShapeIndex]*axialShape[columnShapeIndex];
                }
                else
                {
                    value = tangent.hessian[1][1]*
                        potentialShape[rowShapeIndex]*potentialShape[columnShapeIndex];
                }
                if (!std::isfinite(value))
                    invalidInput("an element tangent contains a non-finite entry.");
                elementTangent[row][column] = value;
            }
        }

        for (std::size_t row=0;row<4;row++)
        {
            const std::size_t globalRow = globalDofs[row];
            const std::size_t freeRow = system.globalToFree[globalRow];
            if (freeRow == noFreeDof)
                continue;
            for (std::size_t column=0;column<4;column++)
            {
                const std::size_t globalColumn = globalDofs[column];
                const double value = elementTangent[row][column];
                const std::size_t freeColumn = system.globalToFree[globalColumn];
                if (freeColumn == noFreeDof)
                {
                    const std::size_t prescribedNode = globalColumn % nNodes;
                    const double prescribedValue = globalColumn < nNodes
                        ? input.prescribedDisplacements[prescribedNode]
                        : input.prescribedPotentials[prescribedNode];
                    system.rhs[freeRow] -= value*prescribedValue;
                }
                else if (value != 0.0)
                {
                    system.matrix.addValue(uint(freeRow),uint(freeColumn),value);
                }
                if (!std::isfinite(system.rhs[freeRow]))
                    invalidInput("prescribed-value elimination produced a non-finite right-hand side.");
            }
        }
    }

    system.diagonalScale.resize(system.freeToGlobal.size(),0.0);
    for (std::size_t row=0;row<system.freeToGlobal.size();row++)
    {
        const double diagonal = system.matrix.findValue(uint(row),uint(row));
        if (!std::isfinite(diagonal) || diagonal == 0.0)
            invalidInput("a free DOF has a zero or non-finite tangent diagonal.");
        system.diagonalScale[row] = 1.0/std::sqrt(std::abs(diagonal));
        if (!std::isfinite(system.diagonalScale[row]) || system.diagonalScale[row] == 0.0)
            invalidInput("symmetric diagonal scaling is non-finite or zero.");
    }
    return system;
}

std::vector<double> prescribedValueVector(const RStructuralPiezoCoupled1DInput &input)
{
    const std::size_t nNodes = input.coordinates.size();
    std::vector<double> values(2*nNodes,0.0);
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (input.displacementPrescribed[node])
            values[node] = input.prescribedDisplacements[node];
        if (input.potentialPrescribed[node])
            values[nNodes+node] = input.prescribedPotentials[node];
    }
    return values;
}

RStructuralPiezoCoupled1DResult recoverResult(
    const RStructuralPiezoCoupled1DInput &input,
    const ValidatedInput &validated,
    const std::vector<double> &fullSolution)
{
    const std::size_t nNodes = input.coordinates.size();
    RStructuralPiezoCoupled1DResult result;
    result.displacements.assign(fullSolution.begin(),fullSolution.begin()+std::ptrdiff_t(nNodes));
    result.potentials.assign(fullSolution.begin()+std::ptrdiff_t(nNodes),fullSolution.end());
    result.internalForces.assign(nNodes,0.0);
    result.conductorCharges.assign(nNodes,0.0);
    result.memberAxialForces.resize(input.members.size(),0.0);
    result.memberRightCharges.resize(input.members.size(),0.0);

    for (std::size_t memberIndex=0;memberIndex<input.members.size();memberIndex++)
    {
        const RStructuralPiezoCoupled1DMember &member = input.members[memberIndex];
        const double direction = validated.directions[memberIndex];
        const double length = validated.lengths[memberIndex];
        const double elongation = direction*(
            result.displacements[member.node1]-result.displacements[member.node0]);
        const double voltage = result.potentials[member.node1]-result.potentials[member.node0];
        if (!std::isfinite(elongation) || !std::isfinite(voltage))
            failedSolve("recovered elongation or voltage is non-finite.");

        const RStructuralPiezoKernelResult response = RStructuralPiezoKernel::evaluate(
            elongation,voltage,member.youngModulus,member.coupling,
            member.permittivity,member.area,length);
        result.memberAxialForces[memberIndex] = response.mechanicalForce;
        result.memberRightCharges[memberIndex] = -response.dielectricDisplacementFlux;
        result.enthalpy += response.enthalpy;
        if (!std::isfinite(result.enthalpy))
            failedSolve("total electric enthalpy is non-finite.");

        const std::array<double,2> axialShape = {{-direction,direction}};
        const double leftCharge = response.dielectricDisplacementFlux;
        const double rightCharge = -response.dielectricDisplacementFlux;
        result.internalForces[member.node0] += axialShape[0]*response.mechanicalForce;
        result.internalForces[member.node1] += axialShape[1]*response.mechanicalForce;
        result.conductorCharges[member.node0] += leftCharge;
        result.conductorCharges[member.node1] += rightCharge;
        if (!std::isfinite(result.internalForces[member.node0]) ||
            !std::isfinite(result.internalForces[member.node1]) ||
            !std::isfinite(result.conductorCharges[member.node0]) ||
            !std::isfinite(result.conductorCharges[member.node1]))
        {
            failedSolve("nodal force or conductor charge recovery is non-finite.");
        }
    }

    std::vector<double> mechanicalResidual;
    std::vector<double> electricResidual;
    std::vector<double> freeMechanicalLoads;
    std::vector<double> freeMechanicalInternal;
    std::vector<double> freePrescribedCharges;
    std::vector<double> freeComputedCharges;
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!input.displacementPrescribed[node])
        {
            mechanicalResidual.push_back(result.internalForces[node]-input.mechanicalLoads[node]);
            freeMechanicalLoads.push_back(input.mechanicalLoads[node]);
            freeMechanicalInternal.push_back(result.internalForces[node]);
        }
        if (!input.potentialPrescribed[node])
        {
            electricResidual.push_back(result.conductorCharges[node]-input.prescribedConductorCharges[node]);
            freePrescribedCharges.push_back(input.prescribedConductorCharges[node]);
            freeComputedCharges.push_back(result.conductorCharges[node]);
        }
    }
    result.mechanicalResidualNorm = stableNorm(mechanicalResidual,"mechanical residual");
    result.electricResidualNorm = stableNorm(electricResidual,"electrical charge residual");
    const double mechanicalScale = std::max(
        stableNorm(freeMechanicalLoads,"free mechanical loads"),
        stableNorm(freeMechanicalInternal,"free mechanical internal forces"));
    const double electricScale = std::max(
        stableNorm(freePrescribedCharges,"free prescribed conductor charges"),
        stableNorm(freeComputedCharges,"free computed conductor charges"));
    const double mechanicalTolerance = input.mechanicalAbsoluteTolerance+
        input.mechanicalRelativeTolerance*mechanicalScale;
    const double electricTolerance = input.electricAbsoluteTolerance+
        input.electricRelativeTolerance*electricScale;
    if (!std::isfinite(mechanicalTolerance) || !std::isfinite(electricTolerance))
        failedSolve("physical residual tolerance overflowed.");
    if (result.mechanicalResidualNorm > mechanicalTolerance)
        failedSolve("fresh mechanical force residual failed its N-based certificate.");
    if (result.electricResidualNorm > electricTolerance)
        failedSolve("fresh electrical charge residual failed its C-based certificate.");

    // Element terminal charges cancel pairwise. Check that recovered nodal
    // conductor charges retain that exact charge balance in each component.
    std::vector<double> componentCharge(validated.componentCount,0.0);
    std::vector<double> componentAbsoluteCharge(validated.componentCount,0.0);
    std::vector<std::size_t> componentNodes(validated.componentCount,0);
    std::vector<std::size_t> componentMembers(validated.componentCount,0);
    for (std::size_t node=0;node<nNodes;node++)
    {
        const uint component = validated.componentByNode[node];
        componentCharge[component] += result.conductorCharges[node];
        componentAbsoluteCharge[component] += std::abs(result.conductorCharges[node]);
        ++componentNodes[component];
        if (!std::isfinite(componentCharge[component]) ||
            !std::isfinite(componentAbsoluteCharge[component]))
        {
            failedSolve("component charge balance is non-finite.");
        }
    }
    for (std::size_t member=0;member<input.members.size();member++)
        ++componentMembers[validated.componentByNode[input.members[member].node0]];
    for (std::size_t component=0;component<componentCharge.size();component++)
    {
        const double operationCount = double(componentNodes[component]+componentMembers[component]+1);
        const double roundoffBound = 4.0*operationCount*std::numeric_limits<double>::epsilon()*
            componentAbsoluteCharge[component]+std::numeric_limits<double>::denorm_min();
        if (std::abs(componentCharge[component]) > roundoffBound)
            failedSolve("recovered conductor charge does not balance over a connected component.");
    }
    return result;
}

void certifySparseResiduals(const FreeSystem &system,
                            const std::vector<double> &fullSolution,
                            const RStructuralPiezoCoupled1DInput &input,
                            const RStructuralPiezoCoupled1DResult &result)
{
    if (system.freeToGlobal.empty())
        return;

    RRVector freeSolution(uint(system.freeToGlobal.size()),0.0);
    for (std::size_t index=0;index<system.freeToGlobal.size();index++)
        freeSolution[uint(index)] = fullSolution[system.freeToGlobal[index]];
    RRVector product;
    RSparseMatrix::mlt(system.matrix,freeSolution,product);
    if (product.getNRows() != freeSolution.getNRows())
        failedSolve("fresh sparse residual has an invalid dimension.");

    const std::size_t nNodes = input.coordinates.size();
    std::vector<double> forceResidual;
    std::vector<double> chargeResidual;
    std::vector<double> freeMechanicalLoads;
    std::vector<double> freeMechanicalInternal;
    std::vector<double> freePrescribedCharges;
    std::vector<double> freeComputedCharges;
    for (std::size_t row=0;row<system.freeToGlobal.size();row++)
    {
        const std::size_t globalDof = system.freeToGlobal[row];
        const double residual = system.rhs[row]-product[uint(row)];
        if (!std::isfinite(residual))
            failedSolve("fresh sparse matrix residual is non-finite.");
        if (globalDof < nNodes)
        {
            forceResidual.push_back(residual);
            freeMechanicalLoads.push_back(input.mechanicalLoads[globalDof]);
            freeMechanicalInternal.push_back(result.internalForces[globalDof]);
        }
        else
        {
            const std::size_t node = globalDof-nNodes;
            chargeResidual.push_back(residual);
            freePrescribedCharges.push_back(input.prescribedConductorCharges[node]);
            freeComputedCharges.push_back(result.conductorCharges[node]);
        }
    }
    const double sparseForceNorm = stableNorm(forceResidual,"fresh sparse mechanical residual");
    const double sparseChargeNorm = stableNorm(chargeResidual,"fresh sparse electrical residual");
    const double forceScale = std::max(
        stableNorm(freeMechanicalLoads,"free mechanical loads"),
        stableNorm(freeMechanicalInternal,"free mechanical internal forces"));
    const double chargeScale = std::max(
        stableNorm(freePrescribedCharges,"free prescribed conductor charges"),
        stableNorm(freeComputedCharges,"free computed conductor charges"));
    const double forceTolerance = input.mechanicalAbsoluteTolerance+
        input.mechanicalRelativeTolerance*forceScale;
    const double chargeTolerance = input.electricAbsoluteTolerance+
        input.electricRelativeTolerance*chargeScale;
    if (!std::isfinite(forceTolerance) || !std::isfinite(chargeTolerance) ||
        sparseForceNorm > forceTolerance || sparseChargeNorm > chargeTolerance)
        failedSolve("fresh sparse residual failed its separate N/C certificates.");
}
}

RStructuralPiezoCoupled1DResult RStructuralPiezoCoupled1D::solve(
    const RStructuralPiezoCoupled1DInput &input)
{
    const ValidatedInput validated = validateInput(input);
    FreeSystem system = assembleFreeSystem(input,validated);
    const std::size_t nNodes = input.coordinates.size();
    std::vector<double> fullSolution = prescribedValueVector(input);

    if (!system.freeToGlobal.empty())
    {
        const std::size_t nFree = system.freeToGlobal.size();
        RSparseMatrix scaledMatrix;
        scaledMatrix.setNRows(uint(nFree));
        RRVector scaledRhs(uint(nFree),0.0);
        for (std::size_t row=0;row<nFree;row++)
        {
            scaledRhs[uint(row)] = system.diagonalScale[row]*system.rhs[row];
            if (!std::isfinite(scaledRhs[uint(row)]))
                invalidInput("symmetric scaling produced a non-finite right-hand side.");
            const RSparseVector<double> &entries = system.matrix.getVector(uint(row));
            for (uint position=0;position<entries.size();position++)
            {
                const uint column = entries.getIndex(position);
                const double value = entries.getValue(position)*
                    system.diagonalScale[row]*system.diagonalScale[column];
                if (!std::isfinite(value))
                    invalidInput("symmetric scaling produced a non-finite tangent entry.");
                if (value != 0.0)
                    scaledMatrix.addValue(uint(row),column,value);
            }
        }

        RMatrixSolverConf solverConf(RMatrixSolverConf::GMRES);
        solverConf.setNOuterIterations(input.gmresCycles);
        solverConf.setNInnerIterations(unsigned(std::min<std::size_t>(input.gmresRestart,nFree)));
        solverConf.setSolverCvgValue(input.gmresRelativeTolerance);
        solverConf.setOutputFrequency(input.gmresCycles);
        RMatrixSolver solver(solverConf);
        solver.disableConvergenceLogFile();
        RRVector scaledSolution;
        solver.solve(scaledMatrix,scaledRhs,scaledSolution,R_MATRIX_PRECONDITIONER_NONE,1);
        if (scaledSolution.getNRows() != nFree)
            failedSolve("GMRES returned a solution with an invalid dimension.");
        for (std::size_t free=0;free<nFree;free++)
        {
            const double value = system.diagonalScale[free]*scaledSolution[uint(free)];
            if (!std::isfinite(value))
                failedSolve("GMRES returned a non-finite solution.");
            fullSolution[system.freeToGlobal[free]] = value;
        }
    }

    for (double value : fullSolution)
        if (!std::isfinite(value))
            failedSolve("the full coupled solution contains a non-finite value.");

    RStructuralPiezoCoupled1DResult result = recoverResult(input,validated,fullSolution);
    certifySparseResiduals(system,fullSolution,input,result);
    return result;
}
