#include "rstructural_cable_dynamics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_rvector.h>
#include <rml_sparse_matrix.h>

namespace
{
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();
constexpr double armijoFraction = 1.0e-4;
constexpr double lineSearchContraction = 0.5;

struct DynamicEvaluation
{
    std::vector<double> internalForces;
    std::vector<double> residual;
    std::vector<double> freeResidual;
    std::vector<double> elementTensions;
    RSparseMatrix effectiveTangent;
    double residualNorm = 0.0;
    double internalForceNorm = 0.0;
    double inertiaForceNorm = 0.0;
    double externalForceNorm = 0.0;
    double strainEnergy = 0.0;
    double freePotential = 0.0;
    double potentialScale = 0.0;
};

double stableNorm(const std::vector<double> &values)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            throw std::runtime_error("Cable dynamics encountered a non-finite vector value");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        throw std::runtime_error("Cable dynamics vector norm is non-finite");
    return norm;
}

double dot(const std::vector<double> &left, const std::vector<double> &right)
{
    if (left.size() != right.size())
        throw std::logic_error("Cable dynamics vector dimensions do not match");
    double value = 0.0;
    for (std::size_t i=0;i<left.size();i++)
    {
        value += left[i]*right[i];
        if (!std::isfinite(value))
            throw std::runtime_error("Cable dynamics dot product is non-finite");
    }
    return value;
}

void validateInputs(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<double> &nodalLumpedMass,
    const std::vector<double> &previousDisplacements,
    const std::vector<double> &previousVelocities,
    const std::vector<double> &previousAccelerations,
    const std::vector<bool> &prescribedMask,
    const std::vector<double> &prescribedDisplacements,
    const std::vector<double> &deadNodalLoads,
    const RStructuralCableDynamicsOptions &options)
{
    if (referenceNodes.empty() ||
        referenceNodes.size() > std::numeric_limits<std::size_t>::max()/3)
        throw std::invalid_argument("Cable dynamics requires a valid nonempty node array");
    const std::size_t dofCount = referenceNodes.size()*3;
    if (dofCount > std::numeric_limits<uint>::max() ||
        nodalLumpedMass.size() != referenceNodes.size() ||
        previousDisplacements.size() != dofCount || previousVelocities.size() != dofCount ||
        previousAccelerations.size() != dofCount || prescribedMask.size() != dofCount ||
        prescribedDisplacements.size() != dofCount || deadNodalLoads.size() != dofCount)
        throw std::invalid_argument("Cable dynamics histories, constraints, loads, and masses have invalid dimensions");

    if (!std::isfinite(options.dt) || options.dt <= 0.0 ||
        !std::isfinite(options.beta) || options.beta <= 0.0 ||
        !std::isfinite(options.gamma) || options.gamma < 0.5 ||
        options.beta < 0.25*(options.gamma+0.5)*(options.gamma+0.5) ||
        options.maxIterations == 0 || options.maxLineSearchIterations == 0 ||
        !std::isfinite(options.absForceTol) || options.absForceTol < 0.0 ||
        !std::isfinite(options.relForceTol) || options.relForceTol < 0.0 ||
        (options.absForceTol == 0.0 && options.relForceTol == 0.0))
        throw std::invalid_argument("Invalid cable dynamics Newmark or convergence options");

    for (const auto &node : referenceNodes)
    {
        for (double value : node)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Cable dynamics reference coordinates must be finite");
        }
    }
    for (double mass : nodalLumpedMass)
    {
        if (!std::isfinite(mass) || mass < 0.0)
            throw std::invalid_argument("Cable dynamics lumped masses must be finite and nonnegative");
    }
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!std::isfinite(previousDisplacements[dof]) ||
            !std::isfinite(previousVelocities[dof]) ||
            !std::isfinite(previousAccelerations[dof]) ||
            !std::isfinite(prescribedDisplacements[dof]) ||
            !std::isfinite(deadNodalLoads[dof]))
            throw std::invalid_argument("Cable dynamics displacement, velocity, acceleration, and load data must be finite");
    }
    for (const auto &member : members)
    {
        if (member.node1 >= referenceNodes.size() || member.node2 >= referenceNodes.size() ||
            member.node1 == member.node2 || !std::isfinite(member.referenceLength) ||
            member.referenceLength <= 0.0 || !std::isfinite(member.initialTension) ||
            member.initialTension < 0.0 || !std::isfinite(member.axialRigidity) ||
            member.axialRigidity <= 0.0)
            throw std::invalid_argument("Cable dynamics member has invalid nodes or material data");
        const double axialStiffness = member.axialRigidity/member.referenceLength;
        if (!std::isfinite(axialStiffness) || axialStiffness <= 0.0)
            throw std::invalid_argument("Cable dynamics member EA/L0 must be positive and finite");
    }
}

RStructuralCableVector3 currentPosition(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<double> &displacements,
    std::size_t node)
{
    RStructuralCableVector3 position{};
    for (std::size_t axis=0;axis<3;axis++)
    {
        const std::size_t dof = node*3+axis;
        position[axis] = referenceNodes[node][axis]+displacements[dof];
        if (!std::isfinite(position[axis]))
            throw std::invalid_argument("Cable dynamics current coordinates are non-finite");
    }
    return position;
}

DynamicEvaluation evaluateNetwork(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<double> &nodalLumpedMass,
    const std::vector<double> &displacements,
    const std::vector<double> &predictor,
    const std::vector<double> &deadNodalLoads,
    const std::vector<std::size_t> &globalToFree,
    std::size_t freeDofCount,
    double inverseBetaDtSquared,
    bool assembleTangent)
{
    DynamicEvaluation state;
    state.internalForces.assign(displacements.size(),0.0);
    state.residual.assign(displacements.size(),0.0);
    state.freeResidual.assign(freeDofCount,0.0);
    state.elementTensions.reserve(members.size());
    if (assembleTangent) state.effectiveTangent.setNRows(uint(freeDofCount));

    for (const auto &member : members)
    {
        const auto response = RStructuralCableKernel::evaluate(
            currentPosition(referenceNodes,displacements,member.node1),
            currentPosition(referenceNodes,displacements,member.node2),
            member.referenceLength,member.initialTension,member.axialRigidity);
        state.elementTensions.push_back(response.tension);
        state.strainEnergy += response.energy;
        state.potentialScale += response.energy;
        if (!std::isfinite(state.strainEnergy) || !std::isfinite(state.potentialScale))
            throw std::runtime_error("Cable dynamics strain energy is non-finite");

        const std::array<std::size_t,6> globalDofs = {{
            member.node1*3,member.node1*3+1,member.node1*3+2,
            member.node2*3,member.node2*3+1,member.node2*3+2
        }};
        for (std::size_t local=0;local<6;local++)
        {
            const std::size_t global = globalDofs[local];
            state.internalForces[global] += response.internalForce[local];
            if (!std::isfinite(state.internalForces[global]))
                throw std::runtime_error("Cable dynamics internal force is non-finite");
        }

        if (assembleTangent)
        {
            for (std::size_t row=0;row<6;row++)
            {
                const std::size_t freeRow = globalToFree[globalDofs[row]];
                if (freeRow == noFreeDof) continue;
                for (std::size_t column=0;column<6;column++)
                {
                    const std::size_t freeColumn = globalToFree[globalDofs[column]];
                    if (freeColumn == noFreeDof) continue;
                    const double value = response.tangent[row][column];
                    if (value == 0.0) continue;
                    const double sum = state.effectiveTangent.findValue(
                        uint(freeRow),uint(freeColumn))+value;
                    if (!std::isfinite(sum))
                        throw std::runtime_error("Cable dynamics tangent assembly is non-finite");
                    state.effectiveTangent.addValue(uint(freeRow),uint(freeColumn),value);
                }
            }
        }
    }

    for (std::size_t dof=0;dof<displacements.size();dof++)
    {
        const std::size_t node = dof/3;
        const double inertia = nodalLumpedMass[node]*inverseBetaDtSquared*
                               (displacements[dof]-predictor[dof]);
        if (!std::isfinite(inertia))
            throw std::runtime_error("Cable dynamics inertial force is non-finite");
        state.residual[dof] = state.internalForces[dof]+inertia-deadNodalLoads[dof];
        if (!std::isfinite(state.residual[dof]))
            throw std::runtime_error("Cable dynamics residual is non-finite");
        const std::size_t freeIndex = globalToFree[dof];
        if (freeIndex != noFreeDof)
        {
            state.freeResidual[freeIndex] = state.residual[dof];
            state.internalForceNorm = std::hypot(state.internalForceNorm,state.internalForces[dof]);
            state.inertiaForceNorm = std::hypot(state.inertiaForceNorm,inertia);
            state.externalForceNorm = std::hypot(state.externalForceNorm,deadNodalLoads[dof]);
            const double inertialPotential = 0.5*nodalLumpedMass[node]*inverseBetaDtSquared*
                (displacements[dof]-predictor[dof])*(displacements[dof]-predictor[dof]);
            const double loadPotential = deadNodalLoads[dof]*displacements[dof];
            state.freePotential += inertialPotential-loadPotential;
            state.potentialScale += std::fabs(inertialPotential)+std::fabs(loadPotential);
            if (assembleTangent && nodalLumpedMass[node] > 0.0)
            {
                const double massTangent = nodalLumpedMass[node]*inverseBetaDtSquared;
                const double diagonal = state.effectiveTangent.findValue(
                    uint(freeIndex),uint(freeIndex))+massTangent;
                if (!std::isfinite(diagonal))
                    throw std::runtime_error("Cable dynamics mass tangent is non-finite");
                state.effectiveTangent.addValue(uint(freeIndex),uint(freeIndex),massTangent);
            }
        }
    }
    state.freePotential += state.strainEnergy;
    state.residualNorm = stableNorm(state.freeResidual);
    if (!std::isfinite(state.internalForceNorm) || !std::isfinite(state.inertiaForceNorm) ||
        !std::isfinite(state.externalForceNorm) || !std::isfinite(state.freePotential) ||
        !std::isfinite(state.potentialScale))
        throw std::runtime_error("Cable dynamics force scale or effective potential is non-finite");
    return state;
}

std::vector<double> sparseMultiply(const RSparseMatrix &matrix,
                                  const std::vector<double> &vector)
{
    RRVector x(uint(vector.size()),0.0);
    for (std::size_t i=0;i<vector.size();i++) x[uint(i)] = vector[i];
    RRVector y;
    RSparseMatrix::mlt(matrix,x,y);
    return std::vector<double>(y.begin(),y.end());
}

std::vector<double> solveJacobiCg(const RSparseMatrix &matrix,
                                  const std::vector<double> &rhs)
{
    const std::size_t n = rhs.size();
    std::vector<double> solution(n,0.0);
    if (n == 0) return solution;

    std::vector<double> diagonal(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        diagonal[i] = matrix.findValue(uint(i),uint(i));
        if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0)
            throw std::runtime_error("Cable effective tangent is singular or has a nonpositive diagonal");
    }

    std::vector<double> residual = rhs;
    const double rhsNorm = stableNorm(rhs);
    if (rhsNorm == 0.0) return solution;
    const double linearTolerance = 1.0e-10*rhsNorm;

    std::vector<double> preconditioned(n,0.0);
    std::vector<double> direction(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        preconditioned[i] = residual[i]/diagonal[i];
        direction[i] = preconditioned[i];
        if (!std::isfinite(preconditioned[i]))
            throw std::runtime_error("Cable effective tangent Jacobi preconditioner is non-finite");
    }
    double rho = dot(residual,preconditioned);
    if (!std::isfinite(rho) || rho <= 0.0)
        throw std::runtime_error("Cable effective tangent has no positive preconditioned residual");

    const std::size_t maxCgIterations = n >= 2500 ? 10000 : std::max<std::size_t>(100,4*n);
    for (std::size_t iteration=0;iteration<maxCgIterations;iteration++)
    {
        const std::vector<double> product = sparseMultiply(matrix,direction);
        const double curvature = dot(direction,product);
        if (!std::isfinite(curvature) || curvature <= 0.0)
            throw std::runtime_error("Cable effective tangent is singular or has nonpositive CG curvature");
        const double alpha = rho/curvature;
        if (!std::isfinite(alpha))
            throw std::runtime_error("Cable effective tangent CG step is non-finite");
        for (std::size_t i=0;i<n;i++)
        {
            solution[i] += alpha*direction[i];
            residual[i] -= alpha*product[i];
            if (!std::isfinite(solution[i]) || !std::isfinite(residual[i]))
                throw std::runtime_error("Cable effective tangent CG iterate is non-finite");
        }
        if (stableNorm(residual) <= linearTolerance) break;
        for (std::size_t i=0;i<n;i++)
        {
            preconditioned[i] = residual[i]/diagonal[i];
            if (!std::isfinite(preconditioned[i]))
                throw std::runtime_error("Cable effective tangent CG residual is non-finite");
        }
        const double nextRho = dot(residual,preconditioned);
        if (!std::isfinite(nextRho) || nextRho <= 0.0)
            throw std::runtime_error("Cable effective tangent CG lost positive residual norm");
        const double beta = nextRho/rho;
        if (!std::isfinite(beta))
            throw std::runtime_error("Cable effective tangent CG recurrence is non-finite");
        for (std::size_t i=0;i<n;i++) direction[i] = preconditioned[i]+beta*direction[i];
        rho = nextRho;
    }

    const std::vector<double> product = sparseMultiply(matrix,solution);
    std::vector<double> freshResidual(n,0.0);
    for (std::size_t i=0;i<n;i++) freshResidual[i] = rhs[i]-product[i];
    if (stableNorm(freshResidual) > 10.0*linearTolerance)
        throw std::runtime_error("Cable effective tangent CG failed fresh residual certification");
    return solution;
}

double kineticEnergy(const std::vector<double> &nodalLumpedMass,
                     const std::vector<double> &velocities)
{
    double energy = 0.0;
    for (std::size_t dof=0;dof<velocities.size();dof++)
    {
        const double component = 0.5*nodalLumpedMass[dof/3]*velocities[dof]*velocities[dof];
        energy += component;
        if (!std::isfinite(component) || !std::isfinite(energy))
            throw std::runtime_error("Cable dynamics kinetic energy is non-finite");
    }
    return energy;
}
}

RStructuralCableDynamicsResult RStructuralCableDynamicsSolver::step(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<double> &nodalLumpedMass,
    const std::vector<double> &previousDisplacements,
    const std::vector<double> &previousVelocities,
    const std::vector<double> &previousAccelerations,
    const std::vector<bool> &prescribedMask,
    const std::vector<double> &prescribedDisplacements,
    const std::vector<double> &deadNodalLoads,
    const RStructuralCableDynamicsOptions &options)
{
    validateInputs(referenceNodes,members,nodalLumpedMass,previousDisplacements,
                   previousVelocities,previousAccelerations,prescribedMask,
                   prescribedDisplacements,deadNodalLoads,options);
    const std::size_t dofCount = referenceNodes.size()*3;
    const double dtSquared = options.dt*options.dt;
    const double betaDtSquared = options.beta*dtSquared;
    const double inverseBetaDtSquared = 1.0/betaDtSquared;
    if (!std::isfinite(dtSquared) || !std::isfinite(betaDtSquared) ||
        betaDtSquared <= 0.0 || !std::isfinite(inverseBetaDtSquared))
        throw std::invalid_argument("Cable dynamics Newmark time scale is non-finite");

    std::vector<double> predictor(dofCount,0.0);
    std::vector<double> velocityPredictor(dofCount,0.0);
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        predictor[dof] = previousDisplacements[dof]+options.dt*previousVelocities[dof]+
            dtSquared*(0.5-options.beta)*previousAccelerations[dof];
        velocityPredictor[dof] = previousVelocities[dof]+
            options.dt*(1.0-options.gamma)*previousAccelerations[dof];
        if (!std::isfinite(predictor[dof]) || !std::isfinite(velocityPredictor[dof]))
            throw std::invalid_argument("Cable dynamics Newmark predictor is non-finite");
    }

    std::vector<std::size_t> globalToFree(dofCount,noFreeDof);
    std::size_t freeDofCount = 0;
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!prescribedMask[dof]) globalToFree[dof] = freeDofCount++;
    }

    std::vector<double> displacements = predictor;
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (prescribedMask[dof]) displacements[dof] = prescribedDisplacements[dof];
    }
    std::size_t iterations = 0;
    for (;;)
    {
        const DynamicEvaluation state = evaluateNetwork(
            referenceNodes,members,nodalLumpedMass,displacements,predictor,
            deadNodalLoads,globalToFree,freeDofCount,inverseBetaDtSquared,true);
        const double forceScale = std::max({state.internalForceNorm,
                                            state.inertiaForceNorm,
                                            state.externalForceNorm});
        const double forceTolerance = options.absForceTol+options.relForceTol*forceScale;
        if (!std::isfinite(forceTolerance))
            throw std::runtime_error("Cable dynamics force tolerance is non-finite");
        if (state.residualNorm <= forceTolerance)
        {
            // Rebuild forces from the converged candidate and certify the
            // nonlinear residual again before exposing any result.
            const DynamicEvaluation fresh = evaluateNetwork(
                referenceNodes,members,nodalLumpedMass,displacements,predictor,
                deadNodalLoads,globalToFree,freeDofCount,inverseBetaDtSquared,false);
            const double freshScale = std::max({fresh.internalForceNorm,
                                                fresh.inertiaForceNorm,
                                                fresh.externalForceNorm});
            const double freshTolerance = options.absForceTol+options.relForceTol*freshScale;
            if (!std::isfinite(freshTolerance) || fresh.residualNorm > freshTolerance)
                throw std::runtime_error("Cable dynamics failed fresh force residual certification");

            RStructuralCableDynamicsResult result;
            result.displacements = displacements;
            result.accelerations.resize(dofCount,0.0);
            result.velocities.resize(dofCount,0.0);
            result.internalForces = fresh.internalForces;
            result.elementTensions = fresh.elementTensions;
            result.residualNorm = fresh.residualNorm;
            result.strainEnergy = fresh.strainEnergy;
            result.iterations = iterations;
            result.reactions = fresh.residual;
            for (std::size_t dof=0;dof<dofCount;dof++)
            {
                result.accelerations[dof] =
                    (displacements[dof]-predictor[dof])*inverseBetaDtSquared;
                result.velocities[dof] = velocityPredictor[dof]+
                    options.gamma*options.dt*result.accelerations[dof];
                if (!std::isfinite(result.accelerations[dof]) ||
                    !std::isfinite(result.velocities[dof]))
                    throw std::runtime_error("Cable dynamics candidate kinematics are non-finite");
                result.externalWorkIncrement += deadNodalLoads[dof]*
                    (displacements[dof]-previousDisplacements[dof]);
                if (!std::isfinite(result.externalWorkIncrement))
                    throw std::runtime_error("Cable dynamics external work is non-finite");
            }

            const std::vector<double> noLoads(dofCount,0.0);
            const DynamicEvaluation previous = evaluateNetwork(
                referenceNodes,members,nodalLumpedMass,previousDisplacements,
                previousDisplacements,noLoads,globalToFree,freeDofCount,0.0,false);
            result.previousStrainEnergy = previous.strainEnergy;
            result.previousKineticEnergy = kineticEnergy(nodalLumpedMass,previousVelocities);
            result.kineticEnergy = kineticEnergy(nodalLumpedMass,result.velocities);
            for (std::size_t dof=0;dof<dofCount;dof++)
            {
                if (!prescribedMask[dof]) continue;
                const double previousSupportForce =
                    previous.internalForces[dof]+
                    nodalLumpedMass[dof/3]*previousAccelerations[dof]-deadNodalLoads[dof];
                const double currentSupportForce = fresh.residual[dof];
                result.supportWorkIncrement += 0.5*(previousSupportForce+currentSupportForce)*
                    (displacements[dof]-previousDisplacements[dof]);
                if (!std::isfinite(result.supportWorkIncrement))
                    throw std::runtime_error("Cable dynamics prescribed-support work is non-finite");
            }
            result.energyBalanceError =
                (result.strainEnergy+result.kineticEnergy)-
                (result.previousStrainEnergy+result.previousKineticEnergy)-
                result.externalWorkIncrement-result.supportWorkIncrement;
            if (!std::isfinite(result.energyBalanceError))
                throw std::runtime_error("Cable dynamics energy balance diagnostic is non-finite");
            return result;
        }
        if (iterations >= options.maxIterations)
            throw std::runtime_error("Cable dynamic step failed to converge within the Newton iteration limit");

        std::vector<double> rhs(state.freeResidual.size(),0.0);
        for (std::size_t i=0;i<rhs.size();i++) rhs[i] = -state.freeResidual[i];
        const std::vector<double> increment = solveJacobiCg(state.effectiveTangent,rhs);
        const double descentSlope = dot(state.freeResidual,increment);
        if (!std::isfinite(descentSlope) || descentSlope >= 0.0)
            throw std::runtime_error("Cable dynamic Newton direction is not a descent direction");

        bool accepted = false;
        double stepLength = 1.0;
        for (std::size_t line=0;line<options.maxLineSearchIterations;line++)
        {
            std::vector<double> trial = displacements;
            for (std::size_t dof=0;dof<dofCount;dof++)
            {
                const std::size_t freeIndex = globalToFree[dof];
                if (freeIndex != noFreeDof)
                    trial[dof] += stepLength*increment[freeIndex];
            }
            const DynamicEvaluation trialState = evaluateNetwork(
                referenceNodes,members,nodalLumpedMass,trial,predictor,
                deadNodalLoads,globalToFree,freeDofCount,inverseBetaDtSquared,false);
            const double trialTolerance = options.absForceTol+
                options.relForceTol*std::max({trialState.internalForceNorm,
                                               trialState.inertiaForceNorm,
                                               trialState.externalForceNorm});
            const bool armijoSatisfied = trialState.freePotential <=
                state.freePotential+armijoFraction*stepLength*descentSlope;
            const double roundoffBound = 32.0*std::numeric_limits<double>::epsilon()*
                std::max(state.potentialScale,trialState.potentialScale);
            const bool unresolvedPotential = std::fabs(
                trialState.freePotential-state.freePotential) <= roundoffBound;
            const bool materiallyImproved = trialState.residualNorm <=
                0.9*state.residualNorm || trialState.residualNorm <= trialTolerance;
            if (armijoSatisfied || (unresolvedPotential && materiallyImproved))
            {
                displacements.swap(trial);
                accepted = true;
                break;
            }
            stepLength *= lineSearchContraction;
        }
        if (!accepted)
            throw std::runtime_error("Cable dynamic Newton line search failed to reduce the effective potential");
        iterations++;
    }
}
