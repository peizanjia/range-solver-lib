#include "rstructural_dynamics_integrator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "rmatrixsolver.h"

namespace
{
void addScaled(RRVector &destination, const RRVector &source, double scale)
{
    for (uint i=0; i<destination.size(); ++i)
    {
        destination[i] += scale*source[i];
    }
}

void addScaled(RSparseMatrix &destination, const RSparseMatrix &source, double scale)
{
    for (uint i=0; i<source.getNRows(); ++i)
    {
        const auto &row = source.getVector(i);
        for (uint j=0; j<row.size(); ++j)
        {
            destination.addValue(i,row.getIndex(j),scale*row.getValue(j));
        }
    }
}

RRVector product(const RSparseMatrix &matrix, const RRVector &vector)
{
    RRVector result;
    RSparseMatrix::mlt(matrix,vector,result);
    return result;
}

double quadraticEnergy(const RSparseMatrix &matrix, const RRVector &vector)
{
    return 0.5*RRVector::dot(vector,product(matrix,vector));
}

bool finiteVector(const RRVector &vector)
{
    for (double value : vector)
    {
        if (!std::isfinite(value)) return false;
    }
    return true;
}
}

RStructuralDynamicsIntegrator::RStructuralDynamicsIntegrator(const RSparseMatrix &mass,
                                                             const RSparseMatrix &damping,
                                                             const RSparseMatrix &stiffness,
                                                             const RMatrixSolverConf &solverConf,
                                                             Options options,
                                                             ForceCallback force)
    : mass(mass), damping(damping), stiffness(stiffness), solverConf(solverConf),
      options(options), force(std::move(force))
{
    const uint n = mass.getNRows();
    if (n == 0 || damping.getNRows() != n || stiffness.getNRows() != n)
        throw std::invalid_argument("Dynamics matrices must have the same nonzero row count");
    for (const RSparseMatrix *matrix : {&this->mass,&this->damping,&this->stiffness})
    {
        for (uint i=0; i<n; ++i)
        {
            const auto &row = matrix->getVector(i);
            for (uint j=0; j<row.size(); ++j)
            {
                if (row.getIndex(j) >= n || !std::isfinite(row.getValue(j)))
                    throw std::invalid_argument("Dynamics matrix has an invalid entry");
            }
        }
    }
    if (options.rhoInfinity < 0.0 || options.rhoInfinity > 1.0 || !std::isfinite(options.rhoInfinity))
        throw std::invalid_argument("Generalized-alpha spectral radius must be in [0,1]");
    if (!this->force)
        throw std::invalid_argument("Dynamics force callback is required");
    if (solverConf.getType() != RMatrixSolverConf::CG && solverConf.getType() != RMatrixSolverConf::GMRES)
        throw std::invalid_argument("A linear solver configuration is required");
    if (options.method == Method::CentralDifference)
    {
        double smallestMassMargin = std::numeric_limits<double>::infinity();
        double largestStiffnessRow = 0.0;
        for (uint row=0;row<n;row++)
        {
            double massMargin = this->mass.findValue(row,row);
            double stiffnessRow = 0.0;
            for (uint position=0;position<this->mass.getNColumns(row);position++)
            {
                if (this->mass.getVector(row).getIndex(position) != row)
                {
                    massMargin -= std::fabs(this->mass.getValue(row,position));
                }
            }
            for (uint position=0;position<this->stiffness.getNColumns(row);position++)
            {
                stiffnessRow += std::fabs(this->stiffness.getValue(row,position));
            }
            smallestMassMargin = std::min(smallestMassMargin,massMargin);
            largestStiffnessRow = std::max(largestStiffnessRow,stiffnessRow);
        }
        if (smallestMassMargin <= 0.0)
            throw std::invalid_argument("Central difference needs a positive strictly diagonally dominant mass matrix for a certified time-step limit; select lumped mass or an implicit method");
        this->centralDifferenceStepLimit = largestStiffnessRow > 0.0
            ? 2.0*std::sqrt(smallestMassMargin/largestStiffnessRow)
            : std::numeric_limits<double>::infinity();
    }
}

RRVector RStructuralDynamicsIntegrator::evaluateForce(double time) const
{
    RRVector result(this->mass.getNRows(),0.0);
    this->force(time,this->accepted,result);
    if (result.size() != this->mass.getNRows() || !finiteVector(result))
        throw std::invalid_argument("Force callback returned an invalid vector");
    return result;
}

const RSparseMatrix &RStructuralDynamicsIntegrator::operatorForStep(double dt)
{
    if (this->stepOperatorDt == dt)
    {
        return this->stepOperator;
    }
    RSparseMatrix assembled;
    assembled.setNRows(this->mass.getNRows());
    if (this->options.method == Method::CentralDifference)
    {
        addScaled(assembled,this->mass,1.0);
        addScaled(assembled,this->damping,0.5*dt);
    }
    else
    {
        double alphaM = 0.0, alphaF = 0.0;
        if (this->options.method == Method::GeneralizedAlpha)
        {
            const double rho = this->options.rhoInfinity;
            alphaM = (2.0*rho-1.0)/(rho+1.0);
            alphaF = rho/(rho+1.0);
        }
        const double gamma = 0.5-alphaM+alphaF;
        const double beta = 0.25*std::pow(1.0-alphaM+alphaF,2.0);
        addScaled(assembled,this->mass,1.0-alphaM);
        addScaled(assembled,this->damping,(1.0-alphaF)*gamma*dt);
        addScaled(assembled,this->stiffness,(1.0-alphaF)*beta*dt*dt);
    }
    this->stepOperator = std::move(assembled);
    this->stepOperatorDt = dt;
    return this->stepOperator;
}

RRVector RStructuralDynamicsIntegrator::solve(const RSparseMatrix &matrix, const RRVector &rhs) const
{
    RMatrixSolver linearSolver(this->solverConf);
    linearSolver.disableConvergenceLogFile();
    RRVector solution;
    linearSolver.solve(matrix,rhs,solution);
    if (!finiteVector(solution))
        throw std::runtime_error("Dynamics linear solve returned non-finite values");
    RRVector residual = product(matrix,solution);
    addScaled(residual,rhs,-1.0);
    if (RRVector::euclideanNorm(residual) >
        1.0e-8*RRVector::euclideanNorm(rhs)+1.0e-12)
        throw std::runtime_error("Dynamics linear solve did not converge");
    return solution;
}

void RStructuralDynamicsIntegrator::initialize(const RRVector &q0, const RRVector &v0, double time0)
{
    if (q0.size() != this->mass.getNRows() || v0.size() != this->mass.getNRows() ||
        !finiteVector(q0) || !finiteVector(v0) || !std::isfinite(time0))
        throw std::invalid_argument("Initial dynamics state has invalid dimensions or values");
    this->accepted = {};
    this->accepted.time = time0;
    this->accepted.q = q0;
    this->accepted.v = v0;
    this->accepted.externalForce = this->evaluateForce(time0);
    RRVector rhs = this->accepted.externalForce;
    addScaled(rhs,product(this->damping,v0),-1.0);
    addScaled(rhs,product(this->stiffness,q0),-1.0);
    this->accepted.a = this->solve(this->mass,rhs);
    this->accepted.kineticEnergy = quadraticEnergy(this->mass,v0);
    this->accepted.strainEnergy = quadraticEnergy(this->stiffness,q0);
    this->initialMechanicalEnergy = this->accepted.kineticEnergy+this->accepted.strainEnergy;
    this->fixedStep = 0.0;
    this->candidateReady = false;
    this->initialized = true;
}

void RStructuralDynamicsIntegrator::refreshAcceptedForce()
{
    if (!this->initialized) throw std::logic_error("Dynamics integrator is not initialized");
    if (this->candidateReady)
        throw std::logic_error("Cannot refresh accepted force while a trial step is pending");
    RRVector refreshed = this->evaluateForce(this->accepted.time);
    RRVector rhs = refreshed;
    addScaled(rhs,product(this->damping,this->accepted.v),-1.0);
    addScaled(rhs,product(this->stiffness,this->accepted.q),-1.0);
    RRVector acceleration = this->solve(this->mass,rhs);
    this->accepted.externalForce = std::move(refreshed);
    this->accepted.a = std::move(acceleration);
}

const RStructuralDynamicsState &RStructuralDynamicsIntegrator::trialStep(double dt)
{
    if (!this->initialized) throw std::logic_error("Dynamics integrator is not initialized");
    if (!std::isfinite(dt) || dt <= 0.0)
        throw std::invalid_argument("Dynamics time step must be positive and finite");
    if (this->options.method == Method::CentralDifference &&
        dt >= 0.95*this->centralDifferenceStepLimit)
        throw std::invalid_argument("Central difference time step exceeds the conservative stability limit");
    if (this->fixedStep > 0.0 && std::abs(dt-this->fixedStep) >
        1.0e-12*std::max(dt,this->fixedStep))
        throw std::invalid_argument("Dynamics integrator requires a fixed time step");

    // Never expose a partially evaluated candidate if a callback or solve throws.
    this->candidateReady = false;
    RStructuralDynamicsState next = this->accepted;
    next.time += dt;
    if (this->options.method == Method::CentralDifference)
    {
        // Central-difference half-step velocity; damping uses a centered term.
        RRVector previousHalf = this->accepted.v;
        addScaled(previousHalf,this->accepted.a,-0.5*dt);
        RRVector rhs = product(this->mass,previousHalf);
        addScaled(rhs,product(this->damping,previousHalf),-0.5*dt);
        addScaled(rhs,this->accepted.externalForce,dt);
        addScaled(rhs,product(this->stiffness,this->accepted.q),-dt);
        const RSparseMatrix &halfMatrix = this->operatorForStep(dt);
        RRVector nextHalf = this->solve(halfMatrix,rhs);
        next.q = this->accepted.q;
        addScaled(next.q,nextHalf,dt);
        next.externalForce = this->evaluateForce(next.time);
        rhs = next.externalForce;
        addScaled(rhs,product(this->stiffness,next.q),-1.0);
        addScaled(rhs,product(this->damping,nextHalf),-1.0);
        next.a = this->solve(halfMatrix,rhs);
        next.v = nextHalf;
        addScaled(next.v,next.a,0.5*dt);
    }
    else
    {
        double alphaM = 0.0, alphaF = 0.0;
        if (this->options.method == Method::GeneralizedAlpha)
        {
            const double rho = this->options.rhoInfinity;
            alphaM = (2.0*rho-1.0)/(rho+1.0);
            alphaF = rho/(rho+1.0);
        }
        const double gamma = 0.5-alphaM+alphaF;
        const double beta = 0.25*std::pow(1.0-alphaM+alphaF,2.0);
        RRVector qPredictor = this->accepted.q;
        addScaled(qPredictor,this->accepted.v,dt);
        addScaled(qPredictor,this->accepted.a,dt*dt*(0.5-beta));
        RRVector vPredictor = this->accepted.v;
        addScaled(vPredictor,this->accepted.a,dt*(1.0-gamma));

        // M[(1-am)a+am a_n]+C[(1-af)v+af v_n]+K[(1-af)q+af q_n]=f(t_n+(1-af)dt).
        RRVector effectiveForce = this->evaluateForce(this->accepted.time+(1.0-alphaF)*dt);
        RRVector rhs = effectiveForce;
        addScaled(rhs,product(this->mass,this->accepted.a),-alphaM);
        RRVector weightedVelocity = vPredictor;
        weightedVelocity *= 1.0-alphaF;
        addScaled(weightedVelocity,this->accepted.v,alphaF);
        addScaled(rhs,product(this->damping,weightedVelocity),-1.0);
        RRVector weightedDisplacement = qPredictor;
        weightedDisplacement *= 1.0-alphaF;
        addScaled(weightedDisplacement,this->accepted.q,alphaF);
        addScaled(rhs,product(this->stiffness,weightedDisplacement),-1.0);

        next.a = this->solve(this->operatorForStep(dt),rhs);
        next.q = qPredictor;
        addScaled(next.q,next.a,beta*dt*dt);
        next.v = vPredictor;
        addScaled(next.v,next.a,gamma*dt);
        next.externalForce = alphaF == 0.0 ? std::move(effectiveForce) :
                                              this->evaluateForce(next.time);
    }

    this->updateEnergy(next,dt);
    this->candidate = std::move(next);
    this->candidateReady = true;
    return this->candidate;
}

void RStructuralDynamicsIntegrator::updateEnergy(RStructuralDynamicsState &state, double dt) const
{
    state.kineticEnergy = quadraticEnergy(this->mass,state.v);
    state.strainEnergy = quadraticEnergy(this->stiffness,state.q);
    RRVector delta = state.q;
    addScaled(delta,this->accepted.q,-1.0);
    RRVector meanForce = this->accepted.externalForce;
    addScaled(meanForce,state.externalForce,1.0);
    state.externalWork += 0.5*RRVector::dot(meanForce,delta);
    const double initialDampingPower = RRVector::dot(this->accepted.v,product(this->damping,this->accepted.v));
    const double finalDampingPower = RRVector::dot(state.v,product(this->damping,state.v));
    state.dampingDissipation += 0.5*dt*(initialDampingPower+finalDampingPower);
    state.energyBalanceResidual = state.kineticEnergy+state.strainEnergy-
                                  this->initialMechanicalEnergy-state.externalWork+
                                  state.dampingDissipation;
}

void RStructuralDynamicsIntegrator::acceptStep()
{
    if (!this->candidateReady) throw std::logic_error("No candidate dynamics step to accept");
    if (this->fixedStep == 0.0) this->fixedStep = this->candidate.time-this->accepted.time;
    this->accepted = std::move(this->candidate);
    this->candidate = {};
    this->candidateReady = false;
}

void RStructuralDynamicsIntegrator::rejectStep()
{
    this->candidate = {};
    this->candidateReady = false;
}

const RStructuralDynamicsState &RStructuralDynamicsIntegrator::acceptedState() const
{
    if (!this->initialized) throw std::logic_error("Dynamics integrator is not initialized");
    return this->accepted;
}

const RStructuralDynamicsState &RStructuralDynamicsIntegrator::candidateState() const
{
    if (!this->candidateReady) throw std::logic_error("No candidate dynamics step");
    return this->candidate;
}

RStructuralDynamicsPowerSample RStructuralDynamicsIntegrator::acceptedPowerSample() const
{
    if (!this->initialized) throw std::logic_error("Dynamics integrator is not initialized");
    RStructuralDynamicsPowerSample sample;
    sample.time = this->accepted.time;
    sample.effectiveForcePower = RRVector::dot(this->accepted.externalForce,this->accepted.v);
    sample.viscousDissipationPower = RRVector::dot(
        this->accepted.v,product(this->damping,this->accepted.v));
    sample.cumulativeEffectiveWork = this->accepted.externalWork;
    sample.cumulativeDampingDissipation = this->accepted.dampingDissipation;
    sample.energyBalanceResidual = this->accepted.energyBalanceResidual;
    if (!std::isfinite(sample.effectiveForcePower) ||
        !std::isfinite(sample.viscousDissipationPower))
        throw std::runtime_error("Dynamics power sample is non-finite");
    return sample;
}

RStructuralDynamicsCheckpoint RStructuralDynamicsIntegrator::checkpoint() const
{
    if (!this->initialized) throw std::logic_error("Dynamics integrator is not initialized");
    return {this->accepted,this->fixedStep,this->initialMechanicalEnergy};
}

void RStructuralDynamicsIntegrator::restore(const RStructuralDynamicsCheckpoint &checkpoint)
{
    if (checkpoint.state.q.size() != this->mass.getNRows() ||
        checkpoint.state.v.size() != this->mass.getNRows() ||
        checkpoint.state.a.size() != this->mass.getNRows() ||
        checkpoint.state.externalForce.size() != this->mass.getNRows() ||
        !finiteVector(checkpoint.state.q) || !finiteVector(checkpoint.state.v) ||
        !finiteVector(checkpoint.state.a) || !finiteVector(checkpoint.state.externalForce) ||
        !std::isfinite(checkpoint.state.time) ||
        !std::isfinite(checkpoint.state.externalWork) ||
        !std::isfinite(checkpoint.state.dampingDissipation) ||
        !std::isfinite(checkpoint.initialMechanicalEnergy) ||
        !std::isfinite(checkpoint.fixedStep) || checkpoint.fixedStep < 0.0)
        throw std::invalid_argument("Invalid dynamics checkpoint");
    this->accepted = checkpoint.state;
    this->accepted.kineticEnergy = quadraticEnergy(this->mass,this->accepted.v);
    this->accepted.strainEnergy = quadraticEnergy(this->stiffness,this->accepted.q);
    this->accepted.energyBalanceResidual = this->accepted.kineticEnergy+
        this->accepted.strainEnergy-checkpoint.initialMechanicalEnergy-
        this->accepted.externalWork+this->accepted.dampingDissipation;
    this->fixedStep = checkpoint.fixedStep;
    this->initialMechanicalEnergy = checkpoint.initialMechanicalEnergy;
    this->candidate = {};
    this->candidateReady = false;
    this->initialized = true;
}
