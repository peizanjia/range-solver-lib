#ifndef RSTRUCTURAL_DYNAMICS_INTEGRATOR_H
#define RSTRUCTURAL_DYNAMICS_INTEGRATOR_H

#include <functional>

#include <rml_matrix_solver_conf.h>
#include <rml_sparse_matrix.h>

//! Reduced free-DOF state for M a + C v + K q = f(t).
struct RStructuralDynamicsState
{
    double time = 0.0;
    RRVector q;
    RRVector v;
    RRVector a;
    RRVector externalForce;
    double kineticEnergy = 0.0;
    double strainEnergy = 0.0;
    double externalWork = 0.0;
    double dampingDissipation = 0.0;
    //! E(t)-E(0)-W_external+D_damping. Generalized-alpha can dissipate numerically.
    double energyBalanceResidual = 0.0;
};

//! Fully committed state; restoring it discards any pending trial step.
struct RStructuralDynamicsCheckpoint
{
    RStructuralDynamicsState state;
    double fixedStep = 0.0;
    double initialMechanicalEnergy = 0.0;
};

//! Read-only power ledger for the accepted reduced free-DOF state.
//! Both powers are endpoint rates, not step-average work/dt. Effective work
//! accumulates 0.5*(f_n+f_n+1) dot (q_n+1-q_n); damping dissipation uses a
//! trapezoid of endpoint v^T*C*v rates. Units follow the supplied matrices,
//! coordinates and time (the FE caller currently uses scaled solver units).
//! The force includes support-equivalent contributions, so neither rate is a
//! controller-specific or whole moving-support physical power measurement.
struct RStructuralDynamicsPowerSample
{
    double time = 0.0;
    double effectiveForcePower = 0.0;
    double viscousDissipationPower = 0.0;
    double cumulativeEffectiveWork = 0.0;
    double cumulativeDampingDissipation = 0.0;
    double energyBalanceResidual = 0.0;
};

class RStructuralDynamicsIntegrator
{
    public:
        //! Central difference is conditionally stable; the caller must choose dt below its CFL limit.
        enum class Method { NewmarkAverageAcceleration, GeneralizedAlpha, CentralDifference };

        struct Options
        {
            Method method = Method::NewmarkAverageAcceleration;
            //! High-frequency spectral radius in [0,1]; 1 removes algorithmic damping.
            double rhoInfinity = 1.0;
        };

        //! Must be deterministic for the same accepted state and time. The caller owns any controller state.
        using ForceCallback = std::function<void(double time,
                                                 const RStructuralDynamicsState &accepted,
                                                 RRVector &force)>;

        RStructuralDynamicsIntegrator(const RSparseMatrix &mass,
                                      const RSparseMatrix &damping,
                                      const RSparseMatrix &stiffness,
                                      const RMatrixSolverConf &solverConf,
                                      Options options,
                                      ForceCallback force);

        //! Solves M a0=f(t0)-C v0-K q0.
        void initialize(const RRVector &q0, const RRVector &v0, double time0 = 0.0);

        //! Re-evaluate a discontinuous, held control load at the accepted
        //! boundary. Call before trialStep() when the controller changes its
        //! command; failure leaves the accepted state unchanged. The caller
        //! must restore its outer transaction if the trial is rejected.
        void refreshAcceptedForce();

        //! Computes a candidate without changing the committed state. Once accepted, dt is fixed.
        const RStructuralDynamicsState &trialStep(double dt);
        void acceptStep();
        void rejectStep();

        bool hasCandidate() const { return candidateReady; }
        const RStructuralDynamicsState &acceptedState() const;
        const RStructuralDynamicsState &candidateState() const;
        //! Derived from accepted state and matrices; trial/reject never publishes it.
        RStructuralDynamicsPowerSample acceptedPowerSample() const;
        RStructuralDynamicsCheckpoint checkpoint() const;
        void restore(const RStructuralDynamicsCheckpoint &checkpoint);

    private:
        RSparseMatrix mass;
        RSparseMatrix damping;
        RSparseMatrix stiffness;
        RMatrixSolverConf solverConf;
        Options options;
        ForceCallback force;
        RStructuralDynamicsState accepted;
        RStructuralDynamicsState candidate;
        double fixedStep = 0.0;
        double initialMechanicalEnergy = 0.0;
        //! Certified conservative limit from a strictly diagonally dominant mass matrix.
        double centralDifferenceStepLimit = 0.0;
        //! Fixed-operator dynamics assembles the dt-dependent sparse system once.
        //! Force/controller changes affect only the RHS, never this matrix.
        RSparseMatrix stepOperator;
        double stepOperatorDt = 0.0;
        bool initialized = false;
        bool candidateReady = false;

        RRVector evaluateForce(double time) const;
        const RSparseMatrix &operatorForStep(double dt);
        RRVector solve(const RSparseMatrix &matrix, const RRVector &rhs) const;
        void updateEnergy(RStructuralDynamicsState &state, double dt) const;
};

#endif // RSTRUCTURAL_DYNAMICS_INTEGRATOR_H
