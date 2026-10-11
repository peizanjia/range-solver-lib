#ifndef RSTRUCTURAL_TRUSS_DYNAMICS_H
#define RSTRUCTURAL_TRUSS_DYNAMICS_H

#include <QByteArray>
#include <rml_structural_interaction.h>
#include "rstructural_truss_assembly.h"

struct RStructuralTrussDynamicsSystem
{
    std::vector<RStructuralTrussVector3> referenceNodes;
    std::vector<RStructuralTrussMember> members;
    std::vector<double> nodalLumpedMass; //!< One fixed isotropic mass per node.
    std::vector<bool> prescribedMask; //!< Full node-order xyz DOFs.
    //! Required full-size array. Free and curve-owned entries must be zero.
    std::vector<double> fixedPrescribedDisplacements;
    std::vector<RStructuralPrescribedMotion> prescribedMotions;
};

struct RStructuralTrussDynamicsInputs
{
    std::vector<RStructuralTrussInput> memberInputs;
    std::vector<double> externalForces;
};

struct RStructuralTrussDynamicsOptions
{
    std::size_t maxIterations = 50;
    std::size_t maxLineSearchIterations = 24;
    double absForceTol = 1.0e-9; //!< [0,1e-8], bounded by accepted-state certificate.
    double relForceTol = 1.0e-9; //!< [0,1e-8], bounded by accepted-state certificate.
    double linearRelativeTolerance = 1.0e-10;
};

struct RStructuralTrussDynamicsState
{
    double time = 0.0;
    std::vector<double> displacements, velocities, accelerations;
    RStructuralTrussDynamicsInputs appliedInputs;
    std::vector<double> reactions; //!< Full M*a+fint-fexternal, including free residual.
    double storedEnergy = 0.0;
    double kineticEnergy = 0.0;
    double cumulativeExternalWork = 0.0;
    double cumulativeSupportWork = 0.0;
    double cumulativeElectricalWork = 0.0;
    double cumulativeLengthWork = 0.0;
    double initialTotalEnergy = 0.0;
    QByteArray systemIdentity; //!< SHA256 compatibility binding; not authentication.
};

struct RStructuralTrussDynamicsWork
{
    double externalWork = 0.0, supportWork = 0.0;
    double electricalWork = 0.0, lengthWork = 0.0;
};

struct RStructuralTrussDynamicsResult
{
    RStructuralTrussDynamicsState candidate;
    RStructuralTrussDynamicsWork jumpWork, stepWork;
    double residualNorm = 0.0;
    std::size_t iterations = 0;
    double energyBalanceError = 0.0; //!< Cumulative physical T+U work diagnostic.
};

//! Pure, fixed-mass mixed TRUSS1 Newmark beta=1/4, gamma=1/2 driver.
//! Constructor copies and freezes geometry/laws/mass/supports/analytic curves.
//! Commands are held over each step: voltage, conductor charge, cable L0 and
//! external force. No time-varying L0 path exists in this API. V/Q switching,
//! support impulses, extrapolation, rotations, MPC and shared circuits are absent.
//! initialize requires prescribed q/v to agree with frozen analytic values.
//! refresh and trialStep revalidate the complete accepted left state and return
//! independent candidates; failures throw without modifying caller data.
class RStructuralTrussDynamics
{
public:
    explicit RStructuralTrussDynamics(const RStructuralTrussDynamicsSystem &system);
    RStructuralTrussDynamicsState initialize(double time,
        const std::vector<double> &displacements, const std::vector<double> &velocities,
        const RStructuralTrussDynamicsInputs &inputs) const;
    RStructuralTrussDynamicsResult refresh(const RStructuralTrussDynamicsState &accepted,
        const RStructuralTrussDynamicsInputs &newHeldInputs) const;
    RStructuralTrussDynamicsResult trialStep(const RStructuralTrussDynamicsState &accepted,
        const RStructuralTrussDynamicsInputs &newHeldInputs, double dt,
        const RStructuralTrussDynamicsOptions &options = {}) const;
    const QByteArray &systemIdentity() const { return identity; }

private:
    RStructuralTrussDynamicsSystem system;
    QByteArray identity;
    std::vector<std::size_t> freeDofs;
    void prescribed(double time, std::vector<double> &q,
                    std::vector<double> &v, std::vector<double> &a) const;
    void validateInputs(const RStructuralTrussDynamicsInputs &inputs) const;
    RStructuralTrussAssemblyResult validateState(const RStructuralTrussDynamicsState &state) const;
    void observe(RStructuralTrussDynamicsState &state,
                 const RStructuralTrussAssemblyResult &assembly) const;
};

#endif // RSTRUCTURAL_TRUSS_DYNAMICS_H
