#ifndef RSTRUCTURAL_PIEZO_LAYERED_DYNAMICS_H
#define RSTRUCTURAL_PIEZO_LAYERED_DYNAMICS_H

#include "rstructural_piezo_layered_model_adapter.h"
#include "rstructural_dynamics_integrator.h"
#include <functional>

struct RStructuralPiezoLayeredDynamicsReadout
{
    QString portId,unit;
    RStructuralPiezoReadoutMode mode=RStructuralPiezoReadoutMode::None;
    bool actionRecord=false;
    double value=0.0,deltaValue=0.0;
};

//! Complete independent native state, also the checkpoint payload. Physical
//! vectors use sorted mesh node order, six global components per node. The
//! integrator coordinates are eta=zAbsolute-z0; initialZ is never an increment.
struct RStructuralPiezoLayeredDynamicsState
{
    QString sourceFingerprint,forcingIdentity;
    QByteArray driverIdentity,integrityDigest; //!< Compatibility/tamper seal, not authentication.
    //! state.dampingDissipation inside this checkpoint is the separate
    //! trapezoid endpoint-rate diagnostic, rather than midpoint discrete work.
    RStructuralDynamicsCheckpoint integrator;
    std::vector<double> heldActuatorVoltages; //!< Exact actuatorPortIds() order; sensor boundaries are immutable.
    std::vector<double> initialActuatorVoltages; //!< Initial absolute electrical-energy baseline.
    std::vector<double> perturbationForces; //!< Stateless callback's full physical endpoint force, excluding f0.
    std::vector<double> physicalDisplacements,physicalVelocities,physicalAccelerations,physicalReactions;
    RStructuralPiezoElectricalResponse electrical;
    std::vector<RStructuralPiezoLayeredDynamicsReadout> readouts;
    double kineticEnergy=0.0,storedEnergy=0.0,initialEnergy=0.0;
    double backgroundExternalWork=0.0,electricalWork=0.0,viscousDissipation=0.0; //!< Viscous work uses dt*vMean^T*C*vMean.
    double energyBalanceResidual=0.0,freeForceResidual=0.0;
};
using RStructuralPiezoLayeredDynamicsCheckpoint=RStructuralPiezoLayeredDynamicsState;

//! Actual linear layered FE plant, using the production Newmark average-
//! acceleration integrator and compiled consistent M/C/electrical tangent.
//! Value-returning operations leave caller data unchanged even on exceptions.
//! A local integrator owns each operation's scratch, so no pending candidate
//! or callback failure can consume an accepted state. Fixed affine supports
//! have zero velocity and do no work. No model publication/controller state.
class RStructuralPiezoLayeredDynamics
{
public:
    //! Must be stateless/deterministic for the same time, six global force/
    //! moment components per mesh node; excludes compiled backgroundForce.
    using BackgroundForce=std::function<std::vector<double>(double)>;
    RStructuralPiezoLayeredDynamics(const RStructuralPiezoLayeredModelAdapter &adapter,
        const RMatrixSolverConf &solverConf,double initialTime,
        BackgroundForce force,const QString &forcingIdentity);

    const std::vector<QString> &actuatorPortIds() const {return electricalSystem.condensed().actuatorPortIds;}
    RStructuralPiezoLayeredDynamicsState initialize() const;
    RStructuralPiezoLayeredDynamicsState initialize(const std::vector<double> &heldVa) const;
    RStructuralPiezoLayeredDynamicsState trial(const RStructuralPiezoLayeredDynamicsState &accepted,
        const std::vector<double> &heldVa,double dt) const;
    //! Original canonical grid, bound to this driver's initialTime. Uses the
    //! representable interval for physics and nominalDt for fixed-step identity.
    //! Optional intervalForce supplies this trial interval while the constructor
    //! force still certifies the preceding accepted boundary from its ledger.
    RStructuralPiezoLayeredDynamicsState trialOnGrid(const RStructuralPiezoLayeredDynamicsState &accepted,
        const std::vector<double> &heldVa,std::uint64_t stepIndex,double nominalDt,BackgroundForce intervalForce={}) const;
    RStructuralPiezoLayeredDynamicsCheckpoint checkpoint(const RStructuralPiezoLayeredDynamicsState &accepted) const;
    RStructuralPiezoLayeredDynamicsState restore(const RStructuralPiezoLayeredDynamicsCheckpoint &checkpoint) const;

private:
    RStructuralPiezoLayeredMeshOperators mesh;
    RStructuralPiezoLayeredCompiledData data;
    RStructuralPiezoElectricalSystem electricalSystem;
    RMatrixSolverConf solverConf;
    RSparseMatrix mass,damping,stiffness;
    double initialTime;
    BackgroundForce force;
    QString forcingIdentity;
    QByteArray identity;
    std::vector<double> z0,initialEta;
    double initialReducedEnergy=0.0;

    std::vector<double> forcing(double time) const;
    RStructuralPiezoLayeredDynamicsState trialInterval(const RStructuralPiezoLayeredDynamicsState &accepted,
        const std::vector<double> &heldVa,double dt,double rightTime,
        std::uint64_t gridIndex,double nominalDt,const BackgroundForce &intervalForce={}) const;
    RRVector rightHandSide(const std::vector<double> &physicalForce,const std::vector<double> &heldVa) const;
    RStructuralPiezoLayeredDynamicsState derive(const RStructuralDynamicsCheckpoint &checkpoint,
        const std::vector<double> &heldVa,const std::vector<double> &physicalForce) const;
    void validate(const RStructuralPiezoLayeredDynamicsState &state) const;
    QByteArray seal(const RStructuralPiezoLayeredDynamicsState &state) const;
};

#endif
