#ifndef RSTRUCTURAL_PIEZO_NATIVE_READOUT_H
#define RSTRUCTURAL_PIEZO_NATIVE_READOUT_H

#include <rml_structural_piezo_readout_setup.h>
#include "rstructural_truss_dynamics_model_adapter.h"

struct RStructuralPiezoNativeElectricalResidual
{
    QString portId;
    double voltageResidual = 0.0, chargeResidual = 0.0;
    double voltageScale = 0.0, chargeScale = 0.0;
};

struct RStructuralPiezoNativeReference
{
    bool explicitWorkpoint = false;
    bool equilibriumCertified = false;
    double time = 0.0;
    std::vector<double> displacements;
    RStructuralTrussDynamicsInputs inputs;
    RStructuralTrussAssemblyResult assembly; //!< Actual full tangent, including geometric terms and sensor condensation.
    double freeForceResidualNorm = 0.0, forceScale = 0.0, forceTolerance = 0.0;
    std::vector<RStructuralPiezoNativeElectricalResidual> electricalResiduals;
};

struct RStructuralPiezoNativeReadoutSample
{
    QString portId;
    RStructuralPiezoReadoutMode mode = RStructuralPiezoReadoutMode::None;
    QString unit;
    bool actionRecord = false; //!< None echoes the held voltage action, not a deformation measurement.
    double voltage = 0.0, charge = 0.0;
    double deltaVoltage = 0.0, deltaCharge = 0.0;
    double value = 0.0, deltaValue = 0.0;
};

//! Frozen native FE readout compiler. No results field supplies its reference.
//! Observation validates a complete native state and never changes it. The
//! initial FE signal reference is not called an equilibrium without an explicit
//! workpoint and successful force/electrical/prescribed-position certificates.
class RStructuralPiezoNativeReadout
{
public:
    RStructuralPiezoNativeReadout(const RStructuralTrussDynamicsModelAdapter &adapter,
                                 const RStructuralPiezoReadoutSetup &setup);
    //! Approve only a physical C change to the frozen native plant. All other
    //! geometry, laws, masses and supports must retain the adapter identity.
    RStructuralPiezoNativeReadout(const RStructuralTrussDynamicsModelAdapter &adapter,
                                 const RStructuralPiezoReadoutSetup &setup,
                                 const RStructuralTrussDynamicsSystem &approvedSystem);
    const RStructuralPiezoNativeReference &reference() const { return frozenReference; }
    const RStructuralPiezoReadoutSetup &setup() const { return catalog; }
    const std::map<QString,std::size_t> &portIndices() const { return readoutByPortId; }
    std::vector<RStructuralPiezoNativeReadoutSample> observe(const RStructuralTrussDynamicsState &state) const;
    //! Integration fast path: caller has validated state through driver.refresh
    //! and freshly evaluated this exact state with this frozen system. The
    //! supplied assembly is trusted; sensor boundaries are still checked here.
    std::vector<RStructuralPiezoNativeReadoutSample> observe(const RStructuralTrussDynamicsState &state,
                                                            const RStructuralTrussAssemblyResult &freshAssembly) const;

private:
    RStructuralTrussDynamicsSystem system;
    RStructuralTrussDynamics driver;
    RStructuralPiezoReadoutSetup catalog;
    std::vector<std::size_t> memberIndices;
    std::map<QString,std::size_t> readoutByPortId;
    RStructuralPiezoNativeReference frozenReference;
    void validateBoundaries(const RStructuralTrussDynamicsInputs &inputs) const;
};

#endif
