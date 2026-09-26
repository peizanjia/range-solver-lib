#ifndef RSTRUCTURAL_THERMO_MONOLITHIC_1D_FE_ADAPTER_H
#define RSTRUCTURAL_THERMO_MONOLITHIC_1D_FE_ADAPTER_H

#include <QHash>

#include "rstructural_thermo_monolithic_1d.h"

class RModel;

//! Explicit physical Kelvin-Voigt viscosity in Pa*s keyed by the current
//! zero-based RModel element position. These keys are not stable across mesh
//! reordering; callers must rebuild the map after a topology change.
using RStructuralThermoMonolithicFEViscosityByElement = QHash<unsigned int,double>;

//! Direct-call bridge from a bounded, straight-collinear RModel TRUSS1 mesh to the
//! constant-coefficient monolithic thermoelastic kernel.
//!
//! This adapter is not registered as an RSolver task implementation: callers
//! invoke advance explicitly. The model task tree is still required to declare
//! exactly one Heat -> Stress iteration with leaf tasks. Loads, prescribed
//! values and heat inputs are supplied as arrays in current RModel node order;
//! force values use N, prescribed displacements use m, heat powers use W, and
//! prescribed temperatures use absolute K. Force and prescribed-displacement
//! values are signed scalar components along the inferred positive axis. That
//! axis is the straight geometry direction whose dominant global component is
//! positive; ties are resolved X, then Y, then Z. Every array must have one
//! entry per RModel node. Model boundary and initial conditions are rejected.
//! Only neutral gravity and ambient-temperature environment entries are accepted.
//! Every element's material E, alpha, density, heat capacity and conductivity
//! must each be constant.
//! The mesh is limited to 64 nodes and 1024 members, and each element must
//! belong to exactly one RLine material group. Geometry must be straight and
//! collinear within max(1e-10 times the node span, 16 machine epsilons times
//! the largest absolute coordinate). Accepted displacement must have no
//! transverse component beyond max(1e-12 m, 1e-10 times its magnitude).
//!
//! A successful call publishes candidate nodal Kelvin temperature and global
//! three-component displacement along the inferred axis, plus six element
//! stress components and stress invariants. For consistency with the existing
//! line-element convention, R_VARIABLE_STRESS_X is the element-local axial
//! total stress including the Kelvin-Voigt contribution; it is not rotated into
//! global tensor components. Other local normal/shear components are zero and
//! the uniaxial von Mises value is abs(stressX). The returned result carries
//! per-step energy diagnostics; RModel has no dedicated thermoelastic energy
//! history field. The RModel time solver must be disabled, and computedTime
//! advances by options.timeStepSeconds.
class RStructuralThermoMonolithic1DFEAdapter
{
    public:
        //! Advance one direct macrostep. No caller-visible model field changes
        //! unless the kernel and all candidate output validation succeed.
        static RStructuralThermoMonolithic1DResult advance(
            RModel &model,
            const RStructuralThermoMonolithicFEViscosityByElement &viscosityPaSecondsByElement,
            const RStructuralThermoMonolithic1DConditions &conditions,
            const RStructuralThermoMonolithic1DOptions &options,
            double referenceTemperatureKelvin = 293.15);
};

#endif // RSTRUCTURAL_THERMO_MONOLITHIC_1D_FE_ADAPTER_H
