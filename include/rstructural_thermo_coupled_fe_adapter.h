#ifndef RSTRUCTURAL_THERMO_COUPLED_FE_ADAPTER_H
#define RSTRUCTURAL_THERMO_COUPLED_FE_ADAPTER_H

#include <QHash>

#include "rstructural_thermo_coupled_1d.h"

class RModel;

//! Explicit physical Kelvin-Voigt input for one TRUSS1 model element.
//!
//! This legacy entry point keys the coefficient by the element's current
//! zero-based model index. Use advanceByPhysicalId for persisted model IDs.
//! The adapter never infers this coefficient from a
//! fluid material's dynamic viscosity or from structural Rayleigh damping.
using RStructuralThermoFEViscosityByElement = QHash<unsigned int,double>;

//! Bounded bridge between an RModel and the certified 1-D thermo-mechanical
//! kernel. Supports only straight x-aligned TRUSS1 meshes with one RLine group
//! (area and material) per element. Density, heat capacity, and conductivity
//! must be temperature-independent; Young modulus and thermal expansion may
//! use the material's temperature tables. Loads, supports, and nodal heaters
//! are supplied explicitly as kernel conditions using the zero-based RModel
//! node index rendered as a decimal QString. The RModel time scheduler must be
//! disabled. options.dt is the physical step interval; each accepted call
//! advances RTimeSolver::computedTime by dt while preserving the disabled
//! scheduler's cursor. RModel result compTime is not available in the current
//! model-library implementation.
//!
//! A successful step atomically publishes nodal temperature and 3-component
//! displacement, plus all six element stress fields. R_VARIABLE_STRESS_X is
//! the axial total force per area, including Kelvin-Voigt force; the other
//! normal and shear components are zero under the supported uniaxial rod
//! assumption. This is quasi-static mechanics with backward-Euler thermal
//! capacity and physical viscous dissipation feedback; it is not a general
//! transient, beam, or full thermoelastic first-law solver.
class RStructuralThermoCoupledFEAdapter
{
    public:
        //! Advance one macro increment. All coefficients and the candidate
        //! result are validated before the caller's RModel is replaced.
        static RStructuralThermoCoupled1DResult advance(
            RModel &model,
            const RStructuralThermoFEViscosityByElement &viscosityPaSecondsByElement,
            const RStructuralThermoCoupled1DConditions &conditions,
            const RStructuralThermoCoupled1DOptions &options,
            double referenceTemperatureKelvin = 293.15);

        //! Accept persisted physical node/element IDs and return node records
        //! under those IDs. Missing, padded or duplicate IDs fail before solve.
        static RStructuralThermoCoupled1DResult advanceByPhysicalId(
            RModel &model,
            const RStructuralThermoFEViscosityByElement &viscosityPaSecondsByPhysicalElement,
            const RStructuralThermoCoupled1DConditions &conditionsByPhysicalNode,
            const RStructuralThermoCoupled1DOptions &options,
            double referenceTemperatureKelvin = 293.15);
};

#endif // RSTRUCTURAL_THERMO_COUPLED_FE_ADAPTER_H
