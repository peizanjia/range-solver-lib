#ifndef RSTRUCTURAL_THERMO_COUPLED_1D_H
#define RSTRUCTURAL_THERMO_COUPLED_1D_H

#include <QString>

#include <functional>
#include <vector>

//! Accepted nodal state, keyed by caller-assigned stable physical node ID.
struct RStructuralThermoCoupled1DNodeState
{
    QString nodeId;
    double temperatureKelvin = 0.0;
    double axialDisplacement = 0.0;
};

//! Small-strain 1-D rod. Temperature-dependent properties use absolute kelvin.
struct RStructuralThermoCoupled1DMember
{
    QString elementId;
    QString node0Id;
    QString node1Id;
    double area = 0.0;                 //!< m^2
    double length = 0.0;               //!< m
    double conductivity = 0.0;         //!< W/(m K)
    double density = 0.0;              //!< kg/m^3
    double specificHeat = 0.0;         //!< J/(kg K)
    double internalViscosity = 0.0;    //!< Pa s; physical Kelvin-Voigt loss
    std::function<double(double)> youngModulusAtKelvin; //!< Pa
    std::function<double(double)> expansionAtKelvin;    //!< 1/K
};

//! 1-D physical mesh and constant thermal data.
struct RStructuralThermoCoupled1DModel
{
    std::vector<QString> nodeIds;
    std::vector<RStructuralThermoCoupled1DMember> members;
    double referenceTemperatureKelvin = 293.15;
};

struct RStructuralThermoCoupled1DForce
{
    QString nodeId;
    double forceNewtons = 0.0;
};

struct RStructuralThermoCoupled1DDisplacementBC
{
    QString nodeId;
    double displacementMetres = 0.0;
};

//! Heater command and its configured nonnegative power bound.
struct RStructuralThermoCoupled1DHeater
{
    QString nodeId;
    double powerWatts = 0.0;
    double maximumPowerWatts = 0.0;
};

struct RStructuralThermoCoupled1DConditions
{
    std::vector<RStructuralThermoCoupled1DForce> mechanicalForces;
    std::vector<RStructuralThermoCoupled1DDisplacementBC> prescribedDisplacements;
    std::vector<RStructuralThermoCoupled1DHeater> heaters;
};

struct RStructuralThermoCoupled1DOptions
{
    double dt = 0.0; //!< s
    //! Feed step-increment Kelvin-Voigt dissipation back into the heat equation.
    bool twoWayDissipationHeat = false;
    unsigned int maxCouplingIterations = 40;
    unsigned int maxLinearIterations = 2000;
    double thermalAbsoluteResidualTolerance = 1.0e-10; //!< W
    double thermalRelativeResidualTolerance = 1.0e-10;
    double mechanicalAbsoluteResidualTolerance = 1.0e-10; //!< N
    double mechanicalRelativeResidualTolerance = 1.0e-10;
    double energyAbsoluteTolerance = 1.0e-10; //!< J
    double energyRelativeTolerance = 1.0e-10;
    double couplingPowerAbsoluteTolerance = 1.0e-10; //!< W
    double couplingPowerRelativeTolerance = 1.0e-10;
};

struct RStructuralThermoCoupled1DNodeResult
{
    QString nodeId;
    double internalForceNewtons = 0.0;
    //! Required support force on prescribed-displacement nodes; zero otherwise.
    double supportReactionNewtons = 0.0;
};

struct RStructuralThermoCoupled1DResult
{
    //! Candidate state in model.nodeIds order. Caller commits only after return.
    std::vector<RStructuralThermoCoupled1DNodeState> nextState;
    std::vector<RStructuralThermoCoupled1DNodeResult> nodes;
    std::vector<double> memberAxialForcesNewtons;
    std::vector<double> memberElasticAxialForcesNewtons;
    std::vector<double> memberViscousAxialForcesNewtons;
    //! Dissipated work per member; must be nonnegative and is measured in J=N m.
    std::vector<double> memberDissipatedEnergyJoules;
    std::vector<double> memberDissipationPowerWatts;
    double thermalResidualNormWatts = 0.0;
    double mechanicalResidualNormNewtons = 0.0;
    double couplingPowerResidualWatts = 0.0;
    double thermalEnergyChangeJoules = 0.0;
    double heaterEnergyJoules = 0.0;
    double viscousDissipationEnergyJoules = 0.0;
    //! Dissipation work actually fed to the thermal equation (zero in one-way mode).
    double dissipationHeatEnergyJoules = 0.0;
    //! Delta thermal energy minus heater and physical viscous-dissipation work.
    double thermalEnergyBalanceResidualJoules = 0.0;
    //! Step work diagnostics. External loads are assumed constant over this increment.
    double mechanicalExternalWorkJoules = 0.0;
    double mechanicalSupportWorkJoules = 0.0;
    double elasticStrainEnergyChangeJoules = 0.0;
    //! Diagnostic only: reversible thermoelastic work is intentionally omitted.
    double mechanicalWorkBalanceResidualJoules = 0.0;
    unsigned int couplingIterations = 0;
};

//! Side-effect-free partitioned backward-Euler thermo-structural rod step.
//!
//! Thermal capacity is lumped by halves to each endpoint, the ends are
//! thermally insulated, and conductivity, density, and heat capacity are
//! constant per member. Mechanical equilibrium uses E(Tmean) and the
//! adaptively integrated free thermal strain integral(alpha(T) dT) from the model
//! reference temperature to the member mean temperature. Temperature fields
//! and displacements are joined by stable node IDs, never by vector position.
//! Material callbacks must be deterministic and side-effect-free.
//!
//! Optional two-way coupling iterates only physical Kelvin-Voigt dissipation
//! from the accepted-to-candidate displacement increment into the heat solve.
//! It is a partitioned quasi-static mechanical update, not transient structural
//! dynamics. The reversible thermoelastic heat term T0 beta:epsilon_dot is
//! omitted; monolithic coupling and beam/patch physics are outside this API.
//! The work outputs are diagnostics, not a global first-law certificate:
//! reversible thermoelastic work is not modeled. External force work assumes
//! the supplied nodal loads are constant during this increment. Support work
//! uses trapezoidal elastic reactions and backward-Euler viscous reactions.
//!
//! Returns a fully certified candidate or throws without modifying acceptedState.
//! The caller commits by replacing its accepted state only after this function
//! returns successfully.
class RStructuralThermoCoupled1D
{
    public:
        static RStructuralThermoCoupled1DResult step(
            const RStructuralThermoCoupled1DModel &model,
            const std::vector<RStructuralThermoCoupled1DNodeState> &acceptedState,
            const RStructuralThermoCoupled1DConditions &conditions,
            const RStructuralThermoCoupled1DOptions &options);
};

#endif // RSTRUCTURAL_THERMO_COUPLED_1D_H
