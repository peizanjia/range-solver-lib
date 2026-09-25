#ifndef RSTRUCTURAL_THERMO_MONOLITHIC_1D_H
#define RSTRUCTURAL_THERMO_MONOLITHIC_1D_H

#include <cstddef>
#include <vector>

//! Constant-property, small-strain TRUSS1 member for the monolithic 1-D step.
//!
//! Connectivity must point in increasing global x. All values use SI units.
struct RStructuralThermoMonolithic1DMember
{
    std::size_t node0 = 0;
    std::size_t node1 = 0;
    double areaSquareMetres = 0.0;
    double youngModulusPascals = 0.0;
    double expansionPerKelvin = 0.0;
    double densityKilogramsPerCubicMetre = 0.0;
    double specificHeatJoulesPerKilogramKelvin = 0.0;
    double conductivityWattsPerMetreKelvin = 0.0;
    double kelvinVoigtViscosityPascalSeconds = 0.0;
};

//! One-dimensional mesh and its single linearization/reference temperature.
struct RStructuralThermoMonolithic1DModel
{
    std::vector<double> nodeCoordinatesMetres;
    std::vector<RStructuralThermoMonolithic1DMember> members;
    double referenceTemperatureKelvin = 293.15;
};

//! Accepted nodal temperature and axial displacement.
struct RStructuralThermoMonolithic1DState
{
    std::vector<double> temperatureKelvin;
    std::vector<double> axialDisplacementMetres;
};

//! Nodal loads and Dirichlet values at the end of the trial step.
struct RStructuralThermoMonolithic1DConditions
{
    std::vector<double> mechanicalForceNewtons;
    std::vector<bool> displacementPrescribed;
    std::vector<double> prescribedDisplacementMetres;
    std::vector<double> heatPowerWatts;
    std::vector<bool> temperaturePrescribed;
    std::vector<double> prescribedTemperatureKelvin;
};

struct RStructuralThermoMonolithic1DOptions
{
    double timeStepSeconds = 0.0;
    //! Add physical Kelvin-Voigt dissipation to the heat equation.
    bool heatKelvinVoigtDissipation = false;
    unsigned int maximumNewtonIterations = 30;
    unsigned int maximumLineSearchSteps = 32;
    double mechanicalAbsoluteResidualToleranceNewtons = 1.0e-9;
    double mechanicalRelativeResidualTolerance = 1.0e-10;
    double thermalAbsoluteResidualToleranceWatts = 1.0e-9;
    double thermalRelativeResidualTolerance = 1.0e-10;
};

//! Candidate results with physical residuals freshly reassembled from the candidate.
struct RStructuralThermoMonolithic1DResult
{
    //! Candidate state; the caller commits only after this function returns.
    RStructuralThermoMonolithic1DState nextState;
    std::vector<double> memberAxialStrain;
    std::vector<double> memberAxialStressPascals;
    std::vector<double> memberAxialForceNewtons;
    std::vector<double> memberDissipatedEnergyJoules;
    double mechanicalResidualNormNewtons = 0.0;
    double thermalResidualNormWatts = 0.0;
    double thermalEnergyChangeJoules = 0.0;
    double heaterEnergyJoules = 0.0;
    double reversibleThermoelasticEnergyJoules = 0.0;
    double viscousDissipationEnergyJoules = 0.0;
    double dissipationHeatEnergyJoules = 0.0;
    //! Heat supplied to the model by prescribed-temperature boundaries (positive in).
    double thermalBoundaryEnergyJoules = 0.0;
    //! Certified balance including the heat supplied/extracted by prescribed-T nodes.
    double thermalEnergyBalanceResidualJoules = 0.0;
    unsigned int newtonIterations = 0;
};

//! Candidate-only monolithic finite-element step for a constant-coefficient rod.
//!
//! The kernel solves the coupled quasi-static displacement and backward-Euler
//! temperature residuals together. TRUSS1 uses linear shape functions, a
//! consistent heat-capacity matrix, and the standard conduction matrix. The
//! reversible heat term is Tref E alpha epsilon_dot; optional heat feedback is
//! the nonnegative Kelvin-Voigt power eta epsilon_dot^2. Temperature-dependent
//! E, alpha, density, heat capacity or conductivity are deliberately rejected
//! by this API: those require a thermodynamically consistent free-energy and
//! entropy model beyond this linearized constant-property law.
//! For each member, sigma = E(epsilon-alpha(Tmean-Tref))+eta*epsilon_dot;
//! the assembled residual is R_u = integral(B^T sigma A) - f and
//! R_T = C*(T-Tn)/dt + K*T + Tref*E*alpha*A*B*(u-un)/dt - Q - Qd.
//!
//! The nonsymmetric coupled Jacobian is solved by equilibrated dense pivoted
//! elimination, not conjugate gradients. The model is bounded to 64 nodes
//! (128 coupled unknowns). A result is returned only after fresh mechanical
//! [N] and thermal [W] residual checks and strictly positive absolute Kelvin;
//! acceptedState and all input objects remain unchanged on failure. Input is
//! also capped at 1024 members, 100 Newton steps and 64 line-search trials.
class RStructuralThermoMonolithic1D
{
    public:
        static RStructuralThermoMonolithic1DResult step(
            const RStructuralThermoMonolithic1DModel &model,
            const RStructuralThermoMonolithic1DState &acceptedState,
            const RStructuralThermoMonolithic1DConditions &conditions,
            const RStructuralThermoMonolithic1DOptions &options);
};

#endif // RSTRUCTURAL_THERMO_MONOLITHIC_1D_H
