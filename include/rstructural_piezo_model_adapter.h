#ifndef RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
#define RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H

#include <QtGlobal>

#include <vector>

#include <rml_structural_piezo_setup.h>

class RModel;

//! Converged static TRUSS1 piezoelectric FE response.
struct RStructuralPiezoModelResult
{
    //! Full model translations, flattened xyz by node, in solver units.
    std::vector<double> displacements;
    //! Full model internal nodal forces, flattened xyz by node, in solver units.
    std::vector<double> internalForces;
    //! Signed member axial forces, aligned with sorted TRUSS1 members by element ID.
    std::vector<double> memberAxialForces;
    //! Right-electrode conductor free charge in coulombs, aligned with the same
    //! sorted members. V1 requires exactly one configured port per member.
    std::vector<double> memberRightCharges;
    //! Euclidean physical free-DOF residual norm in solver force units.
    //! In NodalPotential mode this is the mechanical residual in N.
    double residualNorm = 0.0;
    //! Sum of 1-D rod electric enthalpies in solver energy units.
    double enthalpy = 0.0;
    //! Electrical formulation used to produce this result. Publication rejects
    //! results whose formulation no longer matches the current model setup.
    RStructuralPiezoElectricalMode electricalMode =
        RStructuralPiezoElectricalMode::LegacyTerminalVoltage;
    //! Nodal electric potential in V; populated only in NodalPotential mode.
    std::vector<double> nodalPotentials;
    //! Physical conductor charge in C, including prescribed-potential
    //! electrode reactions; populated only in NodalPotential mode.
    std::vector<double> nodalConductorCharges;
    //! Fresh free-node conductor-charge residual norm in C.
    double electricResidualNorm = 0.0;
};

//! Validates, solves and publishes bounded static TRUSS1 piezoelectric rods.
//!
//! LegacyTerminalVoltage preserves the version 1.3.14 terminal voltage model.
//! NodalPotential uses coupled 1-D mechanical/electric nodal fields and voltage
//! or physical conductor-charge boundary conditions. Neither mode represents
//! dynamics, control, or beam/patch kinematics. Residual certificates do not
//! certify uniqueness or stability when the member network contains a mechanism.
class RStructuralPiezoModelAdapter
{
    public:
        //! Validate piezo setup binding and supported static truss topology.
        static void validate(const RModel &model);

        //! Solve without mutating model result variables.
        static RStructuralPiezoModelResult solve(const RModel &model);

        //! Clear prior structural/electric potential/charge outputs before a repeated solve.
        static void clearResults(RModel &model);

        //! Publish only a dimensionally consistent finite converged result.
        static void publish(RModel &model,
                            const RStructuralPiezoModelResult &result);
};

#endif // RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
