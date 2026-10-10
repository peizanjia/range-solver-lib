#ifndef RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
#define RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H

#include <QtGlobal>
#include <QString>

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
    //! sorted members. Passive member entries are zero placeholders, not electrical
    //! measurements; consult memberElectricalPortAvailable before using them.
    std::vector<double> memberRightCharges;
    //! True exactly for members bound to an active electrical port, in sorted
    //! element-ID order. Publication verifies this identity against the setup.
    std::vector<bool> memberElectricalPortAvailable;
    //! Euclidean physical free-DOF residual norm in solver force units.
    //! In NodalPotential mode this is the mechanical residual in N.
    double residualNorm = 0.0;
    //! Sum of active rod electric enthalpies and passive rod strain energies
    //! in solver energy units.
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

    private:
        //! Adapter-owned consistency certificate, not trusted authentication.
        //! Hand-built or modified public results cannot impersonate solved data.
        QString inputFingerprint;
        QString resultSeal;
        friend class RStructuralPiezoModelAdapter;
};

//! Validates, solves and publishes bounded static TRUSS1 piezoelectric rods.
//!
//! LegacyTerminalVoltage supports explicitly classified active and passive 3-D
//! trusses. Empty passive classification retains the all-active legacy model.
//! Element charge publication retains full-model zero placeholders: only active
//! port elements are electrical outputs. This is a static foundation only.
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

        //! Publish only an unmodified adapter-owned converged result whose exact
        //! physical inputs still match this model. Repeated publication is legal.
        //! The seal enforces local adapter consistency, not trusted authentication.
        static void publish(RModel &model,
                            const RStructuralPiezoModelResult &result);
};

#endif // RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
