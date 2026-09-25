#ifndef RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
#define RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H

#include <QtGlobal>

#include <vector>

class RModel;

//! Converged static voltage-controlled 1-D terminal-rod FE response.
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
    double residualNorm = 0.0;
    //! Sum of 1-D rod electric enthalpies in solver energy units.
    double enthalpy = 0.0;
};

//! Validates, solves and publishes bounded static voltage-controlled TRUSS1 rods.
//!
//! Only the prescribed-voltage terminal rod model is represented. The adapter
//! does not add electric-potential unknowns, open-circuit constraints, dynamics,
//! control, beam/patch kinematics, or a coupled-field FE capability claim.
//! A force-residual certificate also does not certify uniqueness or stability
//! when the mechanical member network contains a mechanism.
class RStructuralPiezoModelAdapter
{
    public:
        //! Validate piezo setup binding and supported static truss topology.
        static void validate(const RModel &model);

        //! Solve without mutating model result variables.
        static RStructuralPiezoModelResult solve(const RModel &model);

        //! Clear prior structural/electric charge outputs before a repeated solve.
        static void clearResults(RModel &model);

        //! Publish only a dimensionally consistent finite converged result.
        static void publish(RModel &model,
                            const RStructuralPiezoModelResult &result);
};

#endif // RSTRUCTURAL_PIEZO_MODEL_ADAPTER_H
