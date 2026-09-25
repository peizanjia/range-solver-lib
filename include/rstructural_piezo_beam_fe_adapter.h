#ifndef RSTRUCTURAL_PIEZO_BEAM_FE_ADAPTER_H
#define RSTRUCTURAL_PIEZO_BEAM_FE_ADAPTER_H

#include <vector>

#include <QString>
#include <QtGlobal>

class RModel;

//! Electrical boundary mode for independent, through-thickness beam patches.
enum class RStructuralPiezoBeamElectricalMode
{
    PrescribedVoltage,
    PrescribedCharge
};

//! One bonded axial piezo patch on a native BEAM1 element.
struct RStructuralPiezoBeamPatch
{
    uint elementId = 0;
    //! Patch cross-sectional area in m^2.
    double area = 0.0;
    //! Electrode separation in m.
    double thickness = 0.0;
    //! Patch centroid offset along local beam Y from the host neutral axis, in m.
    double centroidOffsetY = 0.0;
    //! Patch axial Young modulus in Pa.
    double youngModulus = 0.0;
    //! Signed axial piezoelectric coupling e31 in C/m^2.
    double coupling = 0.0;
    //! Constant-strain permittivity in F/m.
    double permittivity = 0.0;
    //! V=phi(+localY)-phi(-localY) in voltage mode; the conjugate terminal
    //! charge Q=C*V-g^T*u in charge mode. Units are V and C respectively.
    double prescribedValue = 0.0;
};

//! Bounded input for a small-strain, static planar beam/patch solve.
struct RStructuralPiezoBeamFEInput
{
    RStructuralPiezoBeamElectricalMode electricalMode =
        RStructuralPiezoBeamElectricalMode::PrescribedVoltage;
    //! Each BEAM1 patch is an electrically isolated electrode pair.
    std::vector<RStructuralPiezoBeamPatch> patches;
    //! Physical residual gates, in N and C.
    double mechanicalAbsoluteTolerance = 1.0e-8;
    double mechanicalRelativeTolerance = 1.0e-9;
    double electricAbsoluteTolerance = 1.0e-14;
    double electricRelativeTolerance = 1.0e-9;
};

//! Converged response in model-node / input-patch order.
struct RStructuralPiezoBeamFEResult
{
    //! Global translations, flattened xyz by node, in m.
    std::vector<double> displacements;
    //! Global rotations, flattened xyz by node, in rad.
    std::vector<double> rotations;
    //! Recovered full six-component nodal internal resultants, N and N m.
    std::vector<double> internalForces;
    //! Patch electrode potential differences, in V.
    std::vector<double> patchVoltages;
    //! Terminal charges under the declared electric-enthalpy sign convention, in C.
    std::vector<double> patchCharges;
    //! Total mechanical/electric enthalpy, in J.
    double enthalpy = 0.0;
    //! Fresh norm of free in-plane nodal force/moment residuals.
    double mechanicalResidualNorm = 0.0;
    //! Fresh norm of free patch-charge residuals; zero for voltage control.
    double electricResidualNorm = 0.0;
    //! Exact model/input binding checked again at publication time.
    QString modelFingerprint;
};

//! Direct-call, multi-element BEAM1 piezoelectric patch finite-element adapter.
//!
//! The adapter reads host geometry, E, A and Iz from the RModel and adds the
//! axial patch layer and electromechanical coupling to a planar XY frame.
//! It currently accepts straight +X BEAM1 chains with local Y aligned to +Y,
//! point-group mechanical conditions, zero gravity on patched members, and
//! electrically isolated patch pairs. In-plane line gravity is assembled only
//! for unpatched host members; out-of-plane gravity is rejected.
//! It does not wire this formulation into RSolver or RStructuralPiezoSetup.
class RStructuralPiezoBeamFEAdapter
{
    public:
        static void validate(const RModel &model,
                             const RStructuralPiezoBeamFEInput &input);

        //! Solve without mutating model result variables.
        static RStructuralPiezoBeamFEResult solve(
            const RModel &model,
            const RStructuralPiezoBeamFEInput &input);

        //! Publish a dimensionally consistent, finite converged candidate.
        static void publish(RModel &model,
                           const RStructuralPiezoBeamFEInput &input,
                           const RStructuralPiezoBeamFEResult &result);
};

#endif // RSTRUCTURAL_PIEZO_BEAM_FE_ADAPTER_H
