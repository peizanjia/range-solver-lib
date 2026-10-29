#ifndef RSTRUCTURAL_PIEZO_AXIAL_OPERATORS_H
#define RSTRUCTURAL_PIEZO_AXIAL_OPERATORS_H

#include "rstructural_piezo_native_readout.h"

class RModel;

struct RStructuralPiezoAxialDof
{
    uint physicalNodeId = 0;
    uint axis = 0; //!< 0/1/2 = global X/Y/Z translation [m].
    std::size_t fullIndex = 0;
    bool prescribed = false;
};

struct RStructuralPiezoAxialInput
{
    QString portId; //!< Empty for a physical nodal force input.
    QString unit; //!< V or N; the current classified native catalog rejects charge actuators.
    uint physicalNodeId = 0, axis = 0;
    std::size_t memberIndex = 0; //!< Meaningful only for an electrical input.
};

struct RStructuralPiezoAxialOutput
{
    QString portId, unit;
    RStructuralPiezoReadoutMode mode = RStructuralPiezoReadoutMode::None;
    bool actionRecord = false;
    double value0 = 0.0;
};

//! Dense row-major SI matrices. All port columns use orderedPortIds; B/D
//! columns use inputs, and L*/D rows use outputs. Incremental equation:
//! M*dqdd + C*dqd + K*dq = B*du, dy = Lq*dq+Lv*dqd+La*dqdd+D*du.
//! K = Kmaterial + Kgeometric + KfixedCharge; Kshort omits condensation.
//! Geometric stiffness includes the actual electrical/preload axial force.
struct RStructuralPiezoAxialOperatorBlock
{
    std::size_t dofs = 0, ports = 0, inputs = 0, outputs = 0;
    std::vector<double> M, C, K, Kshort, Kmaterial, Kgeometric, KfixedCharge;
    std::vector<double> G, Ce, B, Lq, Lv, La, D;
};

//! Actual frozen mixed TRUSS1/POINT native FE working-point package, not a
//! replacement plant. An explicit reference must have the native certificate;
//! an implicit zero-displacement source reference is independently checked.
//! Free coordinates are selection only; no rotation DOFs or shared circuits.
//! C is zero because the current native dynamics is physically undamped.
//! Designing a different C here does not install it into native time integration.
//! This factory prepares numerical design/export; runtime design integration is
//! pending. Moving supports, cable kinks and unstable/mechanism references fail.
//! Scales with nonrepresentable nonzero contributions or a locally unresolved
//! native tangent are rejected rather than exporting erased cross couplings.
//! Construction never modifies the source; source identity is compatibility data.
class RStructuralPiezoAxialOperators
{
public:
    static RStructuralPiezoAxialOperators compile(const RModel &model);

    QString sourceFingerprint;
    QString certificateAuthority;
    QByteArray nativeSystemIdentity, readoutIdentity;
    RStructuralPiezoNativeReference reference;
    RStructuralTrussDynamicsSystem nativeSystem;
    std::vector<RStructuralPiezoAxialDof> fullDofs, freeDofs;
    std::vector<std::size_t> freeToFull;
    std::vector<uint> orderedMemberPhysicalIds;
    std::vector<QString> orderedPortIds;
    std::vector<RStructuralPiezoAxialInput> inputs;
    std::vector<RStructuralPiezoAxialOutput> outputs;
    RStructuralPiezoAxialOperatorBlock full, free;
    std::vector<double> q0, f0, internalForce0, reactions0, V0, Q0;
    std::vector<double> initialDeltaQ, initialVelocity, initialAcceleration;
    std::vector<double> freeInitialDeltaQ, freeInitialVelocity, freeInitialAcceleration;
    double minimumNormalizedStiffnessPivot = 0.0;
};

#endif // RSTRUCTURAL_PIEZO_AXIAL_OPERATORS_H
