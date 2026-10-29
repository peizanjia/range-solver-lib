#ifndef RSTRUCTURAL_PIEZO_LAYERED_MODEL_ADAPTER_H
#define RSTRUCTURAL_PIEZO_LAYERED_MODEL_ADAPTER_H
#include "rstructural_piezo_layered_mesh.h"
#include "rstructural_piezo_electrical_system.h"
#include <rml_structural_piezo_layered_setup.h>
#include <memory>

struct RStructuralPiezoLayeredCompiledData
{
    uint fullDofs=0, freeDofs=0, constraintRank=0, redundantConstraints=0;
    //! q = T*z+p; row-major fullDofs by freeDofs, Euclidean orthonormal columns.
    std::vector<double> T,p;
    //! Complete consistent generalized operators; M/K/G/Ce/P are in mesh().
    std::vector<double> fullDamping, mass, stiffness, damping, input;
    //! Absolute certified reference and absolute initial coordinates, SI.
    std::vector<double> q0, actuatorVoltages0, referenceForce, backgroundForce;
    std::vector<double> initialZ, initialVelocity, initialAcceleration;
    //! Incremental RHS = backgroundOffset + T^T*f(t) + B*deltaVa.
    //! f(t) here excludes the separately stored source backgroundForce.
    std::vector<double> backgroundOffset, supportReactions;
    double forceResidual=0, electricResidual=0;
    QString sourceFingerprint, dampingAuthority;
    std::vector<RStructuralPiezoLayeredMaterial> materialProvenance;
    std::vector<RStructuralPiezoLayeredPort> portMetadata;
    //! Ordered physical node components; translations m/N, rotations rad/N m.
    std::vector<QString> generalizedCoordinateUnits, generalizedForceUnits;
};

//! Dense native source compiler. Homogeneous planar BEAM1 or flat QUAD1 only;
//! no host material is re-added. The production mesh and electrical system own
//! every FE operator. Static reference requires positive definite free tangent;
//! mechanisms and unstable references are rejected, including with explicit q0.
//! Layered setup is the single support/load/initial-state authority: populated
//! legacy group BC/IC/environment declarations are rejected until migrated.
//! Dense first compiler is bounded to 2048 supported generalized DOFs. This is
//! an admission bound, not a large-model performance qualification.
class RStructuralPiezoLayeredModelAdapter
{
public:
    //! authoritativeFullC is required only for modal strategies and uses the
    //! sorted physical-node supported generalized coordinate order (3 or 5/node).
    RStructuralPiezoLayeredModelAdapter(const RModel &model,
        const RStructuralPiezoLayeredSetup &setup,
        const std::vector<double> &authoritativeFullC={});
    const RStructuralPiezoLayeredMeshOperators &mesh() const {return meshOperators;}
    const RStructuralPiezoElectricalSystem &electrical() const {return *electricalSystem;}
    const RStructuralPiezoLayeredCompiledData &data() const {return compiled;}
    std::vector<double> reduceForce(const std::vector<double> &physical) const;
    std::vector<double> reconstruct(const std::vector<double> &z) const;
    std::vector<double> reduceDisplacement(const std::vector<double> &physical) const;
private:
    RStructuralPiezoLayeredMeshOperators meshOperators;
    RStructuralPiezoLayeredCompiledData compiled;
    std::unique_ptr<RStructuralPiezoElectricalSystem> electricalSystem;
};
#endif
