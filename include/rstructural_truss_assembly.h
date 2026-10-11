#ifndef RSTRUCTURAL_TRUSS_ASSEMBLY_H
#define RSTRUCTURAL_TRUSS_ASSEMBLY_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>
#include <rml_sparse_matrix.h>

using RStructuralTrussVector3 = std::array<double,3>;
using RStructuralTrussVector6 = std::array<double,6>;
using RStructuralTrussMatrix6 = std::array<std::array<double,6>,6>;

struct RStructuralTrussPassiveLaw { double referenceLength = 0.0; double axialRigidity = 0.0; };
struct RStructuralTrussPiezoLaw
{
    double referenceLength = 0.0;
    double axialStiffness = 0.0; //!< N/m.
    double coupling = 0.0; //!< N/V, signed with respect to the independent electrode port.
    double capacitance = 0.0; //!< F, strictly positive.
};
struct RStructuralTrussCableLaw { double initialTension = 0.0; double axialRigidity = 0.0; };
using RStructuralTrussLaw = std::variant<RStructuralTrussPassiveLaw,
    RStructuralTrussPiezoLaw,RStructuralTrussCableLaw>;

struct RStructuralTrussNoInput {};
struct RStructuralTrussVoltage { double value = 0.0; };
struct RStructuralTrussCharge { double value = 0.0; }; //!< Conductor charge Q2 in coulombs.
struct RStructuralTrussCableLength { double value = 0.0; }; //!< Applied L0 in metres.
using RStructuralTrussInput = std::variant<RStructuralTrussNoInput,
    RStructuralTrussVoltage,RStructuralTrussCharge,RStructuralTrussCableLength>;

struct RStructuralTrussElectricalResult
{
    double voltage = 0.0;
    double conductorCharge = 0.0;
    double enthalpy = 0.0;
};

struct RStructuralTrussElementResult
{
    double length = 0.0;
    double axialForce = 0.0; //!< Signed tension; compression is negative.
    double axialTangent = 0.0;
    double storedEnergy = 0.0; //!< Physical U, including electric stored energy.
    double mechanicalPotential = 0.0; //!< H at fixed voltage, U at fixed charge.
    RStructuralTrussVector6 internalForce{};
    RStructuralTrussMatrix6 tangent{};
    std::optional<RStructuralTrussElectricalResult> electrical;
    std::optional<bool> cableTaut; //!< Absent for rods; false at the exact kink.
    std::optional<double> storedEnergyLengthDerivative; //!< Cable dU/dL0 at fixed EA/N0.
};

struct RStructuralTrussMember
{
    std::uint64_t physicalElementId = 0; //!< Nonzero, unique across all member laws.
    std::size_t node1 = 0;
    std::size_t node2 = 0;
    RStructuralTrussLaw law;
};

struct RStructuralTrussAssemblyResult
{
    std::vector<double> internalForces; //!< Full node-order xyz DOFs.
    RSparseMatrix tangent; //!< Full DOFs, including isolated nodes.
    double storedEnergy = 0.0;
    double mechanicalPotential = 0.0;
    std::vector<RStructuralTrussElementResult> memberResponses; //!< Original member order.
};

//! Pure mixed 3-D axial law and common sparse assembly, all quantities SI.
//! Caller supplies effective axial device k/g/Ce; this API does not reduce
//! arbitrary 3-D material tensors or claim lateral boundary equivalence.
//! Small material strain, updated chord direction, fixed topology. Electrodes
//! are independent ports: reversing mechanical endpoints never reverses V/Q2.
//! This incremental core is not standalone dynamics, a beam, a circuit solver,
//! or a complete control workbench. No mass, time integration or controller state.
//! Invalid/nonrepresentable input or result throws std::invalid_argument;
//! candidate construction never mutates caller data or returns partial results.
class RStructuralTrussAssembly
{
public:
    static RStructuralTrussElementResult evaluateElement(
        const RStructuralTrussVector3 &node1Position,
        const RStructuralTrussVector3 &node2Position,
        const RStructuralTrussLaw &law,
        const RStructuralTrussInput &input);

    static RStructuralTrussAssemblyResult evaluate(
        // Reference chord and displacement difference remain separate so tiny
        // extensions are independent of an arbitrary global coordinate origin.
        const std::vector<RStructuralTrussVector3> &referenceNodes,
        const std::vector<double> &displacements,
        const std::vector<RStructuralTrussMember> &members,
        const std::vector<RStructuralTrussInput> &inputs);
};

#endif // RSTRUCTURAL_TRUSS_ASSEMBLY_H
