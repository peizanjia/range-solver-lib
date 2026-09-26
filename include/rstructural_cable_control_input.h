#ifndef RSTRUCTURAL_CABLE_CONTROL_INPUT_H
#define RSTRUCTURAL_CABLE_CONTROL_INPUT_H

#include "rstructural_cable_kernel.h"

using RStructuralCableControlInputMatrix6x3 =
    std::array<std::array<double,3>,6>;

//! Endpoint whose position is treated as a cable control input.
enum class RStructuralCableControlInputEndpoint
{
    Node1,
    Node2
};

//! Linearization for one scalar cable parameter.
struct RStructuralCableScalarInputLinearization
{
    //! -d(internalForce)/d(input), in node 1 xyz then node 2 xyz order.
    RStructuralCableVector6 forceMap{};
    //! d(potential energy)/d(input).
    double potentialDerivative = 0.0;
};

//! Linearization for translation of one cable endpoint in global xyz.
struct RStructuralCableAnchorInputLinearization
{
    //! -d(internalForce)/d(anchor translation), rows are cable force DOFs
    //! and columns are anchor x, y, z translation components.
    RStructuralCableControlInputMatrix6x3 forceMap{};
    //! d(potential energy)/d(anchor translation), in global xyz.
    RStructuralCableVector3 potentialDerivative{};
};

//! Analytic input sensitivities for the tension-only engineering cable law.
//!
//! Positions and referenceLength are in m, initialTension in N, and
//! axialRigidity EA in N. The reported force maps use the control convention
//! B = -d(internalForce)/d(input). Scalar inputs are differentiated with the
//! other cable parameters fixed. Anchor translation moves exactly one current
//! endpoint in global coordinates while the other endpoint remains fixed.
//!
//! Strictly slack states return zero sensitivities. At trial tension exactly
//! zero, the tension-only law is not differentiable with respect to these
//! inputs, so each method throws std::domain_error.
class RStructuralCableControlInput
{
    public:
        //! Differentiate with respect to referenceLength, holding EA and
        //! initialTension fixed.
        static RStructuralCableScalarInputLinearization referenceLength(
            const RStructuralCableVector3 &node1Position,
            const RStructuralCableVector3 &node2Position,
            double referenceLength,
            double initialTension,
            double axialRigidity);

        //! Differentiate with respect to initialTension, holding reference
        //! length and EA fixed.
        static RStructuralCableScalarInputLinearization initialTension(
            const RStructuralCableVector3 &node1Position,
            const RStructuralCableVector3 &node2Position,
            double referenceLength,
            double initialTension,
            double axialRigidity);

        //! Differentiate with respect to global xyz translation of endpoint.
        static RStructuralCableAnchorInputLinearization anchorTranslation(
            const RStructuralCableVector3 &node1Position,
            const RStructuralCableVector3 &node2Position,
            double referenceLength,
            double initialTension,
            double axialRigidity,
            RStructuralCableControlInputEndpoint endpoint);
};

#endif // RSTRUCTURAL_CABLE_CONTROL_INPUT_H
