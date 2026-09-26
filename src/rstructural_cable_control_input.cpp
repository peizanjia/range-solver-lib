#include "rstructural_cable_control_input.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace
{
RStructuralCableKernelResult evaluate(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double referenceLength,
    double initialTension,
    double axialRigidity)
{
    return RStructuralCableKernel::evaluate(
        node1Position,node2Position,referenceLength,initialTension,axialRigidity);
}

void requireDifferentiable(const RStructuralCableKernelResult &response)
{
    if (response.trialTension == 0.0)
        throw std::domain_error(
            "Structural cable input derivative is undefined at the tension kink");
}

RStructuralCableVector3 direction(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double length)
{
    RStructuralCableVector3 value{};
    for (std::size_t i=0;i<3;i++)
        value[i] = (node2Position[i]-node1Position[i])/length;
    return value;
}

void requireFinite(const RStructuralCableScalarInputLinearization &result)
{
    if (!std::isfinite(result.potentialDerivative))
        throw std::invalid_argument("Structural cable input derivative is non-finite");
    for (double value : result.forceMap)
    {
        if (!std::isfinite(value))
            throw std::invalid_argument("Structural cable force input map is non-finite");
    }
}

void requireFinite(const RStructuralCableAnchorInputLinearization &result)
{
    for (double value : result.potentialDerivative)
    {
        if (!std::isfinite(value))
            throw std::invalid_argument("Structural cable potential derivative is non-finite");
    }
    for (const auto &row : result.forceMap)
    {
        for (double value : row)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Structural cable force input map is non-finite");
        }
    }
}
}

RStructuralCableScalarInputLinearization RStructuralCableControlInput::referenceLength(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double referenceLength,
    double initialTension,
    double axialRigidity)
{
    const auto response = evaluate(
        node1Position,node2Position,referenceLength,initialTension,axialRigidity);
    requireDifferentiable(response);

    RStructuralCableScalarInputLinearization result;
    if (!response.taut) return result;

    const double axialStiffness = axialRigidity/referenceLength;
    const double tensionDerivative = -axialStiffness*(response.length/referenceLength);
    const auto cableDirection = direction(node1Position,node2Position,response.length);
    for (std::size_t i=0;i<3;i++)
    {
        result.forceMap[i] = cableDirection[i]*tensionDerivative;
        result.forceMap[i+3] = -cableDirection[i]*tensionDerivative;
    }
    result.potentialDerivative =
        (response.tension/axialStiffness)*tensionDerivative+
        (response.tension*response.tension)/(2.0*axialRigidity);

    requireFinite(result);
    return result;
}

RStructuralCableScalarInputLinearization RStructuralCableControlInput::initialTension(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double referenceLength,
    double initialTension,
    double axialRigidity)
{
    const auto response = evaluate(
        node1Position,node2Position,referenceLength,initialTension,axialRigidity);
    requireDifferentiable(response);

    RStructuralCableScalarInputLinearization result;
    if (!response.taut) return result;

    const double axialStiffness = axialRigidity/referenceLength;
    const auto cableDirection = direction(node1Position,node2Position,response.length);
    for (std::size_t i=0;i<3;i++)
    {
        result.forceMap[i] = cableDirection[i];
        result.forceMap[i+3] = -cableDirection[i];
    }
    result.potentialDerivative = response.tension/axialStiffness;

    requireFinite(result);
    return result;
}

RStructuralCableAnchorInputLinearization RStructuralCableControlInput::anchorTranslation(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double referenceLength,
    double initialTension,
    double axialRigidity,
    RStructuralCableControlInputEndpoint endpoint)
{
    if (endpoint != RStructuralCableControlInputEndpoint::Node1 &&
        endpoint != RStructuralCableControlInputEndpoint::Node2)
        throw std::invalid_argument("Invalid structural cable control endpoint");

    const auto response = evaluate(
        node1Position,node2Position,referenceLength,initialTension,axialRigidity);
    requireDifferentiable(response);

    RStructuralCableAnchorInputLinearization result;
    if (!response.taut) return result;

    const std::size_t firstColumn = endpoint == RStructuralCableControlInputEndpoint::Node1 ? 0 : 3;
    for (std::size_t row=0;row<6;row++)
    {
        for (std::size_t column=0;column<3;column++)
        {
            result.forceMap[row][column] = -response.tangent[row][firstColumn+column];
        }
    }
    const std::size_t firstForce = endpoint == RStructuralCableControlInputEndpoint::Node1 ? 0 : 3;
    for (std::size_t i=0;i<3;i++)
        result.potentialDerivative[i] = response.internalForce[firstForce+i];

    requireFinite(result);
    return result;
}
