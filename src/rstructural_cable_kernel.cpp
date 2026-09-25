#include "rstructural_cable_kernel.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace
{
bool finite(const RStructuralCableVector3 &vector)
{
    for (double value : vector)
    {
        if (!std::isfinite(value)) return false;
    }
    return true;
}

bool finite(const RStructuralCableKernelResult &result)
{
    if (!std::isfinite(result.length) || !std::isfinite(result.trialTension) ||
        !std::isfinite(result.tension) || !std::isfinite(result.energy))
        return false;
    for (double value : result.internalForce)
    {
        if (!std::isfinite(value)) return false;
    }
    for (const auto &row : result.tangent)
    {
        for (double value : row)
        {
            if (!std::isfinite(value)) return false;
        }
    }
    return true;
}
}

RStructuralCableKernelResult RStructuralCableKernel::evaluate(
    const RStructuralCableVector3 &node1Position,
    const RStructuralCableVector3 &node2Position,
    double referenceLength,
    double initialTension,
    double axialRigidity)
{
    if (!finite(node1Position) || !finite(node2Position) ||
        !std::isfinite(referenceLength) || referenceLength <= 0.0 ||
        !std::isfinite(initialTension) || initialTension < 0.0 ||
        !std::isfinite(axialRigidity) || axialRigidity <= 0.0)
        throw std::invalid_argument("Invalid structural cable geometry or material parameters");

    const double axialStiffness = axialRigidity/referenceLength;
    if (!std::isfinite(axialStiffness) || axialStiffness <= 0.0)
        throw std::invalid_argument("Structural cable EA/L0 must be positive and finite");

    RStructuralCableVector3 direction{};
    for (std::size_t i=0;i<3;i++)
    {
        direction[i] = node2Position[i]-node1Position[i];
        if (!std::isfinite(direction[i]))
            throw std::invalid_argument("Structural cable coordinate difference is non-finite");
    }
    const double length = std::hypot(std::hypot(direction[0],direction[1]),direction[2]);
    if (!std::isfinite(length) || length <= 0.0)
        throw std::invalid_argument("Structural cable current length must be positive and finite");
    for (double &component : direction) component /= length;

    const double extension = length-referenceLength;
    const double trialTension = initialTension+axialStiffness*extension;
    if (!std::isfinite(extension) || !std::isfinite(trialTension))
        throw std::invalid_argument("Structural cable trial tension is non-finite");

    RStructuralCableKernelResult result;
    result.length = length;
    result.trialTension = trialTension;
    result.taut = trialTension > 0.0;
    if (!result.taut) return result; // Includes the kink: choose the slack tangent.

    result.tension = trialTension;
    const long double energy = 0.5L*static_cast<long double>(result.tension)*
                               static_cast<long double>(result.tension)/
                               static_cast<long double>(axialStiffness);
    if (!std::isfinite(energy) || energy > std::numeric_limits<double>::max())
        throw std::invalid_argument("Structural cable energy is non-finite");
    result.energy = static_cast<double>(energy);

    const double geometricStiffness = result.tension/length;
    if (!std::isfinite(geometricStiffness))
        throw std::invalid_argument("Structural cable geometric stiffness is non-finite");

    std::array<std::array<double,3>,3> directionTangent{};
    for (std::size_t i=0;i<3;i++)
    {
        result.internalForce[i] = -result.tension*direction[i];
        result.internalForce[i+3] = result.tension*direction[i];
        for (std::size_t j=0;j<3;j++)
        {
            const double nn = direction[i]*direction[j];
            directionTangent[i][j] = axialStiffness*nn+
                                     geometricStiffness*((i == j ? 1.0 : 0.0)-nn);
        }
    }

    // Hessian block signs follow internal force [-N n, +N n].
    for (std::size_t i=0;i<3;i++)
    {
        for (std::size_t j=0;j<3;j++)
        {
            const double value = directionTangent[i][j];
            result.tangent[i][j] = value;
            result.tangent[i][j+3] = -value;
            result.tangent[i+3][j] = -value;
            result.tangent[i+3][j+3] = value;
        }
    }

    if (!finite(result))
        throw std::invalid_argument("Structural cable force or tangent is non-finite");
    return result;
}
