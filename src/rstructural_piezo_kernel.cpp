#include "rstructural_piezo_kernel.h"

#include <cmath>
#include <stdexcept>

namespace
{
bool finite(const RStructuralPiezoKernelResult &result)
{
    if (!std::isfinite(result.enthalpy) || !std::isfinite(result.mechanicalForce) ||
        !std::isfinite(result.dielectricDisplacementFlux) || !std::isfinite(result.strain) ||
        !std::isfinite(result.electricField) || !std::isfinite(result.stress) ||
        !std::isfinite(result.electricDisplacement) ||
        !std::isfinite(result.voltageFixedStiffness) ||
        !std::isfinite(result.openCircuitStiffness))
        return false;
    for (const auto &row : result.hessian)
    {
        for (double value : row)
        {
            if (!std::isfinite(value)) return false;
        }
    }
    return true;
}
}

RStructuralPiezoKernelResult RStructuralPiezoKernel::evaluate(
    double elongation,
    double voltageDifference,
    double cE,
    double e,
    double epsS,
    double area,
    double length)
{
    if (!std::isfinite(elongation) || !std::isfinite(voltageDifference) ||
        !std::isfinite(cE) || cE <= 0.0 || !std::isfinite(e) ||
        !std::isfinite(epsS) || epsS <= 0.0 || !std::isfinite(area) || area <= 0.0 ||
        !std::isfinite(length) || length <= 0.0)
        throw std::invalid_argument("Invalid structural piezo rod inputs or material parameters");

    const double areaOverLength = area/length;
    if (!std::isfinite(areaOverLength) || areaOverLength <= 0.0)
        throw std::invalid_argument("Structural piezo A/L must be positive and finite");

    RStructuralPiezoKernelResult result;
    result.strain = elongation/length;
    result.electricField = -voltageDifference/length;
    result.stress = cE*result.strain-e*result.electricField;
    result.electricDisplacement = e*result.strain+epsS*result.electricField;
    result.mechanicalForce = area*result.stress;
    result.dielectricDisplacementFlux = area*result.electricDisplacement;
    result.enthalpy = 0.5*areaOverLength*
                      (cE*elongation*elongation+
                       2.0*e*elongation*voltageDifference-
                       epsS*voltageDifference*voltageDifference);

    result.hessian[0][0] = cE*areaOverLength;
    result.hessian[0][1] = e*areaOverLength;
    result.hessian[1][0] = e*areaOverLength;
    result.hessian[1][1] = -epsS*areaOverLength;
    result.voltageFixedStiffness = result.hessian[0][0];
    result.openCircuitStiffness = areaOverLength*(cE+e*(e/epsS));

    if (!finite(result))
        throw std::invalid_argument("Structural piezo enthalpy, response, or tangent is non-finite");
    return result;
}
