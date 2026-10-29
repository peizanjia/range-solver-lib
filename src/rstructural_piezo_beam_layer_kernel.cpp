#include "rstructural_piezo_element_operators.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace
{
using Row = std::array<long double,6>;

// Keep the binary scale separate until every physical factor is included.
// MSVC long double has double's exponent range: density*width or the volume
// alone can underflow although the final mass/stiffness is representable.
long double scaledProduct(std::initializer_list<long double> factors,
                          std::initializer_list<long double> divisors={})
{
    bool zero = false;
    for (long double factor : factors)
    {
        if (!std::isfinite(factor))
            throw std::invalid_argument("Layered beam operators exceed the finite numeric range");
        zero = zero || factor == 0.0L;
    }
    for (long double divisor : divisors)
        if (!std::isfinite(divisor) || divisor == 0.0L)
            throw std::invalid_argument("Layered beam operators exceed the finite numeric range");
    if (zero) return 0.0L;
    long double mantissa = 1.0L;
    int exponent = 0;
    for (long double factor : factors)
    {
        int factorExponent = 0;
        mantissa *= std::frexp(factor,&factorExponent);
        exponent += factorExponent;
    }
    for (long double divisor : divisors)
    {
        int divisorExponent = 0;
        mantissa /= std::frexp(divisor,&divisorExponent);
        exponent -= divisorExponent;
    }
    return std::scalbn(mantissa,exponent);
}

double checkedDouble(long double value)
{
    const double result = static_cast<double>(value);
    if (!std::isfinite(value) || !std::isfinite(result))
        throw std::invalid_argument("Layered beam operators exceed the finite numeric range");
    return result;
}

void copyOperators(const std::vector<long double> &source,std::vector<double> &target)
{
    target.reserve(source.size());
    for (long double value : source) target.push_back(checkedDouble(value));
}
}

RStructuralPiezoElementOperators RStructuralPiezoBeamLayerKernel::assemble(
    double length,const std::vector<RStructuralPiezoLayer> &layers)
{
    if (!std::isfinite(length) || length <= 0.0 || layers.empty())
        throw std::invalid_argument("Layered beam requires a positive finite length and layers");

    RStructuralPiezoElementOperators result;
    result.dofs = 6;
    for (size_t i=0;i<layers.size();++i)
    {
        const auto &layer = layers[i];
        layer.validate();
        if (layer.coverage[2] != -1.0 || layer.coverage[3] != 1.0)
            throw std::invalid_argument("Beam coverage must use the complete default eta interval");
        for (size_t j=0;j<i;++j)
        {
            const auto &other = layers[j];
            // Widths are centered on the same local y axis, so any positive
            // intersection in x and z would count material volume twice.
            if (std::max(layer.coverage[0],other.coverage[0]) <
                    std::min(layer.coverage[1],other.coverage[1]) &&
                std::max(layer.lowerZ,other.lowerZ) < std::min(layer.upperZ,other.upperZ))
                throw std::invalid_argument("Layered beam material volumes overlap");
        }
        if (!layer.portId.isEmpty() &&
            std::find(result.portIds.begin(),result.portIds.end(),layer.portId) == result.portIds.end())
            result.portIds.push_back(layer.portId);
    }

    const size_t ports = result.portIds.size();
    std::vector<long double> mass(36,0.0L),stiffness(36,0.0L);
    std::vector<long double> coupling(6*ports,0.0L),capacitance(ports*ports,0.0L);
    long double totalMass = 0.0L;
    const long double beamLength = length;
    // Four axial points integrate the degree-six Hermite mass integrand
    // exactly. Two thickness points retain both eccentric and rotary inertia.
    constexpr std::array<long double,4> points{{
        -0.86113631159405257522394648889281L,-0.33998104358485626480266575910324L,
         0.33998104358485626480266575910324L, 0.86113631159405257522394648889281L}};
    constexpr std::array<long double,4> weights{{
        0.347854845137453857373063949221999L,0.652145154862546142626936050778001L,
        0.652145154862546142626936050778001L,0.347854845137453857373063949221999L}};
    constexpr long double zPoint = 0.57735026918962576450914878050196L;

    for (const auto &layer : layers)
    {
        const auto reduced = rStructuralPiezoReduceFreeStress(layer.material,{0});
        // Half-endpoint sums avoid overflow when the offset is large. Do not
        // impose a geometric epsilon: small representable layers are valid.
        const long double zCenter = static_cast<long double>(layer.lowerZ)/2.0L+
                                    static_cast<long double>(layer.upperZ)/2.0L;
        const long double zHalf = static_cast<long double>(layer.upperZ)/2.0L-
                                  static_cast<long double>(layer.lowerZ)/2.0L;
        const long double sLower = (static_cast<long double>(layer.coverage[0])+1.0L)/2.0L;
        const long double sUpper = (static_cast<long double>(layer.coverage[1])+1.0L)/2.0L;
        const long double sCenter = sLower/2.0L+sUpper/2.0L;
        const long double sHalf = (sUpper-sLower)/2.0L;
        const long double width = layer.width;
        const long double density = layer.material.density;
        totalMass += scaledProduct({density,width,zHalf,beamLength,sUpper-sLower,2.0L});
        const bool hasPort = !layer.portId.isEmpty();
        const size_t port = hasPort ? static_cast<size_t>(std::find(
            result.portIds.begin(),result.portIds.end(),layer.portId)-result.portIds.begin()) : 0;
        if (hasPort)
            capacitance[port*ports+port] += scaledProduct(
                {reduced.permittivity,width,beamLength,sUpper-sLower},{zHalf,2.0L});

        for (size_t ix=0;ix<points.size();++ix)
        {
            const long double s = sCenter+sHalf*points[ix];
            const long double s2 = s*s,s3 = s2*s;
            const Row axial{{1.0L-s,0.0L,0.0L,s,0.0L,0.0L}};
            const Row transverse{{0.0L,1.0L-3.0L*s2+2.0L*s3,
                -beamLength*(s-2.0L*s2+s3),0.0L,3.0L*s2-2.0L*s3,
                -beamLength*(-s2+s3)}};
            const Row rotation{{0.0L,(6.0L*s-6.0L*s2)/beamLength,
                1.0L-4.0L*s+3.0L*s2,0.0L,(-6.0L*s+6.0L*s2)/beamLength,
                -2.0L*s+3.0L*s2}};
            const Row axialDerivative{{-1.0L/beamLength,0.0L,0.0L,
                1.0L/beamLength,0.0L,0.0L}};
            const Row rotationDerivative{{0.0L,((6.0L-12.0L*s)/beamLength)/beamLength,
                (-4.0L+6.0L*s)/beamLength,0.0L,
                ((-6.0L+12.0L*s)/beamLength)/beamLength,(-2.0L+6.0L*s)/beamLength}};
            for (int iz : {-1,1})
            {
                const long double z = zCenter+iz*zHalf*zPoint;
                Row longitudinal{},strain{};
                for (size_t i=0;i<6;++i)
                {
                    longitudinal[i] = axial[i]+z*rotation[i];
                    strain[i] = axialDerivative[i]+z*rotationDerivative[i];
                    if (hasPort)
                        // dz/t = 1/2 at either thickness point: cancel it
                        // before multiplication to avoid tiny/large t^2.
                        coupling[i*ports+port] += scaledProduct({width,beamLength,sHalf,
                            weights[ix],strain[i],reduced.coupling[0],
                            static_cast<long double>(layer.polarity),0.5L});
                }
                for (size_t i=0;i<6;++i)
                    for (size_t j=i;j<6;++j)
                    {
                        // Neither an unweighted square nor the volume needs
                        // to fit the floating-point range on its own.
                        const long double m = scaledProduct({width,beamLength,sHalf,weights[ix],
                            zHalf,density,longitudinal[i],longitudinal[j]})+
                            scaledProduct({width,beamLength,sHalf,weights[ix],
                            zHalf,density,transverse[i],transverse[j]});
                        const long double k = scaledProduct({width,beamLength,sHalf,weights[ix],
                            zHalf,reduced.elasticity[0],strain[i],strain[j]});
                        mass[i*6+j] += m;
                        stiffness[i*6+j] += k;
                        if (i != j)
                        {
                            mass[j*6+i] += m;
                            stiffness[j*6+i] += k;
                        }
                    }
            }
        }
    }
    copyOperators(mass,result.mass);
    copyOperators(stiffness,result.stiffness);
    copyOperators(coupling,result.coupling);
    copyOperators(capacitance,result.capacitance);
    result.totalMass = checkedDouble(totalMass);
    result.validate();
    return result;
}
