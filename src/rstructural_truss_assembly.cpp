#include "rstructural_truss_assembly.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace
{
double checked(double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("Structural truss value is non-finite or unrepresentable");
    return value;
}
double positive(double value)
{
    if (checked(value) <= 0.0)
        throw std::invalid_argument("Structural truss parameter must be strictly positive");
    return value;
}

// A two-word significand and an unrestricted binary exponent keep temporary
// products/sums representable without relying on wider MSVC long double.
// fma retains product roundoff when two individually huge products cancel.
struct ScaledValue { double high = 0.0; double low = 0.0; int exponent = 0; };

ScaledValue normalize(double high, double low, int exponent)
{
    const double sum = high+low;
    const double lowPart = sum-high;
    const double tail = (high-(sum-lowPart))+(low-lowPart);
    if (sum == 0.0) return {};
    int shift = 0;
    const double fraction = std::frexp(sum,&shift);
    return {fraction,std::ldexp(tail,-shift),exponent+shift};
}

ScaledValue product(std::initializer_list<double> factors)
{
    ScaledValue result{1.0,0.0,0};
    for (double factor : factors)
    {
        int exponent = 0;
        const double fraction = std::frexp(checked(factor),&exponent);
        if (fraction == 0.0) return {};
        const double high = result.high*fraction;
        const double low = std::fma(result.high,fraction,-high)+result.low*fraction;
        result = normalize(high,low,result.exponent+exponent);
    }
    return result;
}

ScaledValue add(ScaledValue left, ScaledValue right)
{
    if (left.high == 0.0) return right;
    if (right.high == 0.0) return left;
    const int exponent = std::max(left.exponent,right.exponent);
    const double a = std::ldexp(left.high,left.exponent-exponent);
    const double b = std::ldexp(right.high,right.exponent-exponent);
    const double sum = a+b;
    const double bPart = sum-a;
    const double error = (a-(sum-bPart))+(b-bPart);
    const double low = error+std::ldexp(left.low,left.exponent-exponent)+
                             std::ldexp(right.low,right.exponent-exponent);
    return normalize(sum,low,exponent);
}

ScaledValue sum(std::initializer_list<ScaledValue> terms)
{
    std::vector<ScaledValue> ordered(terms);
    std::sort(ordered.begin(),ordered.end(),[](const ScaledValue &a, const ScaledValue &b) {
        return a.exponent > b.exponent;
    });
    ScaledValue result;
    for (const auto &term : ordered) result = add(result,term);
    return result;
}

ScaledValue divide(ScaledValue scaled, double divisor)
{
    int exponent = 0;
    const double fraction = std::frexp(positive(divisor),&exponent);
    const double quotient = scaled.high/fraction;
    const double correction = (std::fma(-quotient,fraction,scaled.high)+scaled.low)/fraction;
    return normalize(quotient,correction,scaled.exponent-exponent);
}

double value(ScaledValue scaled, double divisor = 1.0)
{
    const auto result = divide(scaled,divisor);
    return checked(std::ldexp(result.high+result.low,result.exponent));
}
template<class T> const T &require(const RStructuralTrussInput &input)
{
    const auto *value = std::get_if<T>(&input);
    if (!value) throw std::invalid_argument("Structural truss law/input type mismatch");
    return *value;
}
}

RStructuralTrussElementResult RStructuralTrussAssembly::evaluateElement(
    const RStructuralTrussVector3 &x1, const RStructuralTrussVector3 &x2,
    const RStructuralTrussLaw &law, const RStructuralTrussInput &input)
{
    RStructuralTrussVector3 n{};
    for (std::size_t i=0;i<3;++i) n[i] = checked(checked(x2[i])-checked(x1[i]));
    RStructuralTrussElementResult result;
    result.length = positive(std::hypot(n[0],n[1],n[2]));
    for (double &value : n) value /= result.length;

    if (const auto *passive = std::get_if<RStructuralTrussPassiveLaw>(&law))
    {
        require<RStructuralTrussNoInput>(input);
        const double delta = checked(result.length-positive(passive->referenceLength));
        result.axialTangent = positive(positive(passive->axialRigidity)/passive->referenceLength);
        result.axialForce = value(product({result.axialTangent,delta}));
        result.storedEnergy = value(product({0.5,result.axialTangent,delta,delta}));
        result.mechanicalPotential = result.storedEnergy;
    }
    else if (const auto *piezo = std::get_if<RStructuralTrussPiezoLaw>(&law))
    {
        const double delta = checked(result.length-positive(piezo->referenceLength));
        const double k = positive(piezo->axialStiffness);
        const double g = checked(piezo->coupling);
        const double ce = positive(piezo->capacitance);
        double voltage = 0.0, charge = 0.0;
        result.axialTangent = k;
        const bool voltageDriven = std::holds_alternative<RStructuralTrussVoltage>(input);
        if (voltageDriven)
        {
            voltage = checked(require<RStructuralTrussVoltage>(input).value);
            charge = value(sum({product({ce,voltage}),product({-g,delta})}));
        }
        else
        {
            charge = checked(require<RStructuralTrussCharge>(input).value);
            voltage = value(sum({product({charge}),product({g,delta})}),ce);
            result.axialTangent = positive(value(sum({product({k,ce}),product({g,g})}),ce));
        }
        result.axialForce = value(sum({product({k,delta}),product({g,voltage})}));
        const auto elasticEnergy = product({0.5,k,delta,delta});
        const auto electricEnergy = product({0.5,ce,voltage,voltage});
        const double enthalpy = value(sum({elasticEnergy,product({g,delta,voltage}),
                                           product({-0.5,ce,voltage,voltage})}));
        result.storedEnergy = value(sum({elasticEnergy,electricEnergy}));
        result.mechanicalPotential = voltageDriven ? enthalpy : result.storedEnergy;
        result.electrical = RStructuralTrussElectricalResult{voltage,charge,enthalpy};
    }
    else
    {
        const auto &cable = std::get<RStructuralTrussCableLaw>(law);
        const double l0 = positive(require<RStructuralTrussCableLength>(input).value);
        const double ea = positive(cable.axialRigidity);
        if (checked(cable.initialTension) < 0.0)
            throw std::invalid_argument("Structural truss cable prestress must be nonnegative");
        // Preserve the established tension-only law while avoiding temporary
        // N*N overflow (MSVC long double has double precision/range).
        const double k = positive(ea/l0);
        const double delta = checked(result.length-l0);
        const double trialTension = value(sum({product({cable.initialTension}),product({k,delta})}));
        const bool taut = trialTension > 0.0;
        result.axialForce = taut ? trialTension : 0.0;
        result.axialTangent = taut ? k : 0.0;
        result.storedEnergy = value(product({0.5,result.axialForce,result.axialForce}),k);
        result.mechanicalPotential = result.storedEnergy;
        result.cableTaut = taut;
        result.storedEnergyLengthDerivative = value(sum({
            divide(product({-result.axialForce,result.length}),l0),
            divide(product({0.5,result.axialForce,result.axialForce}),ea)}));
    }

    for (std::size_t i=0;i<3;++i)
    {
        result.internalForce[i] = checked(-result.axialForce*n[i]);
        result.internalForce[i+3] = -result.internalForce[i];
        for (std::size_t j=0;j<3;++j)
        {
            const double nn = n[i]*n[j];
            const double a = value(sum({product({result.axialTangent,nn,result.length}),
                                       product({result.axialForce,(i == j ? 1.0 : 0.0)-nn})}),result.length);
            result.tangent[i][j] = result.tangent[i+3][j+3] = a;
            result.tangent[i][j+3] = result.tangent[i+3][j] = -a;
        }
    }
    return result;
}

RStructuralTrussAssemblyResult RStructuralTrussAssembly::evaluate(
    const std::vector<RStructuralTrussVector3> &referenceNodes,
    const std::vector<double> &displacements,
    const std::vector<RStructuralTrussMember> &members,
    const std::vector<RStructuralTrussInput> &inputs)
{
    if (referenceNodes.size() > std::numeric_limits<uint>::max()/3 ||
        referenceNodes.size() > std::numeric_limits<std::size_t>::max()/3 ||
        displacements.size() != referenceNodes.size()*3 || inputs.size() != members.size())
        throw std::invalid_argument("Structural truss assembly dimensions are invalid");
    std::vector<RStructuralTrussVector3> positions(referenceNodes.size());
    for (std::size_t node=0;node<referenceNodes.size();++node)
        for (std::size_t i=0;i<3;++i)
            positions[node][i] = checked(checked(referenceNodes[node][i])+checked(displacements[3*node+i]));
    std::unordered_set<std::uint64_t> ids;
    for (const auto &member : members)
        if (member.physicalElementId == 0 || !ids.insert(member.physicalElementId).second ||
            member.node1 >= positions.size() || member.node2 >= positions.size() || member.node1 == member.node2)
            throw std::invalid_argument("Structural truss member identity or connectivity is invalid");

    RStructuralTrussAssemblyResult result;
    result.internalForces.resize(displacements.size(),0.0);
    result.tangent.setNRows(uint(displacements.size()));
    result.memberResponses.reserve(members.size());
    for (std::size_t e=0;e<members.size();++e)
    {
        const auto &member = members[e];
        const auto element = evaluateElement(positions[member.node1],positions[member.node2],member.law,inputs[e]);
        result.storedEnergy = checked(result.storedEnergy+element.storedEnergy);
        result.mechanicalPotential = checked(result.mechanicalPotential+element.mechanicalPotential);
        for (std::size_t i=0;i<6;++i)
        {
            const uint row = uint(3*(i<3 ? member.node1 : member.node2)+i%3);
            result.internalForces[row] = checked(result.internalForces[row]+element.internalForce[i]);
            for (std::size_t j=0;j<6;++j)
            {
                const uint column = uint(3*(j<3 ? member.node1 : member.node2)+j%3);
                const double value = element.tangent[i][j];
                if (value != 0.0)
                {
                    checked(result.tangent.findValue(row,column)+value);
                    result.tangent.addValue(row,column,value);
                }
            }
        }
        result.memberResponses.push_back(element);
    }
    return result;
}
