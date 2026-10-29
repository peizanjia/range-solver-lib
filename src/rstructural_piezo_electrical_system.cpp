#include "rstructural_piezo_electrical_system.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
void require(bool condition,const char *message)
{
    if (!condition) throw std::invalid_argument(message);
}

void finite(const std::vector<double> &values)
{
    for (double value : values)
        require(std::isfinite(value),"Non-finite piezo electrical state");
}

// Keep products as mantissa/exponent pairs until the FINAL returned quantity.
// This must work on MSVC, where long double has the same range as double.
// Equal-exponent terms are summed with compensation before lower exponent groups
// are visited, preserving small signals after cancellation of enormous products.
// No absolute capacitance floor and no intermediate product range restriction.
class ScaledSum
{
    struct Term
    {
        long double mantissa;
        int exponent;
    };
    std::vector<Term> terms;

    static Term normalized(long double mantissa,int exponent)
    {
        if (mantissa == 0.0L) return {0.0L,0};
        int shift = 0;
        const long double fraction = std::frexp(mantissa,&shift);
        return {fraction,exponent+shift};
    }

public:
    void addProduct(double a,double b,double divisor=1.0)
    {
        if (a == 0.0 || b == 0.0) return;
        int ae = 0, be = 0, de = 0;
        const double am = std::frexp(a,&ae), bm = std::frexp(b,&be), dm = std::frexp(divisor,&de);
        terms.push_back(normalized(static_cast<long double>(am)*bm/dm,ae+be-de));
    }

    void add(double value) { addProduct(value,1.0); }

    // Form the displacement difference BEFORE multiplication, while retaining
    // its exponent separately (opposite large states may have |dx| > DBL_MAX).
    // This avoids subtracting independently rounded K*q and K*q0 products.
    void addDifferenceProduct(double coefficient,double value,double reference)
    {
        if (coefficient == 0.0 || value == reference) return;
        if (reference == 0.0) { addProduct(coefficient,value); return; }
        if (value == 0.0) { addProduct(-coefficient,reference); return; }
        int ve = 0, re = 0, ce = 0;
        const double vm = std::frexp(value,&ve), rm = std::frexp(reference,&re);
        const double cm = std::frexp(coefficient,&ce);
        const int common = std::max(ve,re);
        const Term difference = normalized(std::scalbn(static_cast<long double>(vm),ve-common)-
                                           std::scalbn(static_cast<long double>(rm),re-common),common);
        terms.push_back(normalized(difference.mantissa*cm,difference.exponent+ce));
    }

    double result(double divisor=1.0) const
    {
        int de = 0;
        const double dm = std::frexp(divisor,&de);
        auto sorted = terms;
        for (auto &term : sorted)
            term = normalized(term.mantissa/dm,term.exponent-de);
        std::sort(sorted.begin(),sorted.end(),[](const Term &a,const Term &b)
                  { return a.exponent > b.exponent; });
        Term total{0.0L,0};
        for (size_t i=0;i<sorted.size();)
        {
            const int exponent = sorted[i].exponent;
            long double sum = 0.0L, correction = 0.0L;
            do
            {
                const long double value = sorted[i++].mantissa;
                const long double next = sum+value;
                correction += std::abs(sum) >= std::abs(value) ? (sum-next)+value : (value-next)+sum;
                sum = next;
            } while (i<sorted.size() && sorted[i].exponent == exponent);
            const Term group = normalized(sum+correction,exponent);
            if (group.mantissa == 0.0L) continue;
            if (total.mantissa == 0.0L) total = group;
            else
            {
                const int common = std::max(total.exponent,group.exponent);
                total = normalized(std::scalbn(total.mantissa,total.exponent-common)+
                                   std::scalbn(group.mantissa,group.exponent-common),common);
            }
        }
        const long double extended = std::scalbn(total.mantissa,total.exponent);
        const double value = static_cast<double>(extended);
        require(std::isfinite(value) && (total.mantissa == 0.0L || value != 0.0),
                "Final piezo electrical result is not representable");
        return value;
    }
};
}

RStructuralPiezoElectricalSystem::RStructuralPiezoElectricalSystem(
    const RStructuralPiezoElementOperators &source,
    const std::vector<RStructuralPiezoElectricalBoundary> &settings,
    const std::vector<double> &referenceQ,const std::vector<double> &referenceVa)
    : operators(source),q0(referenceQ),actuatorVoltages0(referenceVa)
{
    operators.validate();
    const size_t n = operators.dofs, p = operators.portIds.size();
    require(settings.size() == p,"Each electrical port requires exactly one boundary");
    for (size_t i=0;i<p;++i)
        for (size_t j=0;j<p;++j)
            if (i != j)
                require(operators.capacitance[i*p+j] == 0.0,
                        "Shared inter-port electrical circuits are unsupported");

    // Match exact IDs, preserving physical operator ordering, never trim/rename.
    boundaries.resize(p);
    std::vector<bool> assigned(p,false);
    for (const auto &boundary : settings)
    {
        size_t index = p;
        for (size_t i=0;i<p;++i)
            if (boundary.portId == operators.portIds[i]) { index = i; break; }
        require(index != p,"Unknown electrical boundary port ID");
        require(!assigned[index],"Duplicate electrical boundary port ID");
        assigned[index] = true;
        require(std::isfinite(boundary.fixedCharge) &&
                std::isfinite(boundary.feedbackCapacitance) &&
                std::isfinite(boundary.referenceCharge),"Non-finite electrical boundary");
        switch (boundary.role)
        {
            case RStructuralPiezoElectricalRole::VoltageActuator:
            case RStructuralPiezoElectricalRole::ShortCircuitChargeSensor:
                require(boundary.fixedCharge == 0.0 && boundary.feedbackCapacitance == 0.0 &&
                        boundary.referenceCharge == 0.0,"Inapplicable electrical boundary fields");
                break;
            case RStructuralPiezoElectricalRole::FixedChargeOpenCircuit:
                require(boundary.feedbackCapacitance == 0.0 && boundary.referenceCharge == 0.0,
                        "Inapplicable open-circuit fields");
                break;
            case RStructuralPiezoElectricalRole::IdealChargeAmplifier:
                require(boundary.fixedCharge == 0.0 && boundary.feedbackCapacitance > 0.0,
                        "Ideal charge amplifier requires positive Cf and no fixed charge");
                break;
            default:
                throw std::invalid_argument("Unknown piezo electrical role");
        }
        boundaries[index] = boundary;
    }

    condensedOperators.dofs = operators.dofs;
    condensedOperators.mass = operators.mass;
    condensedOperators.stiffness.resize(n*n);
    condensedOperators.fixedChargeForce.resize(n);
    for (const auto &boundary : boundaries)
        if (boundary.role == RStructuralPiezoElectricalRole::VoltageActuator)
            condensedOperators.actuatorPortIds.push_back(boundary.portId);
    const size_t a = condensedOperators.actuatorPortIds.size();
    condensedOperators.coupling.resize(n*a);
    for (size_t i=0;i<n;++i)
    {
        ScaledSum bias;
        size_t actuator = 0;
        for (size_t s=0;s<p;++s)
        {
            if (boundaries[s].role == RStructuralPiezoElectricalRole::VoltageActuator)
                condensedOperators.coupling[i*a+actuator++] = operators.coupling[i*p+s];
            if (boundaries[s].role == RStructuralPiezoElectricalRole::FixedChargeOpenCircuit)
                bias.addProduct(operators.coupling[i*p+s],boundaries[s].fixedCharge,
                                operators.capacitance[s*p+s]);
        }
        condensedOperators.fixedChargeForce[i] = bias.result();
        for (size_t j=0;j<n;++j)
        {
            ScaledSum stiffness;
            stiffness.add(operators.stiffness[i*n+j]);
            for (size_t s=0;s<p;++s)
                if (boundaries[s].role == RStructuralPiezoElectricalRole::FixedChargeOpenCircuit)
                    stiffness.addProduct(operators.coupling[i*p+s],operators.coupling[j*p+s],
                                         operators.capacitance[s*p+s]);
            condensedOperators.stiffness[i*n+j] = stiffness.result();
        }
    }
    initialResponse = absoluteResponse(q0,actuatorVoltages0);
}

const RStructuralPiezoElectricalCondensedOperators &RStructuralPiezoElectricalSystem::condensed() const
{ return condensedOperators; }
const std::vector<double> &RStructuralPiezoElectricalSystem::referenceDisplacement() const
{ return q0; }
const std::vector<double> &RStructuralPiezoElectricalSystem::referenceActuatorVoltages() const
{ return actuatorVoltages0; }
const RStructuralPiezoElectricalResponse &RStructuralPiezoElectricalSystem::referenceResponse() const
{ return initialResponse; }

RStructuralPiezoElectricalResponse RStructuralPiezoElectricalSystem::absoluteResponse(
    const std::vector<double> &q,const std::vector<double> &voltages) const
{
    const size_t n = operators.dofs, p = boundaries.size();
    require(q.size() == n && voltages.size() == condensedOperators.actuatorPortIds.size(),
            "Piezo electrical state dimensions mismatch");
    finite(q);
    finite(voltages);
    RStructuralPiezoElectricalResponse response;
    response.ports.resize(p);
    response.internalForce.resize(n);
    response.deltaInternalForce.assign(n,0.0);
    size_t actuator = 0;
    for (size_t s=0;s<p;++s)
    {
        auto &port = response.ports[s];
        const auto &boundary = boundaries[s];
        port.portId = boundary.portId;
        port.role = boundary.role;
        const double capacitance = operators.capacitance[s*p+s];
        if (boundary.role == RStructuralPiezoElectricalRole::FixedChargeOpenCircuit)
        {
            ScaledSum voltage;
            voltage.add(boundary.fixedCharge);
            for (size_t i=0;i<n;++i) voltage.addProduct(operators.coupling[i*p+s],q[i]);
            port.voltage = voltage.result(capacitance);
            port.charge = boundary.fixedCharge;
        }
        else
        {
            if (boundary.role == RStructuralPiezoElectricalRole::VoltageActuator)
                port.voltage = voltages[actuator++];
            ScaledSum charge;
            charge.addProduct(capacitance,port.voltage);
            for (size_t i=0;i<n;++i) charge.addProduct(-operators.coupling[i*p+s],q[i]);
            port.charge = charge.result();
            if (boundary.role == RStructuralPiezoElectricalRole::IdealChargeAmplifier)
            {
                ScaledSum amplifier;
                amplifier.add(-port.charge);
                amplifier.add(boundary.referenceCharge);
                port.amplifierVoltage = amplifier.result(boundary.feedbackCapacitance);
            }
        }
    }
    for (size_t i=0;i<n;++i)
    {
        ScaledSum force;
        for (size_t j=0;j<n;++j) force.addProduct(operators.stiffness[i*n+j],q[j]);
        for (size_t s=0;s<p;++s) force.addProduct(operators.coupling[i*p+s],response.ports[s].voltage);
        response.internalForce[i] = force.result();
    }
    return response;
}

RStructuralPiezoElectricalResponse RStructuralPiezoElectricalSystem::evaluate(
    const std::vector<double> &q,const std::vector<double> &voltages) const
{
    auto response = absoluteResponse(q,voltages);
    // Compute increments from increment states rather than subtracting biased
    // absolute readouts: large Qref/fixed-charge offsets must not erase signals.
    std::vector<double> dva(voltages.size());
    for (size_t a=0;a<dva.size();++a)
    {
        ScaledSum difference;
        difference.add(voltages[a]); difference.add(-actuatorVoltages0[a]);
        dva[a] = difference.result();
    }
    const size_t n = operators.dofs, p = boundaries.size(), a = dva.size();
    for (size_t i=0;i<response.internalForce.size();++i)
    {
        ScaledSum force;
        for (size_t j=0;j<n;++j)
            force.addDifferenceProduct(condensedOperators.stiffness[i*n+j],q[j],q0[j]);
        for (size_t j=0;j<a;++j)
            force.addProduct(condensedOperators.coupling[i*a+j],dva[j]);
        response.deltaInternalForce[i] = force.result();
    }
    size_t actuator = 0;
    for (size_t s=0;s<response.ports.size();++s)
    {
        auto &port = response.ports[s];
        ScaledSum delta;
        for (size_t i=0;i<n;++i)
        {
            const double sign = port.role == RStructuralPiezoElectricalRole::FixedChargeOpenCircuit ? 1.0 : -1.0;
            delta.addDifferenceProduct(sign*operators.coupling[i*p+s],q[i],q0[i]);
        }
        if (port.role == RStructuralPiezoElectricalRole::FixedChargeOpenCircuit)
            port.deltaVoltage = delta.result(operators.capacitance[s*p+s]);
        else
        {
            if (port.role == RStructuralPiezoElectricalRole::VoltageActuator)
                port.deltaVoltage = dva[actuator++];
            delta.addProduct(operators.capacitance[s*p+s],port.deltaVoltage);
            port.deltaCharge = delta.result();
            if (port.role == RStructuralPiezoElectricalRole::IdealChargeAmplifier)
            {
                ScaledSum amplifier;
                amplifier.add(-port.deltaCharge);
                port.deltaAmplifierVoltage = amplifier.result(boundaries[s].feedbackCapacitance);
            }
        }
    }
    return response;
}

void RStructuralPiezoElectricalSystem::evaluateInto(
    const std::vector<double> &q,const std::vector<double> &voltages,
    RStructuralPiezoElectricalResponse &workspace) const
{
    auto response = evaluate(q,voltages);
    // Vector/string swaps do not allocate: no partial workspace update on errors.
    workspace.internalForce.swap(response.internalForce);
    workspace.deltaInternalForce.swap(response.deltaInternalForce);
    workspace.ports.swap(response.ports);
}
