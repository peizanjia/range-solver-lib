#include "rstructural_piezo_native_readout.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok,const char *message)
{ if (!ok) throw std::invalid_argument(message); }
double finite(double x)
{ require(std::isfinite(x),"Native piezo readout is nonfinite or unrepresentable"); return x; }

// Compensated, exponent-scaled arithmetic keeps small increments separate from
// large biases, including on MSVC where long double has binary64 precision.
struct Number { double hi=0.0,lo=0.0; int exponent=0; };
Number normalize(double hi,double lo,int exponent)
{
    const double s=hi+lo,b=s-hi,tail=(hi-(s-b))+(lo-b);
    if(s==0.0) return {};
    int shift=0; const double h=std::frexp(s,&shift);
    return {h,std::ldexp(tail,-shift),exponent+shift};
}
Number product(std::initializer_list<double> factors)
{
    Number n{1.0,0.0,0};
    for(double x:factors) {
        int e=0; const double f=std::frexp(finite(x),&e);
        if(f==0.0) return {};
        const double h=n.hi*f;
        n=normalize(h,std::fma(n.hi,f,-h)+n.lo*f,n.exponent+e);
    }
    return n;
}
Number sum(std::initializer_list<Number> terms)
{
    std::vector<Number> ordered(terms);
    std::sort(ordered.begin(),ordered.end(),[](const Number &a,const Number &b){return a.exponent>b.exponent;});
    Number n;
    for(const auto &v:ordered) {
        if(v.hi==0.0) continue;
        if(n.hi==0.0) { n=v; continue; }
        const int e=std::max(n.exponent,v.exponent);
        const double a=std::ldexp(n.hi,n.exponent-e),b=std::ldexp(v.hi,v.exponent-e),s=a+b,t=s-a;
        n=normalize(s,(a-(s-t))+(b-t)+std::ldexp(n.lo,n.exponent-e)+std::ldexp(v.lo,v.exponent-e),e);
    }
    return n;
}
double value(Number n,double divisor=1.0)
{
    require(finite(divisor)>0.0,"Native piezo readout divisor must be positive");
    int e=0; const double f=std::frexp(divisor,&e),q=n.hi/f;
    n=normalize(q,(std::fma(-q,f,n.hi)+n.lo)/f,n.exponent-e);
    return finite(std::ldexp(n.hi+n.lo,n.exponent));
}
double lengthIncrement(const RStructuralTrussDynamicsSystem &s,const RStructuralTrussMember &m,
                       const std::vector<double> &q,const std::vector<double> &q0)
{
    RStructuralTrussVector3 c{},d{},next{};
    for(std::size_t a=0;a<3;++a) {
        const auto i=3*m.node1+a,j=3*m.node2+a;
        c[a]=value(sum({product({s.referenceNodes[m.node2][a]}),product({-s.referenceNodes[m.node1][a]}),
                       product({q0[j]}),product({-q0[i]})}));
        d[a]=value(sum({product({q[j]}),product({-q0[j]}),product({-q[i]}),product({q0[i]})}));
        next[a]=value(sum({product({c[a]}),product({d[a]})}));
    }
    const double l0=finite(std::hypot(c[0],c[1],c[2])),l=finite(std::hypot(next[0],next[1],next[2]));
    require(l0>0.0 && l>0.0,"Native piezo readout chord collapsed");
    const double denominator=value(sum({product({.5,l0}),product({.5,l})}));
    // |c+d|-|c|=(c.d+.5*d.d)/(.5*|c+d|+.5*|c|).
    return value(sum({product({c[0],d[0]}),product({.5,d[0],d[0]}),
                      product({c[1],d[1]}),product({.5,d[1],d[1]}),
                      product({c[2],d[2]}),product({.5,d[2],d[2]})}),denominator);
}
double inputValue(const RStructuralTrussInput &input)
{
    if(const auto *v=std::get_if<RStructuralTrussVoltage>(&input)) return finite(v->value);
    if(const auto *q=std::get_if<RStructuralTrussCharge>(&input)) return finite(q->value);
    throw std::invalid_argument("Native piezo readout requires an electrical boundary");
}
bool sameBoundary(const RStructuralTrussInput &a,const RStructuralTrussInput &b)
{ return a.index()==b.index() && inputValue(a)==inputValue(b); }
bool relativeMatch(double residual,double scale)
{
    // No one-volt/one-coulomb absolute floor. Zero scale requires exact zero.
    return scale==0.0 ? residual==0.0 : std::fabs(residual)/scale<=128.0*std::numeric_limits<double>::epsilon();
}
RStructuralTrussDynamicsSystem approvedReadoutSystem(
    const RStructuralTrussDynamicsModelAdapter &adapter,
    const RStructuralTrussDynamicsSystem &approved)
{
    auto original = approved;
    original.physicalDamping = adapter.system().physicalDamping;
    require(RStructuralTrussDynamics(original).systemIdentity()==adapter.driver().systemIdentity(),
            "Approved piezo readout system changes frozen geometry, laws, mass or supports");
    return approved;
}
}

RStructuralPiezoNativeReadout::RStructuralPiezoNativeReadout(
    const RStructuralTrussDynamicsModelAdapter &adapter,const RStructuralPiezoReadoutSetup &setup)
    : RStructuralPiezoNativeReadout(adapter,setup,adapter.system())
{}

RStructuralPiezoNativeReadout::RStructuralPiezoNativeReadout(
    const RStructuralTrussDynamicsModelAdapter &adapter,const RStructuralPiezoReadoutSetup &setup,
    const RStructuralTrussDynamicsSystem &approvedSystem)
    : system(approvedReadoutSystem(adapter,approvedSystem)),driver(system),catalog(setup)
{
    catalog.validate();
    const auto &maps=adapter.maps();
    require(catalog.ports.empty() || catalog.ports.size()==maps.memberByPortId.size(),
            "Native piezo readout catalog must cover every electrical port");
    for(const auto &p:catalog.ports) {
        const auto found=maps.memberByPortId.find(p.portId);
        require(found!=maps.memberByPortId.end(),"Native piezo readout references an unknown port");
        const auto index=found->second;
        require(std::holds_alternative<RStructuralTrussPiezoLaw>(system.members.at(index).law),"Native readout port is not piezoelectric");
        const auto &input=adapter.initialCompleteInputs().memberInputs.at(index);
        if(p.readout==RStructuralPiezoReadoutMode::OpenCircuitVoltage)
            require(std::holds_alternative<RStructuralTrussCharge>(input),"Open sensor requires fixed charge");
        else {
            require(std::holds_alternative<RStructuralTrussVoltage>(input),"Readout requires voltage boundary");
            if(p.role==RStructuralPiezoPortRole::Sensor)
                require(inputValue(input)==0.0,"Short/amplifier sensor requires zero voltage");
        }
        readoutByPortId.emplace(p.portId,memberIndices.size());
        memberIndices.push_back(index);
    }
    const auto initial=adapter.sourceInitialState();
    for(std::size_t i=0;i<catalog.ports.size();++i) if(catalog.ports[i].role==RStructuralPiezoPortRole::Sensor) {
        const auto index=memberIndices[i];
        require(sameBoundary(initial.appliedInputs.memberInputs.at(index),adapter.initialCompleteInputs().memberInputs.at(index)),
                "Native piezo initial sensor boundary differs from frozen source");
    }
    frozenReference.time=adapter.initialTime();
    frozenReference.displacements=initial.displacements;
    frozenReference.inputs=initial.appliedInputs;
    frozenReference.explicitWorkpoint=catalog.workpoint.available;
    if(catalog.workpoint.available) {
        require(catalog.workpoint.nodes.size()==maps.nodeByPhysicalId.size() &&
                catalog.workpoint.ports.size()==maps.memberByPortId.size(),"Native piezo workpoint must cover every physical node and port");
        frozenReference.inputs=adapter.initialCompleteInputs();
        for(const auto &node:catalog.workpoint.nodes) {
            const auto found=maps.nodeByPhysicalId.find(node.physicalNodeId);
            require(found!=maps.nodeByPhysicalId.end(),"Native piezo workpoint references an unknown physical node");
            for(std::size_t axis=0;axis<3;++axis) {
                const auto d=3*found->second+axis;
                frozenReference.displacements[d]=node.q0[axis];
                require(node.f0[axis]==frozenReference.inputs.externalForces[d],"Native piezo workpoint force differs from frozen initial background load");
            }
        }
        for(const auto &port:catalog.workpoint.ports) {
            const auto found=maps.memberByPortId.find(port.portId);
            require(found!=maps.memberByPortId.end(),"Native piezo workpoint references an unknown port");
            const auto &input=frozenReference.inputs.memberInputs.at(found->second);
            require(inputValue(input)==(std::holds_alternative<RStructuralTrussVoltage>(input)?port.V0:port.Q0),
                    "Native piezo workpoint differs from frozen initial electrical boundary");
        }
        // initialize enforces the source analytic prescribed position/velocity.
        driver.initialize(frozenReference.time,frozenReference.displacements,initial.velocities,frozenReference.inputs);
    }
    validateBoundaries(frozenReference.inputs);
    frozenReference.assembly=RStructuralTrussAssembly::evaluate(system.referenceNodes,frozenReference.displacements,
                                                               system.members,frozenReference.inputs.memberInputs);
    double externalNorm=0.0,internalNorm=0.0;
    for(std::size_t d=0;d<system.prescribedMask.size();++d) if(!system.prescribedMask[d]) {
        const double f=frozenReference.inputs.externalForces[d],i=frozenReference.assembly.internalForces[d];
        frozenReference.freeForceResidualNorm=finite(std::hypot(frozenReference.freeForceResidualNorm,finite(i-f)));
        externalNorm=finite(std::hypot(externalNorm,f)); internalNorm=finite(std::hypot(internalNorm,i));
    }
    frozenReference.forceScale=std::max(externalNorm,internalNorm);
    frozenReference.forceTolerance=finite(adapter.options().absForceTol+finite(adapter.options().relForceTol*frozenReference.forceScale));
    if(catalog.workpoint.available) {
        require(frozenReference.freeForceResidualNorm<=frozenReference.forceTolerance,"Native piezo explicit workpoint is not in free-force equilibrium");
        for(const auto &port:catalog.workpoint.ports) {
            const auto index=maps.memberByPortId.at(port.portId);
            const auto &actual=*frozenReference.assembly.memberResponses[index].electrical;
            RStructuralPiezoNativeElectricalResidual residual;
            residual.portId=port.portId;
            residual.voltageResidual=finite(actual.voltage-port.V0);
            residual.chargeResidual=finite(actual.conductorCharge-port.Q0);
            residual.voltageScale=std::max(std::fabs(actual.voltage),std::fabs(port.V0));
            residual.chargeScale=std::max(std::fabs(actual.conductorCharge),std::fabs(port.Q0));
            require(relativeMatch(residual.voltageResidual,residual.voltageScale) &&
                    relativeMatch(residual.chargeResidual,residual.chargeScale),"Native piezo explicit workpoint has an electrical residual");
            frozenReference.electricalResiduals.push_back(residual);
        }
        frozenReference.equilibriumCertified=true;
    }
}

void RStructuralPiezoNativeReadout::validateBoundaries(const RStructuralTrussDynamicsInputs &inputs) const
{
    require(inputs.memberInputs.size()==system.members.size(),"Native piezo readout input count mismatch");
    for(std::size_t i=0;i<catalog.ports.size();++i) if(catalog.ports[i].role==RStructuralPiezoPortRole::Sensor) {
        const auto index=memberIndices[i];
        require(sameBoundary(inputs.memberInputs[index],frozenReference.inputs.memberInputs[index]),
                "Native piezo sensor-only electrical boundary was changed");
    }
}

std::vector<RStructuralPiezoNativeReadoutSample> RStructuralPiezoNativeReadout::observe(const RStructuralTrussDynamicsState &state) const
{
    validateBoundaries(state.appliedInputs);
    // Complete accepted-state validation, including identity and prescribed q/v/a.
    driver.refresh(state,state.appliedInputs);
    const auto assembly=RStructuralTrussAssembly::evaluate(system.referenceNodes,state.displacements,system.members,state.appliedInputs.memberInputs);
    return observe(state,assembly);
}

std::vector<RStructuralPiezoNativeReadoutSample> RStructuralPiezoNativeReadout::observe(
    const RStructuralTrussDynamicsState &state,const RStructuralTrussAssemblyResult &assembly) const
{
    validateBoundaries(state.appliedInputs);
    require(state.systemIdentity==driver.systemIdentity() && state.displacements.size()==system.referenceNodes.size()*3 &&
            assembly.memberResponses.size()==system.members.size(),"Native piezo readout state/assembly identity or dimensions mismatch");
    for(double q:state.displacements) finite(q);
    std::vector<RStructuralPiezoNativeReadoutSample> samples;
    samples.reserve(catalog.ports.size());
    for(std::size_t i=0;i<catalog.ports.size();++i) {
        const auto &p=catalog.ports[i]; const auto index=memberIndices[i];
        const auto &member=system.members[index]; const auto &law=std::get<RStructuralTrussPiezoLaw>(member.law);
        require(assembly.memberResponses[index].electrical.has_value(),"Native piezo readout assembly lacks electrical response");
        const auto &actual=*assembly.memberResponses[index].electrical;
        const auto &baseline=*frozenReference.assembly.memberResponses[index].electrical;
        const double dl=lengthIncrement(system,member,state.displacements,frozenReference.displacements);
        RStructuralPiezoNativeReadoutSample sample;
        sample.portId=p.portId; sample.mode=p.readout; sample.unit=rStructuralPiezoReadoutUnit(p.readout);
        sample.voltage=actual.voltage; sample.charge=actual.conductorCharge;
        const auto &input=state.appliedInputs.memberInputs[index]; const auto &refInput=frozenReference.inputs.memberInputs[index];
        require(input.index()==refInput.index(),"Native piezo readout cannot switch voltage/charge boundary kind");
        if(std::holds_alternative<RStructuralTrussVoltage>(input)) {
            sample.deltaVoltage=value(sum({product({inputValue(input)}),product({-inputValue(refInput)})}));
            sample.deltaCharge=value(sum({product({law.capacitance,sample.deltaVoltage}),product({-law.coupling,dl})}));
        } else {
            sample.deltaCharge=value(sum({product({inputValue(input)}),product({-inputValue(refInput)})}));
            sample.deltaVoltage=value(sum({product({sample.deltaCharge}),product({law.coupling,dl})}),law.capacitance);
        }
        switch(p.readout) {
        case RStructuralPiezoReadoutMode::None:
            sample.actionRecord=true; sample.unit=QStringLiteral("V");
            sample.value=sample.voltage; sample.deltaValue=sample.deltaVoltage; break;
        case RStructuralPiezoReadoutMode::OpenCircuitVoltage:
            sample.value=sample.voltage; sample.deltaValue=sample.deltaVoltage; break;
        case RStructuralPiezoReadoutMode::ShortCircuitCharge:
        case RStructuralPiezoReadoutMode::ActuatorCharge:
            sample.value=sample.charge; sample.deltaValue=sample.deltaCharge; break;
        case RStructuralPiezoReadoutMode::IdealChargeAmplifier:
            // Preserve a small transferred charge independently of large Qref.
            sample.value=value(sum({product({p.VoutBias,p.Cf}),product({double(p.amplifierPolarity),baseline.conductorCharge}),
                                   product({-double(p.amplifierPolarity),p.Qref}),product({double(p.amplifierPolarity),sample.deltaCharge})}),p.Cf);
            sample.deltaValue=value(product({double(p.amplifierPolarity),sample.deltaCharge}),p.Cf); break;
        }
        for(double v:{sample.voltage,sample.charge,sample.deltaVoltage,sample.deltaCharge,sample.value,sample.deltaValue}) finite(v);
        samples.push_back(sample);
    }
    return samples;
}
