#include "rstructural_piezo_layered_dynamics.h"
#include "rparallel_execution_context.h"
#include <QCryptographicHash>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace {
using V=std::vector<double>;
void require(bool ok,const char *message) {if(!ok) throw std::invalid_argument(message);}
double checked(double x) {require(std::isfinite(x),"Nonfinite layered dynamics result"); return x;}
void finite(const V &v) {for(double x:v) checked(x);}
RParallelExecutionConfig serialPolicy() {
    RParallelExecutionConfig policy; policy.setMode(RParallelExecutionMode::SerialDeterministic); return policy;
}
V mv(const V &a,const V &x,size_t rows) {
    require(a.size()==rows*x.size(),"Layered dynamics matrix dimensions mismatch"); V y(rows);
    for(size_t i=0;i<rows;++i) for(size_t j=0;j<x.size();++j) y[i]=checked(y[i]+checked(a[i*x.size()+j]*x[j]));
    return y;
}
V transpose(const V &t,const V &x,size_t columns) {
    require(t.size()==x.size()*columns,"Layered dynamics projection dimensions mismatch"); V y(columns);
    for(size_t i=0;i<x.size();++i) for(size_t j=0;j<columns;++j) y[j]=checked(y[j]+checked(t[i*columns+j]*x[i]));
    return y;
}
double dot(const V &a,const V &b) {
    require(a.size()==b.size(),"Layered dynamics dot dimensions mismatch"); double x=0;
    for(size_t i=0;i<a.size();++i) x=checked(x+checked(a[i]*b[i])); return x;
}
double energy(const V &a,const V &x) {return checked(.5*dot(x,mv(a,x,x.size())));}
RSparseMatrix sparse(const V &a,uint n) {
    require(a.size()==size_t(n)*n,"Layered dynamics sparse dimensions mismatch"); finite(a);
    RSparseMatrix result; result.setNRows(n);
    for(uint i=0;i<n;++i) for(uint j=0;j<n;++j) if(a[size_t(i)*n+j]!=0) result.addValue(i,j,a[size_t(i)*n+j]);
    return result;
}
bool close(double a,double b) {
    checked(a); checked(b); if(a==b) return true;
    const double scale=std::max(std::fabs(a),std::fabs(b));
    return scale>0 && std::fabs(a-b)/scale<=2e-9;
}
bool residualClose(double a,double b,std::initializer_list<double> terms) {
    double scale=0; for(double x:terms) scale=checked(scale+std::fabs(checked(x)));
    return a==b || (scale>0 && std::fabs(checked(a-b))/scale<=2e-9);
}
void equal(const V &a,const V &b) {
    require(a.size()==b.size(),"Layered dynamics derived vector dimensions mismatch");
    for(size_t i=0;i<a.size();++i) require(close(a[i],b[i]),"Layered dynamics derived state mismatch");
}
struct Writer {
    QByteArray bytes;
    void integer(quint64 v) {for(uint i=0;i<8;++i) bytes.append(char((v>>(8*i))&255));}
    void real(double x) {checked(x); quint64 bits; static_assert(sizeof(x)==sizeof(bits),"binary64 required"); std::memcpy(&bits,&x,8); integer(bits);}
    void blob(const QByteArray &b) {integer(quint64(b.size())); bytes.append(b);}
    void text(const QString &s) {blob(s.toUtf8());}
    void vector(const V &v) {integer(v.size()); for(double x:v) real(x);}
    QByteArray hash() const {return QCryptographicHash::hash(bytes,QCryptographicHash::Sha256);}
};
void write(Writer &w,const RStructuralDynamicsCheckpoint &c) {
    const auto &s=c.state; w.real(s.time); w.vector(s.q); w.vector(s.v); w.vector(s.a); w.vector(s.externalForce);
    for(double x:{s.kineticEnergy,s.strainEnergy,s.externalWork,s.dampingDissipation,s.energyBalanceResidual,c.fixedStep,c.initialMechanicalEnergy}) w.real(x);
}
void write(Writer &w,const RStructuralPiezoElectricalResponse &e) {
    w.vector(e.internalForce); w.vector(e.deltaInternalForce); w.integer(e.ports.size());
    for(const auto &p:e.ports) {w.text(p.portId); w.integer(uint(p.role));
        for(double x:{p.voltage,p.charge,p.deltaVoltage,p.deltaCharge,p.amplifierVoltage,p.deltaAmplifierVoltage}) w.real(x);}
}
double voltageWork(const RStructuralPiezoElectricalResponse &left,const RStructuralPiezoElectricalResponse &right,bool jump) {
    require(left.ports.size()==right.ports.size(),"Layered dynamics electrical work dimensions mismatch"); double work=0;
    for(size_t i=0;i<left.ports.size();++i) {
        const auto &l=left.ports[i],&r=right.ports[i];
        require(l.portId==r.portId && l.role==r.role,"Layered dynamics electrical boundary changed");
        if(l.role==RStructuralPiezoElectricalRole::VoltageActuator) {
            const double voltage=jump?checked(.5*l.voltage+.5*r.voltage):r.voltage;
            // Both deltas share the compiler's immutable reference; retain
            // small transferred charge independently of a large absolute Q0.
            work=checked(work+checked(voltage*checked(r.deltaCharge-l.deltaCharge)));
        }
    }
    return work;
}
}

RStructuralPiezoLayeredDynamics::RStructuralPiezoLayeredDynamics(
    const RStructuralPiezoLayeredModelAdapter &adapter,const RMatrixSolverConf &conf,double time,
    BackgroundForce forcingCallback,const QString &forcingName)
    :mesh(adapter.mesh()),data(adapter.data()),electricalSystem(adapter.electrical()),solverConf(conf),
     mass(sparse(data.mass,data.freeDofs)),damping(sparse(data.damping,data.freeDofs)),stiffness(sparse(data.stiffness,data.freeDofs)),
     initialTime(checked(time)),force(std::move(forcingCallback)),forcingIdentity(forcingName)
{
    require(bool(force) && !forcingIdentity.isEmpty() && forcingIdentity.trimmed()==forcingIdentity,"Layered dynamics needs a stateless force callback and declared identity");
    require(data.freeDofs>0,"Layered dynamics needs free supported DOFs");
    require(conf.getType()==RMatrixSolverConf::CG || conf.getType()==RMatrixSolverConf::GMRES,"Layered dynamics needs a supported linear solver");
    V relative=data.q0; for(size_t i=0;i<relative.size();++i) relative[i]=checked(relative[i]-data.p[i]);
    z0=transpose(data.T,relative,data.freeDofs);
    require(data.initialZ.size()==z0.size(),"Layered dynamics initial absolute-coordinate dimensions mismatch");
    initialEta=data.initialZ; for(size_t i=0;i<initialEta.size();++i) initialEta[i]=checked(initialEta[i]-z0[i]);
    initialReducedEnergy=checked(energy(data.mass,data.initialVelocity)+energy(data.stiffness,initialEta));
    Writer w; w.text("range-fea/layered-native-Newmark/v1"); w.text(data.sourceFingerprint); w.text(forcingIdentity); w.real(initialTime);
    w.text("OpenMP/SerialDeterministic/one-worker/v1");
    w.text("linear-solve/Jacobi/symmetric-positive-diagonal-bounded-block-Jacobi-corrections-width-at-most-ten/compensated-original-matrix-componentwise-1e-12/twenty-residual-corrections/v4");
    w.integer(conf.getType()); w.integer(conf.getNInnerIterations()); w.integer(conf.getNOuterIterations()); w.real(conf.getSolverCvgValue());
    for(const auto *v:{&data.T,&data.p,&data.q0,&data.mass,&data.stiffness,&data.damping,&data.fullDamping,&data.input,
                       &data.backgroundOffset,&data.backgroundForce,&data.actuatorVoltages0,&initialEta,&data.initialVelocity}) w.vector(*v);
    identity=w.hash();
}

V RStructuralPiezoLayeredDynamics::forcing(double time) const {
    auto f=force(checked(time)); require(f.size()==6*mesh.nodeIds.size(),"Layered force callback returned invalid physical dimensions"); finite(f);
    mesh.projectPhysical(f,RStructuralPiezoPhysicalDomain::ForceAndMoment); // Reject unsupported force/moment components.
    return f;
}
RRVector RStructuralPiezoLayeredDynamics::rightHandSide(const V &physicalForce,const V &heldVa) const {
    require(heldVa.size()==actuatorPortIds().size(),"Layered voltage command must cover every actuator in compiled order"); finite(heldVa);
    auto f=transpose(data.T,mesh.projectPhysical(physicalForce,RStructuralPiezoPhysicalDomain::ForceAndMoment),data.freeDofs);
    RRVector rhs(data.backgroundOffset);
    for(size_t i=0;i<rhs.size();++i) {
        rhs[i]=checked(rhs[i]+f[i]);
        for(size_t j=0;j<heldVa.size();++j) rhs[i]=checked(rhs[i]+checked(data.input[i*heldVa.size()+j]*checked(heldVa[j]-data.actuatorVoltages0[j])));
    }
    return rhs;
}

RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::derive(
    const RStructuralDynamicsCheckpoint &checkpoint,const V &heldVa,const V &physicalForce) const
{
    const auto &s=checkpoint.state; const size_t n=data.fullDofs,r=data.freeDofs;
    require(s.q.size()==r && s.v.size()==r && s.a.size()==r,"Layered integrator checkpoint dimensions mismatch");
    finite(s.q); finite(s.v); finite(s.a); checked(s.time);
    V q=mv(data.T,s.q,n),v=mv(data.T,s.v,n),a=mv(data.T,s.a,n);
    for(size_t i=0;i<n;++i) q[i]=checked(q[i]+data.q0[i]);
    RStructuralPiezoLayeredDynamicsState state;
    state.sourceFingerprint=data.sourceFingerprint; state.forcingIdentity=forcingIdentity; state.driverIdentity=identity;
    state.integrator=checkpoint; state.heldActuatorVoltages=heldVa; state.perturbationForces=physicalForce;
    state.physicalDisplacements=mesh.reconstructPhysical(q); state.physicalVelocities=mesh.reconstructPhysical(v); state.physicalAccelerations=mesh.reconstructPhysical(a);
    state.electrical=electricalSystem.evaluate(q,heldVa);
    const auto ma=mv(mesh.operators.mass,a,n),cv=mv(data.fullDamping,v,n);
    const auto perturb=mesh.projectPhysical(physicalForce,RStructuralPiezoPhysicalDomain::ForceAndMoment);
    V reaction(n),scale(n);
    for(size_t i=0;i<n;++i) {
        reaction[i]=checked(checked(ma[i]+cv[i])+state.electrical.internalForce[i]-data.backgroundForce[i]-perturb[i]);
        scale[i]=checked(std::fabs(data.backgroundForce[i])+std::fabs(perturb[i]));
        for(size_t j=0;j<n;++j) scale[i]=checked(scale[i]+std::fabs(mesh.operators.mass[i*n+j]*a[j])+std::fabs(data.fullDamping[i*n+j]*v[j])+std::fabs(mesh.operators.stiffness[i*n+j]*q[j]));
        for(size_t j=0;j<state.electrical.ports.size();++j) scale[i]=checked(scale[i]+std::fabs(mesh.operators.coupling[i*state.electrical.ports.size()+j]*state.electrical.ports[j].voltage));
    }
    const auto residual=transpose(data.T,reaction,r);
    for(size_t j=0;j<r;++j) {
        double bound=0; for(size_t i=0;i<n;++i) bound=checked(bound+std::fabs(data.T[i*r+j])*scale[i]);
        if(!(bound==0?residual[j]==0:std::fabs(residual[j])/bound<=2e-8)) {
            const auto reducedMa=mv(data.mass,s.a,r),reducedCv=mv(data.damping,s.v,r),reducedKq=mv(data.stiffness,s.q,r);
            const auto rhs=rightHandSide(physicalForce,heldVa);double reducedNorm=0,rhsNorm=0;
            for(size_t i=0;i<r;++i){reducedNorm=std::hypot(reducedNorm,reducedMa[i]+reducedCv[i]+reducedKq[i]-rhs[i]);rhsNorm=std::hypot(rhsNorm,rhs[i]);}
            throw std::invalid_argument(QString("Layered native free-force residual certificate failed: row=%1 time=%2 residual=%3 bound=%4 reducedResidual=%5 reducedNorm=%6 rhsNorm=%7 solverTolerance=%8")
                .arg(qulonglong(j)).arg(s.time,0,'g',17).arg(residual[j],0,'g',17).arg(bound,0,'g',17)
                .arg(reducedMa[j]+reducedCv[j]+reducedKq[j]-rhs[j],0,'g',17).arg(reducedNorm,0,'g',17).arg(rhsNorm,0,'g',17)
                .arg(solverConf.getSolverCvgValue(),0,'g',17).toStdString());
        }
        state.freeForceResidual=checked(std::hypot(state.freeForceResidual,residual[j]));
    }
    state.physicalReactions=mesh.reconstructPhysical(reaction);
    state.kineticEnergy=energy(mesh.operators.mass,v);
    V voltages; for(const auto &p:state.electrical.ports) voltages.push_back(p.voltage);
    state.storedEnergy=checked(energy(mesh.operators.stiffness,q)+energy(mesh.operators.capacitance,voltages));
    std::map<QString,const RStructuralPiezoReadoutPort*> metadata;
    for(const auto &port:data.portMetadata) metadata.emplace(port.readout.portId,&port.readout);
    for(const auto &p:state.electrical.ports) {
        require(metadata.count(p.portId)!=0,"Layered native readout metadata missing"); const auto &spec=*metadata.at(p.portId);
        RStructuralPiezoLayeredDynamicsReadout readout; readout.portId=p.portId; readout.mode=spec.readout; readout.unit=rStructuralPiezoReadoutUnit(spec.readout);
        switch(spec.readout) {
        case RStructuralPiezoReadoutMode::None:
            readout.actionRecord=true; readout.unit="V"; readout.value=p.voltage; readout.deltaValue=p.deltaVoltage; break;
        case RStructuralPiezoReadoutMode::OpenCircuitVoltage:
            readout.value=p.voltage; readout.deltaValue=p.deltaVoltage; break;
        case RStructuralPiezoReadoutMode::ShortCircuitCharge:
        case RStructuralPiezoReadoutMode::ActuatorCharge:
            readout.value=p.charge; readout.deltaValue=p.deltaCharge; break;
        case RStructuralPiezoReadoutMode::IdealChargeAmplifier:
            readout.value=checked(spec.VoutBias-double(spec.amplifierPolarity)*p.amplifierVoltage);
            readout.deltaValue=checked(-double(spec.amplifierPolarity)*p.deltaAmplifierVoltage); break;
        }
        checked(readout.value); checked(readout.deltaValue); state.readouts.push_back(readout);
    }
    for(const auto *physical:{&state.physicalDisplacements,&state.physicalVelocities,&state.physicalAccelerations,&state.physicalReactions}) finite(*physical);
    return state;
}

RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::initialize() const {return initialize(data.actuatorVoltages0);}
RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::initialize(const V &heldVa) const {
    const RParallelExecutionContext execution(serialPolicy());
    const auto physicalForce=forcing(initialTime);
    const auto rhs=rightHandSide(physicalForce,heldVa);
    RStructuralDynamicsIntegrator integrator(mass,damping,stiffness,solverConf,{},
        [rhs](double,const RStructuralDynamicsState &,RRVector &f){f=rhs;});
    integrator.initialize(RRVector(initialEta),RRVector(data.initialVelocity),initialTime);
    auto state=derive(integrator.checkpoint(),heldVa,physicalForce);
    state.initialActuatorVoltages=heldVa;
    state.initialEnergy=checked(state.kineticEnergy+state.storedEnergy);
    state.integrityDigest=seal(state); validate(state); return state;
}

RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::trial(
    const RStructuralPiezoLayeredDynamicsState &accepted,const V &heldVa,double dt) const
{
    const double leftTime=accepted.integrator.state.time,rightTime=checked(leftTime+dt);
    return trialInterval(accepted,heldVa,dt,rightTime,0,0);
}

RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::trialOnGrid(
    const RStructuralPiezoLayeredDynamicsState &accepted,const V &heldVa,
    std::uint64_t stepIndex,double nominalDt,BackgroundForce intervalForce) const
{
    require(stepIndex>0 && stepIndex<=(std::uint64_t(1)<<53) && checked(nominalDt)>0,
            "Invalid layered canonical grid");
    const double leftTime=stepIndex==1?initialTime:checked(initialTime+double(stepIndex-1)*nominalDt);
    require(std::memcmp(&leftTime,&accepted.integrator.state.time,sizeof(double))==0,
            "Layered accepted time differs from preceding canonical grid point");
    const double rightTime=checked(initialTime+double(stepIndex)*nominalDt);
    return trialInterval(accepted,heldVa,rightTime-leftTime,rightTime,stepIndex,nominalDt,intervalForce);
}

RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::trialInterval(
    const RStructuralPiezoLayeredDynamicsState &accepted,const V &heldVa,
    double dt,double rightTime,std::uint64_t gridIndex,double nominalDt,const BackgroundForce &intervalForce) const
{
    // This dense native route promises bit-identical replay. Keep reductions
    // deterministic throughout scratch solves and restore the host's policy.
    const RParallelExecutionContext execution(serialPolicy());
    validate(accepted); require(checked(dt)>0,"Layered dynamics interval must be positive");
    const double leftTime=accepted.integrator.state.time;
    require(rightTime>leftTime,"Layered dynamics interval does not advance representable time");
    // Evaluate both endpoints independently. Changing voltage refreshes only
    // the right-side acceleration; caller's left sampling state is preserved.
    auto interval=[&](double time){
        if(!intervalForce)return forcing(time);
        auto values=intervalForce(checked(time));require(values.size()==6*mesh.nodeIds.size(),"Layered interval force dimensions mismatch");finite(values);
        mesh.projectPhysical(values,RStructuralPiezoPhysicalDomain::ForceAndMoment);return values;
    };
    const auto leftForce=interval(leftTime),rightForce=interval(rightTime);
    const auto leftRhs=rightHandSide(leftForce,heldVa),rightRhs=rightHandSide(rightForce,heldVa);
    RStructuralDynamicsIntegrator integrator(mass,damping,stiffness,solverConf,{},
        [leftTime,rightTime,leftRhs,rightRhs](double t,const RStructuralDynamicsState &,RRVector &f){
            require(t==leftTime || t==rightTime,"Unexpected Newmark force evaluation time"); f=t==leftTime?leftRhs:rightRhs;
        });
    integrator.restore(accepted.integrator); integrator.refreshAcceptedForce();
    const auto rightLeft=derive(integrator.checkpoint(),heldVa,leftForce);
    if(gridIndex) integrator.trialStepOnGrid(initialTime,gridIndex,nominalDt);
    else integrator.trialStep(dt);
    integrator.acceptStep();
    auto candidate=derive(integrator.checkpoint(),heldVa,rightForce);
    candidate.initialActuatorVoltages=accepted.initialActuatorVoltages; candidate.initialEnergy=accepted.initialEnergy;
    const auto leftF=mesh.projectPhysical(leftForce,RStructuralPiezoPhysicalDomain::ForceAndMoment);
    const auto rightF=mesh.projectPhysical(rightForce,RStructuralPiezoPhysicalDomain::ForceAndMoment);
    V deltaEta(data.freeDofs); for(size_t i=0;i<deltaEta.size();++i) deltaEta[i]=checked(candidate.integrator.state.q[i]-accepted.integrator.state.q[i]);
    const auto dq=mv(data.T,deltaEta,data.fullDofs); V meanForce(data.fullDofs);
    for(size_t i=0;i<dq.size();++i) meanForce[i]=checked(data.backgroundForce[i]+.5*leftF[i]+.5*rightF[i]);
    candidate.backgroundExternalWork=checked(accepted.backgroundExternalWork+dot(meanForce,dq));
    candidate.electricalWork=checked(accepted.electricalWork+voltageWork(accepted.electrical,rightLeft.electrical,true)+voltageWork(rightLeft.electrical,candidate.electrical,false));
    V meanReducedVelocity(data.freeDofs);
    for(size_t i=0;i<meanReducedVelocity.size();++i) meanReducedVelocity[i]=checked(.5*accepted.integrator.state.v[i]+.5*candidate.integrator.state.v[i]);
    const auto meanV=mv(data.T,meanReducedVelocity,data.fullDofs);
    const double dissipation=checked(dt*dot(meanV,mv(data.fullDamping,meanV,meanV.size())));
    require(dissipation>=0,"Layered native discrete viscous work is negative");
    candidate.viscousDissipation=checked(accepted.viscousDissipation+dissipation);
    candidate.energyBalanceResidual=checked(candidate.kineticEnergy+candidate.storedEnergy-candidate.initialEnergy-candidate.backgroundExternalWork-candidate.electricalWork+candidate.viscousDissipation);
    candidate.integrityDigest=seal(candidate); validate(candidate); return candidate;
}

QByteArray RStructuralPiezoLayeredDynamics::seal(const RStructuralPiezoLayeredDynamicsState &s) const {
    Writer w; w.text("range-fea/layered-native-checkpoint/v1"); w.text(s.sourceFingerprint); w.text(s.forcingIdentity); w.blob(s.driverIdentity); write(w,s.integrator);
    for(const auto *v:{&s.heldActuatorVoltages,&s.initialActuatorVoltages,&s.perturbationForces,&s.physicalDisplacements,&s.physicalVelocities,&s.physicalAccelerations,&s.physicalReactions}) w.vector(*v);
    write(w,s.electrical); w.integer(s.readouts.size());
    for(const auto &p:s.readouts) {w.text(p.portId); w.text(p.unit); w.integer(uint(p.mode)); w.integer(p.actionRecord?1:0); w.real(p.value); w.real(p.deltaValue);}
    for(double x:{s.kineticEnergy,s.storedEnergy,s.initialEnergy,s.backgroundExternalWork,s.electricalWork,s.viscousDissipation,s.energyBalanceResidual,s.freeForceResidual}) w.real(x);
    return w.hash();
}

void RStructuralPiezoLayeredDynamics::validate(const RStructuralPiezoLayeredDynamicsState &s) const {
    require(s.sourceFingerprint==data.sourceFingerprint && s.forcingIdentity==forcingIdentity && s.driverIdentity==identity,"Layered native checkpoint source/forcing/driver identity mismatch");
    require(s.integrityDigest.size()==32 && s.integrityDigest==seal(s),"Layered native checkpoint integrity seal mismatch");
    const auto &c=s.integrator; const auto &reduced=c.state;
    require(checked(c.fixedStep)>=0 && checked(reduced.time)>=initialTime,"Layered native checkpoint clock is invalid");
    require(s.heldActuatorVoltages.size()==actuatorPortIds().size() && s.initialActuatorVoltages.size()==actuatorPortIds().size(),"Layered native held actuator boundary count mismatch");
    finite(s.heldActuatorVoltages); finite(s.initialActuatorVoltages);
    const auto physicalForce=forcing(reduced.time); equal(s.perturbationForces,physicalForce);
    equal(reduced.externalForce,rightHandSide(physicalForce,s.heldActuatorVoltages));
    const auto actual=derive(c,s.heldActuatorVoltages,physicalForce);
    equal(s.physicalDisplacements,actual.physicalDisplacements); equal(s.physicalVelocities,actual.physicalVelocities);
    equal(s.physicalAccelerations,actual.physicalAccelerations); equal(s.physicalReactions,actual.physicalReactions);
    Writer suppliedElectrical,actualElectrical; write(suppliedElectrical,s.electrical); write(actualElectrical,actual.electrical);
    require(suppliedElectrical.bytes==actualElectrical.bytes,"Layered native checkpoint electrical response mismatch");
    require(s.readouts.size()==actual.readouts.size(),"Layered native checkpoint readout dimensions mismatch");
    for(size_t i=0;i<s.readouts.size();++i) {
        const auto &x=s.readouts[i],&a=actual.readouts[i];
        require(x.portId==a.portId && x.unit==a.unit && x.mode==a.mode && x.actionRecord==a.actionRecord && close(x.value,a.value) && close(x.deltaValue,a.deltaValue),"Layered native checkpoint readout mismatch");
    }
    const double kinetic=energy(data.mass,reduced.v),strain=energy(data.stiffness,reduced.q);
    require(close(reduced.kineticEnergy,kinetic) && close(reduced.strainEnergy,strain) && close(c.initialMechanicalEnergy,initialReducedEnergy),"Layered native incremental integrator energy mismatch");
    require(residualClose(reduced.energyBalanceResidual,checked(kinetic+strain-c.initialMechanicalEnergy-reduced.externalWork+reduced.dampingDissipation),
        {kinetic,strain,c.initialMechanicalEnergy,reduced.externalWork,reduced.dampingDissipation}),"Layered native integrator trapezoid damping diagnostic mismatch");
    require(reduced.dampingDissipation>=0 && s.viscousDissipation>=0,"Layered native checkpoint has negative damping work");
    require(close(s.kineticEnergy,actual.kineticEnergy) && close(s.storedEnergy,actual.storedEnergy) && close(s.freeForceResidual,actual.freeForceResidual),"Layered native checkpoint physical energy/residual mismatch");
    V initialQ=mv(data.T,initialEta,data.fullDofs),initialV=mv(data.T,data.initialVelocity,data.fullDofs);
    for(size_t i=0;i<initialQ.size();++i) initialQ[i]=checked(initialQ[i]+data.q0[i]);
    const auto initialElectric=electricalSystem.evaluate(initialQ,s.initialActuatorVoltages); V initialVoltages;
    for(const auto &p:initialElectric.ports) initialVoltages.push_back(p.voltage);
    const double initialFullEnergy=checked(energy(mesh.operators.mass,initialV)+energy(mesh.operators.stiffness,initialQ)+energy(mesh.operators.capacitance,initialVoltages));
    require(close(s.initialEnergy,initialFullEnergy),"Layered native checkpoint initial absolute energy mismatch");
    require(residualClose(s.energyBalanceResidual,checked(s.kineticEnergy+s.storedEnergy-s.initialEnergy-s.backgroundExternalWork-s.electricalWork+s.viscousDissipation),
        {s.kineticEnergy,s.storedEnergy,s.initialEnergy,s.backgroundExternalWork,s.electricalWork,s.viscousDissipation}),"Layered native checkpoint full energy ledger mismatch");
    if(c.fixedStep==0) {
        require(reduced.time==initialTime && s.backgroundExternalWork==0 && s.electricalWork==0 && s.viscousDissipation==0 && reduced.externalWork==0 && reduced.dampingDissipation==0,"Layered native initial checkpoint contains consumed time/work");
        equal(reduced.q,initialEta); equal(reduced.v,data.initialVelocity); equal(s.heldActuatorVoltages,s.initialActuatorVoltages);
    } else require(reduced.time>initialTime,"Layered native noninitial checkpoint did not advance time");
}

RStructuralPiezoLayeredDynamicsCheckpoint RStructuralPiezoLayeredDynamics::checkpoint(const RStructuralPiezoLayeredDynamicsState &state) const {validate(state); return state;}
RStructuralPiezoLayeredDynamicsState RStructuralPiezoLayeredDynamics::restore(const RStructuralPiezoLayeredDynamicsCheckpoint &checkpoint) const {validate(checkpoint); return checkpoint;}
