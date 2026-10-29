#include "rstructural_piezo_axial_operators.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok,const char *message)
{ if (!ok) throw std::invalid_argument(message); }
double finite(double x)
{ require(std::isfinite(x),"Axial operator is nonfinite or unrepresentable"); return x; }
std::vector<double> zeros(std::size_t rows,std::size_t columns)
{
    require(columns==0 || rows<=std::numeric_limits<std::size_t>::max()/columns,
            "Axial operator dimensions overflow");
    return std::vector<double>(rows*columns,0.0);
}
RStructuralPiezoAxialOperatorBlock block(std::size_t n,std::size_t p,std::size_t u,std::size_t y)
{
    RStructuralPiezoAxialOperatorBlock b; b.dofs=n; b.ports=p; b.inputs=u; b.outputs=y;
    b.M=b.C=b.K=b.Kshort=b.Kmaterial=b.Kgeometric=b.KfixedCharge=zeros(n,n);
    b.G=zeros(n,p); b.Ce=zeros(p,p); b.B=zeros(n,u);
    b.Lq=b.Lv=b.La=zeros(y,n); b.D=zeros(y,u); return b;
}
std::vector<double> select(const std::vector<double> &x,const std::vector<std::size_t> &map)
{
    std::vector<double> result; result.reserve(map.size());
    for(auto d:map) result.push_back(x.at(d)); return result;
}
void add(double &target,double x) { target=finite(target+finite(x)); }
// Evaluate the complete product/ratio before materializing its binary64
// exponent. Intermediate underflow/overflow must not erase a representable
// cross term (MSVC long double has the same range as double).
double productRatio(std::initializer_list<double> factors,double divisor,const char *context)
{
    require(finite(divisor)>0,"Axial scaled product divisor must be positive");
    bool zero=false;
    for(double factor:factors) { finite(factor); zero=zero || factor==0; }
    if(zero) return 0;
    double fraction=1; int exponent=0;
    for(double factor:factors) {
        int shift=0; fraction*=std::frexp(factor,&shift); exponent+=shift;
    }
    int shift=0; fraction/=std::frexp(divisor,&shift); exponent-=shift;
    const double value=std::scalbn(fraction,exponent);
    require(std::isfinite(value) && value!=0,context);
    return value;
}
bool decompositionMatches(double native,double first,double second,double contributionScale,std::size_t members)
{
    // Per-entry backward bound includes the contributing element terms, so
    // cancellation has a defensible bound without borrowing a large diagonal
    // from another physical subspace. Normalize before computing the error:
    // epsilon*scale itself may underflow for resolvable subnormal entries.
    const double scale=std::max({std::fabs(native),std::fabs(first),std::fabs(second),contributionScale});
    if(scale==0) return native==0 && first==0 && second==0;
    const double tolerance=1024*std::numeric_limits<double>::epsilon()*double(members+2);
    return std::fabs(finite(native/scale-first/scale-second/scale))<=tolerance;
}
// Normalize by the stiffness diagonal before Cholesky. This detects mechanisms
// and unstable references without an arbitrary dimensional stiffness floor.
double certifyStiffness(const RStructuralPiezoAxialOperatorBlock &b)
{
    const auto n=b.dofs;
    require(n>0,"Axial design requires at least one free mechanical DOF");
    auto l=zeros(n,n); std::vector<double> scale(n);
    for(std::size_t i=0;i<n;++i) {
        require(b.K[i*n+i]>0,"Axial working point has a mechanism or unstable stiffness");
        scale[i]=std::sqrt(b.K[i*n+i]);
    }
    double minimum=1;
    const double tolerance=1024*std::numeric_limits<double>::epsilon()*double(n);
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<=i;++j) {
        double v=finite((b.K[i*n+j]/scale[i])/scale[j]);
        for(std::size_t k=0;k<j;++k) v=finite(v-l[i*n+k]*l[j*n+k]);
        if(i==j) {
            require(v>tolerance,"Axial working point has a mechanism, instability or unresolved conditioning");
            minimum=std::min(minimum,v); l[i*n+j]=std::sqrt(v);
        } else l[i*n+j]=finite(v/l[j*n+j]);
    }
    return minimum;
}
}

RStructuralPiezoAxialOperators RStructuralPiezoAxialOperators::compile(const RModel &model)
{
    const RStructuralTrussDynamicsModelAdapter adapter(model);
    const auto &readout=adapter.piezoReadout();
    const auto &ref=readout.reference(); const auto &s=adapter.system(); const auto &maps=adapter.maps();
    require(ref.equilibriumCertified || (!ref.explicitWorkpoint &&
            std::all_of(ref.displacements.begin(),ref.displacements.end(),[](double q){return q==0;}) &&
            ref.freeForceResidualNorm<=ref.forceTolerance),
            "Axial operators require a certified native working point or equilibrated zero source reference");
    require(s.prescribedMotions.empty(),"Axial stationary design does not support moving prescribed supports");
    RStructuralPiezoAxialOperators result;
    result.sourceFingerprint=adapter.sourceFingerprint(); result.nativeSystemIdentity=adapter.driver().systemIdentity();
    result.readoutIdentity=readout.setup().identityBytes(); result.reference=ref; result.nativeSystem=s;
    result.reference.equilibriumCertified=true;
    result.certificateAuthority=ref.explicitWorkpoint?QStringLiteral("native-explicit-workpoint-certificate"):
        QStringLiteral("native-source-reference-with-operator-compiler-equilibrium-check");
    result.orderedMemberPhysicalIds=maps.orderedMemberPhysicalIds; result.orderedPortIds=maps.orderedPortIds;
    result.q0=ref.displacements; result.f0=ref.inputs.externalForces;
    result.internalForce0=ref.assembly.internalForces; result.reactions0=result.internalForce0;
    for(std::size_t d=0;d<result.q0.size();++d) {
        result.reactions0[d]=finite(result.reactions0[d]-result.f0[d]);
        RStructuralPiezoAxialDof dof{maps.orderedPhysicalNodeIds.at(d/3),uint(d%3),d,s.prescribedMask[d]};
        result.fullDofs.push_back(dof);
        if(!dof.prescribed) { result.freeDofs.push_back(dof); result.freeToFull.push_back(d); }
    }
    const auto p=result.orderedPortIds.size(),n=result.q0.size();
    std::map<QString,std::size_t> portColumns;
    for(std::size_t j=0;j<p;++j) portColumns.emplace(result.orderedPortIds[j],j);
    std::map<QString,const RStructuralPiezoReadoutPort*> catalog;
    for(const auto &port:readout.setup().ports) catalog.emplace(port.portId,&port);
    require(catalog.size()==p,"Axial operators require every isolated port to be classified");
    std::map<QString,std::size_t> electricalInputs;
    result.V0.resize(p); result.Q0.resize(p);
    for(std::size_t j=0;j<p;++j) {
        const auto &id=result.orderedPortIds[j]; const auto m=maps.memberByPortId.at(id);
        const auto &e=*ref.assembly.memberResponses.at(m).electrical;
        result.V0[j]=e.voltage; result.Q0[j]=e.conductorCharge;
        if(catalog.at(id)->role!=RStructuralPiezoPortRole::Sensor) {
            require(std::holds_alternative<RStructuralTrussVoltage>(ref.inputs.memberInputs[m]),
                    "Current classified native catalog supports voltage actuators only");
            electricalInputs.emplace(id,result.inputs.size());
            result.inputs.push_back({id,QStringLiteral("V"),0,0,m});
        }
    }
    for(const auto &d:result.freeDofs)
        result.inputs.push_back({QString(),QStringLiteral("N"),d.physicalNodeId,d.axis,0});
    const auto referenceState=adapter.driver().initialize(ref.time,ref.displacements,std::vector<double>(n,0),ref.inputs);
    const auto samples=readout.observe(referenceState,ref.assembly);
    for(const auto &sample:samples)
        result.outputs.push_back({sample.portId,sample.unit,sample.mode,sample.actionRecord,sample.value});
    result.full=block(n,p,result.inputs.size(),result.outputs.size()); auto &b=result.full;
    auto shortContributionScale=zeros(n,n),chargeContributionScale=zeros(n,n);
    for(std::size_t d=0;d<n;++d) b.M[d*n+d]=s.nodalLumpedMass.at(d/3);
    // Decomposition uses the actual chord/force and uncondensed material law.
    // The independently compiled native tangent below remains authoritative.
    for(std::size_t m=0;m<s.members.size();++m) {
        const auto &member=s.members[m]; const auto &response=ref.assembly.memberResponses[m];
        RStructuralTrussVector3 direction{};
        for(std::size_t a=0;a<3;++a) direction[a]=finite((finite(s.referenceNodes[member.node2][a]-s.referenceNodes[member.node1][a])+
            finite(ref.displacements[3*member.node2+a]-ref.displacements[3*member.node1+a]))/response.length);
        double materialTangent=response.axialTangent;
        if(const auto *law=std::get_if<RStructuralTrussPiezoLaw>(&member.law)) materialTangent=law->axialStiffness;
        if(const auto *law=std::get_if<RStructuralTrussCableLaw>(&member.law)) {
            const auto l0=std::get<RStructuralTrussCableLength>(ref.inputs.memberInputs[m]).value;
            const double trial=finite(law->initialTension+finite((law->axialRigidity/l0)*finite(response.length-l0)));
            const double scale=std::max({std::fabs(law->initialTension),std::fabs((law->axialRigidity/l0)*(response.length-l0)),law->axialRigidity});
            require(std::fabs(trial)>128*std::numeric_limits<double>::epsilon()*scale,
                    "Axial working point lies at an unresolved cable slack/taut kink");
        }
        for(std::size_t i=0;i<6;++i) for(std::size_t j=i;j<6;++j) {
            const auto row=3*(i<3?member.node1:member.node2)+i%3,column=3*(j<3?member.node1:member.node2)+j%3;
            const double sign=(i<3)==(j<3)?1.0:-1.0,nn=direction[i%3]*direction[j%3];
            const double material=productRatio({sign,materialTangent,direction[i%3],direction[j%3]},1,
                "Axial material contribution is outside binary64 representability");
            const double geometric=i%3==j%3?
                productRatio({sign,response.axialForce,1-nn},response.length,
                    "Axial geometric contribution is outside binary64 representability"):
                productRatio({-sign,response.axialForce,direction[i%3],direction[j%3]},response.length,
                    "Axial geometric contribution is outside binary64 representability");
            const double scale=std::max(std::fabs(material),std::fabs(geometric));
            add(b.Kmaterial[row*n+column],material); add(b.Kgeometric[row*n+column],geometric);
            shortContributionScale[row*n+column]=std::max(shortContributionScale[row*n+column],scale);
            if(row!=column) {
                add(b.Kmaterial[column*n+row],material); add(b.Kgeometric[column*n+row],geometric);
                shortContributionScale[column*n+row]=std::max(shortContributionScale[column*n+row],scale);
            }
        }
    }
    auto shortInputs=ref.inputs.memberInputs;
    for(std::size_t j=0;j<p;++j) {
        const auto m=maps.memberByPortId.at(result.orderedPortIds[j]); const auto &member=s.members[m];
        const auto &law=std::get<RStructuralTrussPiezoLaw>(member.law);
        const auto &response=ref.assembly.memberResponses[m];
        b.Ce[j*p+j]=law.capacitance; shortInputs[m]=RStructuralTrussVoltage{result.V0[j]};
        for(std::size_t a=0;a<3;++a) {
            const double direction=finite((finite(s.referenceNodes[member.node2][a]-s.referenceNodes[member.node1][a])+
                finite(ref.displacements[3*member.node2+a]-ref.displacements[3*member.node1+a]))/response.length);
            const double coupling=productRatio({law.coupling,direction},1,
                "Axial coupling contribution is outside binary64 representability");
            b.G[(3*member.node1+a)*p+j]=-coupling;
            b.G[(3*member.node2+a)*p+j]=coupling;
        }
        if(std::holds_alternative<RStructuralTrussCharge>(ref.inputs.memberInputs[m]))
            for(std::size_t i=0;i<n;++i) for(std::size_t k=i;k<n;++k) {
                const double condensed=productRatio({b.G[i*p+j],b.G[k*p+j]},law.capacitance,
                    "Axial fixed-charge condensation is outside binary64 representability");
                add(b.KfixedCharge[i*n+k],condensed);
                chargeContributionScale[i*n+k]=std::max(chargeContributionScale[i*n+k],std::fabs(condensed));
            }
        const auto input=electricalInputs.find(result.orderedPortIds[j]);
        if(input!=electricalInputs.end()) {
            for(std::size_t i=0;i<n;++i) b.B[i*b.inputs+input->second]=finite(-b.G[i*p+j]);
        }
    }
    const auto shortAssembly=RStructuralTrussAssembly::evaluate(s.referenceNodes,ref.displacements,s.members,shortInputs);
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j) {
        b.K[i*n+j]=ref.assembly.tangent.findValue(uint(i),uint(j));
        b.Kshort[i*n+j]=shortAssembly.tangent.findValue(uint(i),uint(j));
        if(i>j) {
            b.KfixedCharge[i*n+j]=b.KfixedCharge[j*n+i];
            chargeContributionScale[i*n+j]=chargeContributionScale[j*n+i];
        }
    }
    for(std::size_t i=0;i<n*n;++i) {
        require(decompositionMatches(b.Kshort[i],b.Kmaterial[i],b.Kgeometric[i],shortContributionScale[i],s.members.size()),
                "Native axial short tangent cannot resolve the material/geometric entry at its local scale");
        require(decompositionMatches(b.K[i],b.Kshort[i],b.KfixedCharge[i],chargeContributionScale[i],s.members.size()),
                "Native axial tangent cannot resolve the fixed-charge entry at its local scale");
    }
    for(std::size_t j=0;j<result.freeToFull.size();++j)
        b.B[result.freeToFull[j]*b.inputs+electricalInputs.size()+j]=1;
    for(std::size_t row=0;row<result.outputs.size();++row) {
        const auto &output=result.outputs[row]; const auto j=portColumns.at(output.portId);
        const auto m=maps.memberByPortId.at(output.portId); const auto &port=*catalog.at(output.portId);
        const bool voltage=std::holds_alternative<RStructuralTrussVoltage>(ref.inputs.memberInputs[m]);
        double qFactor=0,qDivisor=1,direct=0;
        switch(output.mode) {
        case RStructuralPiezoReadoutMode::None:
        case RStructuralPiezoReadoutMode::OpenCircuitVoltage:
            qFactor=voltage?0:1; qDivisor=voltage?1:b.Ce[j*p+j]; direct=voltage?1:0; break;
        case RStructuralPiezoReadoutMode::ShortCircuitCharge:
        case RStructuralPiezoReadoutMode::ActuatorCharge:
            qFactor=voltage?-1:0; direct=voltage?b.Ce[j*p+j]:1; break;
        case RStructuralPiezoReadoutMode::IdealChargeAmplifier:
            qFactor=(voltage?-1:0)*double(port.amplifierPolarity); qDivisor=port.Cf; break;
        }
        for(std::size_t d=0;d<n;++d) b.Lq[row*n+d]=productRatio({qFactor,b.G[d*p+j]},qDivisor,
            "Axial readout contribution is outside binary64 representability");
        const auto input=electricalInputs.find(output.portId);
        if(input!=electricalInputs.end()) b.D[row*b.inputs+input->second]=finite(direct);
    }
    auto &f=result.free; const auto nf=result.freeToFull.size(); f=block(nf,p,b.inputs,b.outputs);
    f.Ce=b.Ce; f.D=b.D;
    for(std::size_t i=0;i<nf;++i) {
        const auto d=result.freeToFull[i];
        for(std::size_t j=0;j<nf;++j) {
            const auto e=result.freeToFull[j];
            for(const auto matrices:{std::make_pair(&b.M,&f.M),std::make_pair(&b.C,&f.C),std::make_pair(&b.K,&f.K),
                std::make_pair(&b.Kshort,&f.Kshort),std::make_pair(&b.Kmaterial,&f.Kmaterial),
                std::make_pair(&b.Kgeometric,&f.Kgeometric),std::make_pair(&b.KfixedCharge,&f.KfixedCharge)})
                (*matrices.second)[i*nf+j]=(*matrices.first)[d*n+e];
        }
        require(f.M[i*nf+i]>0,"Axial free DOF has no positive native mass");
        for(std::size_t j=0;j<p;++j) f.G[i*p+j]=b.G[d*p+j];
        for(std::size_t j=0;j<b.inputs;++j) f.B[i*b.inputs+j]=b.B[d*b.inputs+j];
        for(std::size_t j=0;j<b.outputs;++j) f.Lq[j*nf+i]=b.Lq[j*n+d];
    }
    result.minimumNormalizedStiffnessPivot=certifyStiffness(f);
    const auto initial=adapter.sourceInitialState();
    result.initialDeltaQ.resize(n); result.initialVelocity=initial.velocities; result.initialAcceleration=initial.accelerations;
    for(std::size_t d=0;d<n;++d) result.initialDeltaQ[d]=finite(initial.displacements[d]-result.q0[d]);
    result.freeInitialDeltaQ=select(result.initialDeltaQ,result.freeToFull);
    result.freeInitialVelocity=select(result.initialVelocity,result.freeToFull);
    result.freeInitialAcceleration=select(result.initialAcceleration,result.freeToFull);
    return result;
}
