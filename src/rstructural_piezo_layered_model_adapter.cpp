#include "rstructural_piezo_layered_model_adapter.h"
#include <rml_model.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace {
using V=std::vector<double>;
constexpr double tol=1e-10;
void require(bool ok,const char *message) {if(!ok) throw std::invalid_argument(message);}
void finite(const V &a) {for(double x:a) require(std::isfinite(x),"Non-finite layered compiler value");}
V mv(const V &a,const V &x,size_t rows) {
    require(a.size()==rows*x.size(),"Matrix-vector dimension mismatch");
    V y(rows); for(size_t i=0;i<rows;++i) for(size_t j=0;j<x.size();++j) y[i]+=a[i*x.size()+j]*x[j]; finite(y); return y;
}
V transposeProduct(const V &t,const V &x,size_t r) {
    const size_t n=x.size(); require(t.size()==n*r,"Projection dimension mismatch");
    V y(r); for(size_t i=0;i<n;++i) for(size_t j=0;j<r;++j) y[j]+=t[i*r+j]*x[i]; finite(y); return y;
}
V congruence(const V &a,const V &t,size_t n,size_t r) {
    require(a.size()==n*n && t.size()==n*r,"Congruence dimension mismatch");
    V at(n*r),b(r*r);
    for(size_t i=0;i<n;++i) for(size_t j=0;j<r;++j) for(size_t k=0;k<n;++k) at[i*r+j]+=a[i*n+k]*t[k*r+j];
    for(size_t i=0;i<r;++i) for(size_t j=0;j<r;++j) for(size_t k=0;k<n;++k) b[i*r+j]+=t[k*r+i]*at[k*r+j];
    finite(b); return b;
}
void symmetrizeSolidRoundoff(V &matrix,size_t n,const char *message) {
    require(matrix.size()==n*n,"Solid reduced operator dimension mismatch");
    for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j){
        const double a=matrix[i*n+j],b=matrix[j*n+i];
        const double scale=std::max({std::abs(a),std::abs(b),std::sqrt(std::abs(matrix[i*n+i]))*std::sqrt(std::abs(matrix[j*n+j]))});
        const double halfDifference=std::abs(.5*a-.5*b);
        const double defect=scale==0?(halfDifference==0?0:std::numeric_limits<double>::infinity()):halfDifference/std::max(.5*scale,std::numeric_limits<double>::min());
        require(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(matrix[i*n+i])&&std::isfinite(matrix[j*n+j])&&defect<=1e-9,message);
        matrix[i*n+j]=matrix[j*n+i]=.5*a+.5*b;
    }
}
//! Dimensionless diagonal equilibration prevents comparing bending and axial SI
//! pivots against one dimensional absolute floor. Unresolved conditioning rejects.
struct PositiveFactor {
    size_t n; V scale,L;
    PositiveFactor(const V &a,size_t count,bool semidefinite=false):n(count),scale(n),L(n*n) {
        require(a.size()==n*n,"Positive factor dimension mismatch"); finite(a);
        const double threshold=128*std::numeric_limits<double>::epsilon()*std::max(size_t(1),n);
        for(size_t i=0;i<n;++i) {
            require(a[i*n+i]>=0,"Negative operator diagonal");
            require(semidefinite || a[i*n+i]>0,"Nonunique or unstable free operator");
            if(a[i*n+i]==0) for(size_t j=0;j<n;++j)
                require(a[i*n+j]==0 && a[j*n+i]==0,"PSD operator has nonzero zero-diagonal row");
            scale[i]=a[i*n+i]>0?std::sqrt(a[i*n+i]):1;
        }
        for(size_t i=0;i<n;++i) for(size_t j=0;j<=i;++j) {
            double large=std::max(scale[i],scale[j]),small=std::min(scale[i],scale[j]);
            double s=a[i*n+j]/large/small,other=a[j*n+i]/large/small;
            require(std::isfinite(s) && std::isfinite(other),"Unrepresentable equilibrated operator");
            require(std::abs(s-other)<=tol*std::max({1.,std::abs(s),std::abs(other)}),"Operator is not symmetric");
            s=(s+other)/2;
            for(size_t k=0;k<j;++k) s-=L[i*n+k]*L[j*n+k];
            if(i==j) {
                require(s>=-threshold,"Operator is indefinite");
                if(semidefinite && s<=threshold) L[i*n+i]=0;
                else {require(s>threshold,"Mechanism or unresolved free operator conditioning"); L[i*n+i]=std::sqrt(s);}
            } else if(L[j*n+j]>0) L[i*n+j]=s/L[j*n+j];
            else require(std::abs(s)<=threshold,"Indefinite semidefinite operator");
        }
        finite(L);
    }
    V solve(const V &b) const {
        require(b.size()==n,"Solve dimension mismatch"); V x(n);
        for(size_t i=0;i<n;++i) {x[i]=b[i]/scale[i]; for(size_t j=0;j<i;++j) x[i]-=L[i*n+j]*x[j]; x[i]/=L[i*n+i];}
        for(size_t ii=n;ii>0;--ii) {size_t i=ii-1; for(size_t j=i+1;j<n;++j) x[i]-=L[j*n+i]*x[j]; x[i]/=L[i*n+i];}
        for(size_t i=0;i<n;++i) x[i]/=scale[i]; finite(x); return x;
    }
};
size_t nodeIndex(const RStructuralPiezoLayeredMeshOperators &mesh,uint id) {
    auto it=std::lower_bound(mesh.nodeIds.begin(),mesh.nodeIds.end(),id);
    require(it!=mesh.nodeIds.end() && *it==id,"Node field is outside compiled source mesh"); return size_t(it-mesh.nodeIds.begin());
}
V field(const RStructuralPiezoLayeredMeshOperators &mesh,const std::vector<RStructuralPiezoLayeredNodeField> &nodes,RStructuralPiezoPhysicalDomain domain) {
    V physical(6*mesh.nodeIds.size()); std::set<uint> seen;
    for(const auto &node:nodes) {require(seen.insert(node.physicalNodeId).second,"Duplicate source node field"); size_t i=nodeIndex(mesh,node.physicalNodeId);
        std::copy(node.values.begin(),node.values.end(),physical.begin()+6*i);}
    return mesh.projectPhysical(physical,domain);
}
//! Affine RREF of global-axis rows of P. Unsupported homogeneous constraints
//! are redundant, but unsupported nonzero prescriptions are inconsistent.
void constraints(const RStructuralPiezoLayeredMeshOperators &mesh,const RStructuralPiezoLayeredSetup &s,RStructuralPiezoLayeredCompiledData &d) {
    size_t n=d.fullDofs,m=s.supports.size(); V a(m*(n+1)); std::vector<size_t> pivots;
    for(size_t i=0;i<m;++i) {
        size_t row=6*nodeIndex(mesh,s.supports[i].physicalNodeId)+size_t(s.supports[i].dof);
        require(size_t(s.supports[i].dof)<6,"Unknown physical support DOF");
        double scale=0; for(size_t j=0;j<n;++j) scale=std::hypot(scale,mesh.physicalTransform[row*n+j]);
        if(scale==0) {require(s.supports[i].value==0,"Unsupported nonzero fixed DOF"); continue;}
        for(size_t j=0;j<n;++j) a[i*(n+1)+j]=mesh.physicalTransform[row*n+j]/scale;
        a[i*(n+1)+n]=s.supports[i].value/scale;
    }
    finite(a); size_t rank=0;
    for(size_t j=0;j<n && rank<m;++j) {
        size_t best=rank; for(size_t i=rank;i<m;++i) if(std::abs(a[i*(n+1)+j])>std::abs(a[best*(n+1)+j])) best=i;
        if(std::abs(a[best*(n+1)+j])<=tol) continue;
        for(size_t k=0;k<=n;++k) std::swap(a[rank*(n+1)+k],a[best*(n+1)+k]);
        double pivot=a[rank*(n+1)+j]; for(size_t k=0;k<=n;++k) a[rank*(n+1)+k]/=pivot;
        for(size_t i=0;i<m;++i) if(i!=rank) {double factor=a[i*(n+1)+j]; for(size_t k=0;k<=n;++k) a[i*(n+1)+k]-=factor*a[rank*(n+1)+k];}
        pivots.push_back(j); ++rank;
    }
    double rhsScale=0; for(const auto &support:s.supports) rhsScale=std::max(rhsScale,std::abs(support.value));
    for(size_t i=rank;i<m;++i) require(std::abs(a[i*(n+1)+n])<=tol*rhsScale,"Inconsistent affine supports");
    d.constraintRank=uint(rank); d.redundantConstraints=uint(m-rank); d.freeDofs=uint(n-rank); d.p.assign(n,0);
    for(size_t i=0;i<rank;++i) d.p[pivots[i]]=a[i*(n+1)+n];
    const size_t r=d.freeDofs; d.T.assign(n*r,0); size_t column=0;
    for(size_t j=0;j<n;++j) if(std::find(pivots.begin(),pivots.end(),j)==pivots.end()) {
        V v(n); v[j]=1; for(size_t i=0;i<rank;++i) v[pivots[i]]=-a[i*(n+1)+j];
        // Twice reorthogonalized modified Gram-Schmidt.
        for(size_t pass=0;pass<2;++pass) for(size_t k=0;k<column;++k) {double dot=0; for(size_t i=0;i<n;++i) dot+=d.T[i*r+k]*v[i]; for(size_t i=0;i<n;++i) v[i]-=dot*d.T[i*r+k];}
        double norm=0; for(double x:v) norm=std::hypot(norm,x); require(norm>tol,"Constraint nullspace lost rank");
        for(size_t i=0;i<n;++i) d.T[i*r+column]=v[i]/norm; ++column;
    }
    finite(d.T); finite(d.p);
}
V reducedState(const RStructuralPiezoLayeredCompiledData &d,const V &q,bool velocity=false) {
    require(q.size()==d.fullDofs,"State dimension mismatch"); V delta=q;
    if(!velocity) for(size_t i=0;i<delta.size();++i) delta[i]-=d.p[i];
    V z=transposeProduct(d.T,delta,d.freeDofs),back=mv(d.T,z,d.fullDofs);
    double norm=0,residual=0; for(size_t i=0;i<q.size();++i) {norm=std::hypot(norm,std::max(std::abs(q[i]),velocity?0.:std::abs(d.p[i]))); residual=std::hypot(residual,back[i]-delta[i]);}
    require(residual<=tol*norm,"Absolute initial/reference state violates fixed supports"); return z;
}
double residualCertificate(const V &t,const V &residual,const V &scale,size_t r) {
    V projected=transposeProduct(t,residual,r); double maximum=0;
    for(size_t j=0;j<r;++j) {double bound=0; for(size_t i=0;i<residual.size();++i) bound+=std::abs(t[i*r+j])*scale[i];
        require(std::abs(projected[j])<=tol*bound,"Reference force balance failed");
        if(bound>0) maximum=std::max(maximum,std::abs(projected[j])/bound);}
    return maximum;
}
}

RStructuralPiezoLayeredModelAdapter::RStructuralPiezoLayeredModelAdapter(const RModel &model,const RStructuralPiezoLayeredSetup &s,const V &modalC) {
    require(s.enabled,"Layered source is disabled"); s.validate(model);
    const auto noLegacyConditions=[](const RElementGroup &group) {
        require(group.getNBoundaryConditions()==0 && group.getNInitialConditions()==0 && group.getNEnvironmentConditions()==0,
                "Migrate legacy group conditions to the authoritative layered setup before compilation");
    };
    for(uint i=0;i<model.getNPoints();++i) noLegacyConditions(model.getPoint(i));
    for(uint i=0;i<model.getNLines();++i) noLegacyConditions(model.getLine(i));
    for(uint i=0;i<model.getNSurfaces();++i) noLegacyConditions(model.getSurface(i));
    for(uint i=0;i<model.getNVolumes();++i) noLegacyConditions(model.getVolume(i));
    const QString fingerprint=s.computeSourceFingerprint(model);
    require(s.sourceFingerprint.isEmpty() || s.sourceFingerprint==fingerprint,"Stale layered source fingerprint");
    const auto family=s.family==RStructuralPiezoLayeredModelFamily::PlanarBeam?RStructuralPiezoLayeredFamily::Beam:
        (s.family==RStructuralPiezoLayeredModelFamily::Solid3D?RStructuralPiezoLayeredFamily::Solid:RStructuralPiezoLayeredFamily::Plate);
    std::map<QString,RStructuralPiezoContinuumMaterial> materials;
    for(const auto &source:s.materials) {
        RStructuralPiezoContinuumMaterial m;
        if(source.form==RStructuralPiezoMaterialForm::StrainCharge && source.piezoelectric) {
            m=RStructuralPiezoContinuumMaterial::fromStrainCharge(source.elasticityOrCompliance,source.coupling,source.permittivity,source.density);
            m.piezoelectric=source.piezoelectric;
        } else if(source.form==RStructuralPiezoMaterialForm::StrainCharge) {
            PositiveFactor compliance(V(source.elasticityOrCompliance.begin(),source.elasticityOrCompliance.end()),6);
            for(size_t j=0;j<6;++j) {V unit(6); unit[j]=1; const auto column=compliance.solve(unit); for(size_t i=0;i<6;++i) m.elasticity[i*6+j]=column[i];}
            m.density=source.density;
        } else {m.elasticity=source.elasticityOrCompliance; m.coupling=source.coupling; m.permittivity=source.permittivity; m.density=source.density; m.piezoelectric=source.piezoelectric;}
        m.validate(); require(materials.emplace(source.materialId,m).second,"Duplicate material identity");
    }
    std::map<uint,RStructuralPiezoLayeredNode> nodes; std::vector<RStructuralPiezoLayeredMember> members;
    for(const auto &source:s.members) {
        uint position=model.findElementPositionByPhysicalId(source.physicalElementId); require(position<model.getNElements(),"Missing physical element");
        const auto &e=model.getElement(position);
        require(family==RStructuralPiezoLayeredFamily::Solid?
                (e.getType()==R_ELEMENT_TETRA1||e.getType()==R_ELEMENT_HEXA1):
                e.getType()==(family==RStructuralPiezoLayeredFamily::Beam?R_ELEMENT_BEAM1:R_ELEMENT_QUAD1),"Source element family mismatch");
        RStructuralPiezoLayeredMember member; member.physicalId=source.physicalElementId;
        for(uint j=0;j<RElement::getNNodes(e.getType());++j) {
            uint np=e.getNodeId(j); require(np<model.getNNodes(),"Invalid source connectivity"); const auto &node=model.getNode(np); uint id=model.getNodePhysicalId(np);
            nodes[id]={id,{{node.getX(),node.getY(),node.getZ()}}}; member.nodeIds.push_back(id);
        }
        for(const auto &sourceLayer:source.layers) {
            auto found=materials.find(sourceLayer.materialId); require(found!=materials.end(),"Missing source layer material");
            RStructuralPiezoLayer l; l.material=found->second; l.lowerZ=sourceLayer.lowerZ; l.upperZ=sourceLayer.upperZ;
            l.width=sourceLayer.width; l.portId=sourceLayer.portId; l.polarity=sourceLayer.polarity; l.coverage=sourceLayer.coverage; member.layers.push_back(l);
        }
        members.push_back(member);
    }
    std::vector<RStructuralPiezoLayeredNode> ordered; for(const auto &node:nodes) ordered.push_back(node.second);
    require(ordered.size()<=2048/(family==RStructuralPiezoLayeredFamily::Plate?5u:3u),"Dense layered compiler exceeds its 2048 generalized-DOF admission bound");
    RStructuralPiezoPlanarFrame frame; frame.origin=s.frame.origin; frame.ex=s.frame.ex; frame.ey=s.frame.ey; frame.ez=s.frame.ez;
    meshOperators=RStructuralPiezoLayeredMesh::assemble(family,frame,ordered,members);
    auto &d=compiled; const auto &op=meshOperators.operators; size_t n=op.dofs; d.fullDofs=op.dofs;
    d.sourceFingerprint=fingerprint; d.materialProvenance=s.materials; d.portMetadata=s.ports;
    constraints(meshOperators,s,d); const size_t r=d.freeDofs;
    std::map<QString,const RStructuralPiezoLayeredPort*> ports; for(const auto &port:s.ports) require(ports.emplace(port.readout.portId,&port).second,"Duplicate port identity");
    std::vector<RStructuralPiezoElectricalBoundary> boundaries;
    for(const auto &id:op.portIds) {
        auto it=ports.find(id); require(it!=ports.end(),"Unclassified native electrode"); const auto &source=*it->second;
        RStructuralPiezoElectricalBoundary b; b.portId=id;
        if(source.readout.role!=RStructuralPiezoPortRole::Sensor) {b.role=RStructuralPiezoElectricalRole::VoltageActuator; d.actuatorVoltages0.push_back(source.initialValue);}
        else if(source.readout.readout==RStructuralPiezoReadoutMode::OpenCircuitVoltage) {b.role=RStructuralPiezoElectricalRole::FixedChargeOpenCircuit; b.fixedCharge=source.initialValue;}
        else if(source.readout.readout==RStructuralPiezoReadoutMode::ShortCircuitCharge) b.role=RStructuralPiezoElectricalRole::ShortCircuitChargeSensor;
        else {require(source.readout.readout==RStructuralPiezoReadoutMode::IdealChargeAmplifier,"Unsupported sensor boundary"); b.role=RStructuralPiezoElectricalRole::IdealChargeAmplifier; b.feedbackCapacitance=source.readout.Cf; b.referenceCharge=source.readout.Qref;}
        boundaries.push_back(b);
    }
    require(ports.size()==op.portIds.size(),"Declared electrode is absent from source mesh");
    d.backgroundForce=field(meshOperators,s.backgroundForces,RStructuralPiezoPhysicalDomain::ForceAndMoment);
    d.referenceForce=d.backgroundForce;
    std::map<QString,RStructuralPiezoPortWorkpoint> workports;
    if(s.workpoint.available) {
        std::vector<RStructuralPiezoLayeredNodeField> q,f;
        for(const auto &node:s.workpoint.nodes) {q.push_back({node.physicalNodeId,node.q0}); f.push_back({node.physicalNodeId,node.f0});}
        require(q.size()==nodes.size(),"Reference must cover every source node");
        d.q0=field(meshOperators,q,RStructuralPiezoPhysicalDomain::Displacement); d.referenceForce=field(meshOperators,f,RStructuralPiezoPhysicalDomain::ForceAndMoment); reducedState(d,d.q0);
        for(const auto &port:s.workpoint.ports) require(workports.emplace(port.portId,port).second,"Duplicate reference port");
        require(workports.size()==op.portIds.size(),"Reference must cover every source port");
        size_t ai=0; for(const auto &b:boundaries) {require(workports.count(b.portId),"Missing reference port"); if(b.role==RStructuralPiezoElectricalRole::VoltageActuator) d.actuatorVoltages0[ai++]=workports.at(b.portId).V0;}
    } else d.q0=d.p;
    RStructuralPiezoElectricalSystem provisional(op,boundaries,d.q0,d.actuatorVoltages0);
    const auto &effective=provisional.condensed();
    d.mass=congruence(op.mass,d.T,n,r); d.stiffness=congruence(effective.stiffness,d.T,n,r);
    if(family==RStructuralPiezoLayeredFamily::Solid){
        symmetrizeSolidRoundoff(d.mass,r,"Solid reduced mass symmetry defect exceeds roundoff tolerance");
        symmetrizeSolidRoundoff(d.stiffness,r,"Solid reduced stiffness symmetry defect exceeds roundoff tolerance");
    }
    PositiveFactor stiffnessFactor(d.stiffness,r),massFactor(d.mass,r);
    size_t na=d.actuatorVoltages0.size(); d.input.assign(r*na,0);
    for(size_t i=0;i<r;++i) for(size_t j=0;j<na;++j) for(size_t k=0;k<n;++k) d.input[i*na+j]-=d.T[k*r+i]*effective.coupling[k*na+j];
    if(!s.workpoint.available) {
        V rhs=d.referenceForce,kp=mv(effective.stiffness,d.p,n),gv=mv(effective.coupling,d.actuatorVoltages0,n);
        for(size_t i=0;i<n;++i) rhs[i]-=kp[i]+gv[i]+effective.fixedChargeForce[i];
        V z=stiffnessFactor.solve(transposeProduct(d.T,rhs,r)); d.q0=mv(d.T,z,n); for(size_t i=0;i<n;++i) d.q0[i]+=d.p[i];
    }
    electricalSystem.reset(new RStructuralPiezoElectricalSystem(op,boundaries,d.q0,d.actuatorVoltages0));
    const auto &response=electricalSystem->referenceResponse(); V residual(n),scale(n);
    for(size_t i=0;i<n;++i) {
        residual[i]=response.internalForce[i]-d.referenceForce[i]; scale[i]=std::abs(d.referenceForce[i]);
        for(size_t j=0;j<n;++j) scale[i]+=std::abs(op.stiffness[i*n+j]*d.q0[j]);
        for(size_t j=0;j<op.portIds.size();++j) scale[i]+=std::abs(op.coupling[i*op.portIds.size()+j]*response.ports[j].voltage);
    }
    d.forceResidual=residualCertificate(d.T,residual,scale,r); d.supportReactions=residual;
    // Open-circuit response stores the prescribed charge. Independently
    // recompute Ce*V-G^T*q so the certificate also checks the actual FE equation.
    const size_t np=op.portIds.size();
    for(size_t j=0;j<np;++j) {
        double charge=0,bound=std::abs(response.ports[j].charge);
        for(size_t k=0;k<np;++k) {double term=op.capacitance[j*np+k]*response.ports[k].voltage; charge+=term; bound+=std::abs(term);}
        for(size_t i=0;i<n;++i) {double term=op.coupling[i*np+j]*d.q0[i]; charge-=term; bound+=std::abs(term);}
        double error=std::abs(charge-response.ports[j].charge);
        require(std::isfinite(charge) && std::isfinite(bound) && error<=tol*bound,"Native reference electrical equation residual failed");
        if(bound>0) d.electricResidual=std::max(d.electricResidual,error/bound);
    }
    if(s.workpoint.available) for(const auto &port:response.ports) {
        const auto &given=workports.at(port.portId);
        double vscale=std::max(std::abs(given.V0),std::abs(port.voltage)),qscale=std::max(std::abs(given.Q0),std::abs(port.charge));
        double vr=std::abs(given.V0-port.voltage),qr=std::abs(given.Q0-port.charge);
        require(vr<=tol*vscale && qr<=tol*qscale,"Reference electrical V0/Q0 is inconsistent with native operators");
        if(vscale>0) d.electricResidual=std::max(d.electricResidual,vr/vscale); if(qscale>0) d.electricResidual=std::max(d.electricResidual,qr/qscale);
    }
    d.fullDamping.assign(n*n,0);
    switch(s.damping.kind) {
    case RStructuralPiezoDampingKind::None: require(modalC.empty(),"Unexpected damping authority"); d.dampingAuthority="none"; break;
    case RStructuralPiezoDampingKind::Rayleigh:
        require(modalC.empty(),"Unexpected damping authority"); require(s.damping.rayleighAlpha>=0 && s.damping.rayleighBeta>=0,"Negative Rayleigh coefficient");
        for(size_t i=0;i<n*n;++i) d.fullDamping[i]=s.damping.rayleighAlpha*op.mass[i]+s.damping.rayleighBeta*effective.stiffness[i]; d.dampingAuthority="Rayleigh: full M and electrical tangent"; break;
    case RStructuralPiezoDampingKind::Imported:
        require(modalC.empty(),"Unexpected damping authority"); require(s.damping.importedDimension==n,"Imported C must use supported generalized source order"); d.fullDamping=s.damping.importedMatrix; d.dampingAuthority="imported: "+s.damping.sourceId; break;
    default:
        require(!modalC.empty(),"Modal damping requires authoritative full C from validated design backend"); d.fullDamping=modalC;
        d.dampingAuthority="modal backend supplied full C; strategy="+QString::number(uint(s.damping.kind)); break;
    }
    PositiveFactor dampingCertificate(d.fullDamping,n,true); d.damping=congruence(d.fullDamping,d.T,n,r);
    V initialQ=s.initialDisplacements.empty()?d.q0:field(meshOperators,s.initialDisplacements,RStructuralPiezoPhysicalDomain::Displacement);
    V initialV=field(meshOperators,s.initialVelocities,RStructuralPiezoPhysicalDomain::Velocity);
    d.initialZ=reducedState(d,initialQ); d.initialVelocity=reducedState(d,initialV,true);
    auto initialResponse=electricalSystem->evaluate(initialQ,d.actuatorVoltages0); V rhs=d.backgroundForce,cv=mv(d.fullDamping,initialV,n);
    for(size_t i=0;i<n;++i) rhs[i]-=initialResponse.internalForce[i]+cv[i];
    d.initialAcceleration=massFactor.solve(transposeProduct(d.T,rhs,r));
    V offset=d.backgroundForce; for(size_t i=0;i<n;++i) offset[i]-=response.internalForce[i]; d.backgroundOffset=transposeProduct(d.T,offset,r);
    size_t per=family==RStructuralPiezoLayeredFamily::Plate?5:3;
    for(size_t i=0;i<n;++i) {bool rotation=family!=RStructuralPiezoLayeredFamily::Solid&&i%per>=(per==3?2u:3u); d.generalizedCoordinateUnits.push_back(rotation?"rad":"m"); d.generalizedForceUnits.push_back(rotation?"N m":"N");}
    finite(d.input); finite(d.backgroundOffset);
}

std::vector<double> RStructuralPiezoLayeredModelAdapter::reduceForce(const V &physical) const {
    return transposeProduct(compiled.T,meshOperators.projectPhysical(physical,RStructuralPiezoPhysicalDomain::ForceAndMoment),compiled.freeDofs);
}
std::vector<double> RStructuralPiezoLayeredModelAdapter::reconstruct(const V &z) const {
    require(z.size()==compiled.freeDofs,"Reduced state dimension mismatch"); V q=mv(compiled.T,z,compiled.fullDofs);
    for(size_t i=0;i<q.size();++i) q[i]+=compiled.p[i]; return meshOperators.reconstructPhysical(q);
}
std::vector<double> RStructuralPiezoLayeredModelAdapter::reduceDisplacement(const V &physical) const {
    return reducedState(compiled,meshOperators.projectPhysical(physical,RStructuralPiezoPhysicalDomain::Displacement));
}
