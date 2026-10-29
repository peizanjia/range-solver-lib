#include "rstructural_piezo_solid_kernel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace {
using V3=std::array<double,3>;using V=std::vector<double>;
void require(bool ok,const char *s){if(!ok)throw std::invalid_argument(s);}
V3 sub(const V3&a,const V3&b){return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
double dot(const V3&a,const V3&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
V3 cross(const V3&a,const V3&b){return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};}
double norm(const V3&a){return std::hypot(std::hypot(a[0],a[1]),a[2]);}
std::vector<std::vector<size_t>> faces(size_t nn){
    if(nn==4)return {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
    return {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
}
//! Complete separating-axis test for the supported convex affine solids.
bool overlap(const std::vector<V3>&a,const std::vector<V3>&b,double tol){
    std::vector<V3> axes,ea,eb;
    for(const auto *poly:{&a,&b})for(const auto&face:faces(poly->size())){
        axes.push_back(cross(sub((*poly)[face[1]],(*poly)[face[0]]),sub((*poly)[face[2]],(*poly)[face[0]])));
        for(size_t j=0;j<face.size();++j)(poly==&a?ea:eb).push_back(sub((*poly)[face[(j+1)%face.size()]],(*poly)[face[j]]));
    }
    for(const auto&x:ea)for(const auto&y:eb)axes.push_back(cross(x,y));
    double scale=0;for(const auto *poly:{&a,&b})for(const auto&x:*poly)scale=std::max(scale,norm(sub(x,a[0])));
    for(auto axis:axes){
        const double length=norm(axis);if(length==0)continue;for(double&x:axis)x/=length;
        double loA=std::numeric_limits<double>::infinity(),hiA=-loA,loB=loA,hiB=hiA;
        for(const auto&x:a){double p=dot(axis,sub(x,a[0]));loA=std::min(loA,p);hiA=std::max(hiA,p);}
        for(const auto&x:b){double p=dot(axis,sub(x,a[0]));loB=std::min(loB,p);hiB=std::max(hiB,p);}
        if(std::min(hiA,hiB)-std::max(loA,loB)<=tol*scale)return false;
    }
    return true;
}
struct Factor{
    size_t n;V scale,l;
    Factor(const V&a,size_t size):n(size),scale(n),l(n*n){
        for(size_t i=0;i<n;++i){require(std::isfinite(a[i*n+i])&&a[i*n+i]>0,"Solid electrical domain has an ungrounded or degenerate potential");scale[i]=std::sqrt(a[i*n+i]);}
        for(size_t i=0;i<n;++i)for(size_t j=0;j<=i;++j){
            double v=a[i*n+j]/scale[i]/scale[j];for(size_t k=0;k<j;++k)v-=l[i*n+k]*l[j*n+k];
            if(i==j){require(std::isfinite(v)&&v>128*std::numeric_limits<double>::epsilon()*std::max(size_t(1),n),"Solid electrical interior is floating or unresolved");l[i*n+j]=std::sqrt(v);}
            else l[i*n+j]=v/l[j*n+j];
        }
    }
    V solve(V b)const{
        for(size_t i=0;i<n;++i){b[i]/=scale[i];for(size_t j=0;j<i;++j)b[i]-=l[i*n+j]*b[j];b[i]/=l[i*n+i];}
        for(size_t k=n;k>0;--k){size_t i=k-1;for(size_t j=i+1;j<n;++j)b[i]-=l[j*n+i]*b[j];b[i]/=l[i*n+i];}
        for(size_t i=0;i<n;++i){b[i]/=scale[i];require(std::isfinite(b[i]),"Solid electrical condensation overflow");}return b;
    }
};
struct Element{std::vector<size_t> indices;std::vector<V3>x;RStructuralPiezoSolidNodalOperators op;const RStructuralPiezoLayer*source=nullptr;};
void symmetrizeRoundoff(V &matrix,size_t n,const char *name){
    require(matrix.size()==n*n,"Invalid symmetric solid operator shape");
    for(size_t i=0;i<n;++i){
        const double diagonal=matrix[i*n+i];require(std::isfinite(diagonal),"Non-finite solid operator diagonal");
        for(size_t j=i+1;j<n;++j){
            const double a=matrix[i*n+j],b=matrix[j*n+i],otherDiagonal=matrix[j*n+j];
            require(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(otherDiagonal),"Non-finite solid operator entry");
            const double pairScale=std::max({std::abs(a),std::abs(b),std::sqrt(std::abs(diagonal))*std::sqrt(std::abs(otherDiagonal))});
            const double halfDifference=std::abs(.5*a-.5*b);
            const double defect=pairScale==0?(halfDifference==0?0:std::numeric_limits<double>::infinity()):halfDifference/std::max(.5*pairScale,std::numeric_limits<double>::min());
            require(defect<=1e-9,name);
            const double mean=.5*a+.5*b;matrix[i*n+j]=mean;matrix[j*n+i]=mean;
        }
    }
}
}

RStructuralPiezoLayeredMeshOperators rStructuralPiezoAssembleSolidMesh(
    const RStructuralPiezoPlanarFrame &frame,const std::vector<RStructuralPiezoLayeredNode>&nodes,
    const std::vector<RStructuralPiezoLayeredMember>&members,double tol){
    require(std::isfinite(tol)&&tol>=0&&tol<=1e-6,"Invalid solid relative tolerance");
    require(!nodes.empty()&&!members.empty()&&nodes.size()<=2048/3,"Solid dense admission requires 1..682 nodes and physical members");
    for(const auto *a:{&frame.origin,&frame.ex,&frame.ey,&frame.ez})for(double v:*a)require(std::isfinite(v),"Non-finite solid frame");
    require(std::abs(norm(frame.ex)-1)<=tol&&std::abs(norm(frame.ey)-1)<=tol&&std::abs(norm(frame.ez)-1)<=tol&&
            std::abs(dot(frame.ex,frame.ey))<=tol&&std::abs(dot(frame.ex,frame.ez))<=tol&&std::abs(dot(frame.ey,frame.ez))<=tol&&
            norm(sub(cross(frame.ex,frame.ey),frame.ez))<=tol,"Solid frame must be orthonormal and right handed");
    const V3 axes[3]={frame.ex,frame.ey,frame.ez};
    RStructuralPiezoLayeredMeshOperators result;result.family=RStructuralPiezoLayeredFamily::Solid;
    std::map<uint,V3> catalog;std::map<uint,size_t> index;std::vector<V3> local;
    for(const auto&node:nodes){require(node.physicalId!=0&&node.physicalId!=std::numeric_limits<uint>::max(),"Invalid solid node ID");
        for(double x:node.position)require(std::isfinite(x),"Non-finite solid node");require(catalog.emplace(node.physicalId,node.position).second,"Duplicate solid node ID");}
    for(const auto&node:catalog){index[node.first]=local.size();result.nodeIds.push_back(node.first);auto x=sub(node.second,frame.origin);local.push_back({{dot(x,axes[0]),dot(x,axes[1]),dot(x,axes[2])}});for(double v:local.back())require(std::isfinite(v),"Solid local geometry overflow");}
    std::map<uint,const RStructuralPiezoLayeredMember*> ordered;
    for(const auto&m:members){require(m.physicalId!=0&&m.physicalId!=std::numeric_limits<uint>::max()&&ordered.emplace(m.physicalId,&m).second,"Invalid or duplicate solid member ID");}
    std::vector<Element> elements;std::set<size_t> used;std::set<QString> portIds;
    for(const auto&entry:ordered){const auto&m=*entry.second;require(m.nodeIds.size()==4||m.nodeIds.size()==8,"Solid host must be TET4 or affine HEX8");
        require(elements.empty() || m.nodeIds.size()==elements.front().indices.size(),"Mixed TET4/HEX8 transition meshes require interface constraints and are unsupported");
        require(m.layers.size()==1,"Solid host needs exactly one whole-volume material declaration");const auto&source=m.layers.front();source.validate();
        require(source.width==1&&source.coverage==std::array<double,4>{{-1,1,-1,1}},"Solid material occupies its full mesh volume; width/partial coverage unsupported");
        if(!source.material.piezoelectric)require(source.lowerZ==0&&source.upperZ==1&&source.polarity==1,"Passive solid requires neutral electrode planes 0..1 and polarity +1");
        Element e;e.source=&source;std::set<uint> unique;
        for(uint id:m.nodeIds){require(index.count(id)&&unique.insert(id).second,"Invalid solid connectivity");const size_t ni=index.at(id);used.insert(ni);e.indices.push_back(ni);e.x.push_back(local[ni]);}
        e.op=RStructuralPiezoSolidKernel::assemble(e.x,source.material);
        for(const auto&other:elements)require(!overlap(e.x,other.x,tol),"Solid host interiors overlap");
        if(source.material.piezoelectric)portIds.insert(source.portId);
        result.memberIds.push_back(entry.first);result.memberNodeIndices.emplace_back();for(size_t ni:e.indices)result.memberNodeIndices.back().push_back(uint(ni));elements.push_back(std::move(e));
    }
    require(used.size()==nodes.size(),"Solid mesh contains unused nodes");
    auto&op=result.operators;const size_t n=3*nodes.size(),np=portIds.size();op.dofs=uint(n);op.portIds.assign(portIds.begin(),portIds.end());
    op.mass.assign(n*n,0);op.stiffness.assign(n*n,0);op.coupling.assign(n*np,0);op.capacitance.assign(np*np,0);
    // Mechanical and electromechanical assembly is first in common material
    // coordinates; all three translation components are rotated back at end.
    for(const auto&e:elements){const size_t nd=e.op.mechanical.dofs;
        for(size_t i=0;i<nd;++i)for(size_t j=0;j<nd;++j){size_t gi=3*e.indices[i/3]+i%3,gj=3*e.indices[j/3]+j%3;op.mass[gi*n+gj]+=e.op.mechanical.mass[i*nd+j];op.stiffness[gi*n+gj]+=e.op.mechanical.stiffness[i*nd+j];}
        op.totalMass+=e.op.mechanical.totalMass;
    }
    for(size_t p=0;p<np;++p){
        std::set<size_t> electricNodes;std::map<std::vector<size_t>,uint> faceCounts;const RStructuralPiezoLayer*declaration=nullptr;
        for(const auto&e:elements)if(e.source->portId==op.portIds[p]){
            if(!declaration)declaration=e.source;
            require(e.source->lowerZ==declaration->lowerZ&&e.source->upperZ==declaration->upperZ&&e.source->polarity==declaration->polarity,"One solid port requires common electrode planes and polarity");
            electricNodes.insert(e.indices.begin(),e.indices.end());
            for(const auto&f:faces(e.indices.size())){std::vector<size_t> key;for(size_t i:f)key.push_back(e.indices[i]);std::sort(key.begin(),key.end());require(++faceCounts[key]<=2,"Nonmanifold solid electrical face");}
        }
        const auto&d=*declaration;const double h=d.upperZ-d.lowerZ;require(std::isfinite(h)&&h>0,"Unrepresentable solid electrode-plane separation");const double ztol=tol*h;
        std::vector<size_t> en(electricNodes.begin(),electricNodes.end());const size_t ne=en.size();std::map<size_t,size_t> ei;
        for(size_t i=0;i<ne;++i){ei[en[i]]=i;const double z=local[en[i]][2];require(z>=d.lowerZ-ztol&&z<=d.upperZ+ztol,"Piezo solid node lies outside electrode planes");}
        std::vector<size_t> component(ne);for(size_t i=0;i<ne;++i)component[i]=i;
        const auto root=[&](size_t node){while(component[node]!=node)node=component[node];return node;};
        for(const auto&e:elements)if(e.source->portId==op.portIds[p])for(size_t node:e.indices){size_t a=root(ei.at(e.indices[0])),b=root(ei.at(node));component[b]=a;}
        std::map<size_t,std::pair<bool,bool>> electrodes;
        for(size_t i=0;i<ne;++i)electrodes[root(i)];
        for(const auto&face:faceCounts)if(face.second==1){bool l=true,u=true;for(size_t node:face.first){l=l&&std::abs(local[node][2]-d.lowerZ)<=ztol;u=u&&std::abs(local[node][2]-d.upperZ)<=ztol;}auto&flags=electrodes[root(ei.at(face.first[0]))];flags.first=flags.first||l;flags.second=flags.second||u;}
        for(const auto&e:electrodes)require(e.second.first&&e.second.second,"Each connected piezo solid domain requires boundary electrode faces on both planes");
        V g(n*ne),c(ne*ne);
        for(const auto&e:elements)if(e.source->portId==op.portIds[p]){const size_t nn=e.indices.size();for(size_t i=0;i<3*nn;++i)for(size_t j=0;j<nn;++j)g[(3*e.indices[i/3]+i%3)*ne+ei.at(e.indices[j])]+=e.op.coupling[i*nn+j];for(size_t i=0;i<nn;++i)for(size_t j=0;j<nn;++j)c[ei.at(e.indices[i])*ne+ei.at(e.indices[j])]+=e.op.dielectric[i*nn+j];}
        std::vector<size_t> free;V v(ne);
        for(size_t i=0;i<ne;++i){const double z=local[en[i]][2];if(std::abs(z-d.lowerZ)>ztol&&std::abs(z-d.upperZ)>ztol)free.push_back(i);else v[i]=std::abs(z-d.lowerZ)<=ztol?0:double(d.polarity);}
        const size_t nf=free.size();V cf(nf*nf),cv(nf),gf(n*nf);
        for(size_t i=0;i<nf;++i){for(size_t j=0;j<nf;++j)cf[i*nf+j]=c[free[i]*ne+free[j]];for(size_t j=0;j<ne;++j)cv[i]+=c[free[i]*ne+j]*v[j];for(size_t j=0;j<n;++j)gf[j*nf+i]=g[j*ne+free[i]];}
        Factor factor(cf,nf);const auto correction=factor.solve(cv);
        V solvedG(n*nf);for(size_t i=0;i<n;++i){V rhs(nf);for(size_t j=0;j<nf;++j)rhs[j]=gf[i*nf+j];const auto x=factor.solve(rhs);for(size_t j=0;j<nf;++j)solvedG[i*nf+j]=x[j];}
        for(size_t i=0;i<n;++i){for(size_t j=0;j<ne;++j)op.coupling[i*np+p]+=g[i*ne+j]*v[j];for(size_t j=0;j<nf;++j)op.coupling[i*np+p]-=gf[i*nf+j]*correction[j];for(size_t j=0;j<n;++j)for(size_t k=0;k<nf;++k)op.stiffness[i*n+j]+=gf[i*nf+k]*solvedG[j*nf+k];}
        double capacitance=0;for(size_t i=0;i<ne;++i)for(size_t j=0;j<ne;++j)capacitance+=v[i]*c[i*ne+j]*v[j];for(size_t i=0;i<nf;++i)capacitance-=cv[i]*correction[i];op.capacitance[p*np+p]=capacitance;
    }
    // Block congruence: local displacement component a = axis[a] dot u_global.
    for(auto*matrix:{&op.mass,&op.stiffness}){V rotated(n*n);for(size_t i=0;i<nodes.size();++i)for(size_t j=0;j<nodes.size();++j)for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)for(size_t k=0;k<3;++k)for(size_t l=0;l<3;++l)rotated[(3*i+a)*n+3*j+b]+=axes[k][a]*(*matrix)[(3*i+k)*n+3*j+l]*axes[l][b];*matrix=std::move(rotated);}
    symmetrizeRoundoff(op.mass,n,"Solid mass symmetry defect exceeds roundoff tolerance");
    symmetrizeRoundoff(op.stiffness,n,"Solid stiffness symmetry defect exceeds roundoff tolerance");
    V globalG(n*np);for(size_t i=0;i<nodes.size();++i)for(size_t a=0;a<3;++a)for(size_t k=0;k<3;++k)for(size_t p=0;p<np;++p)globalG[(3*i+a)*np+p]+=axes[k][a]*op.coupling[(3*i+k)*np+p];op.coupling=std::move(globalG);
    result.physicalTransform.assign(6*nodes.size()*n,0);for(size_t i=0;i<nodes.size();++i)for(size_t a=0;a<3;++a)result.physicalTransform[(6*i+a)*n+3*i+a]=1;
    op.validate();return result;
}
