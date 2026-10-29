#include "rstructural_piezo_layered_mesh.h"
#include "rstructural_piezo_solid_kernel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace
{
using V3=std::array<double,3>;
using V2=std::array<double,2>;
double dot(const V3 &a,const V3 &b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
V3 difference(const V3 &a,const V3 &b) {return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
double norm(const V3 &a) {return std::hypot(std::hypot(a[0],a[1]),a[2]);}
void finite(const V3 &a) {for(double v:a) if(!std::isfinite(v)) throw std::invalid_argument("Planar mesh requires finite vectors");}
void tolerance(double t) {if(!std::isfinite(t)||t<0||t>1.0e-6) throw std::invalid_argument("Planar relative tolerance must be in [0,1e-6]");}
size_t count(size_t a,size_t b)
{
    if(b && a>std::numeric_limits<size_t>::max()/b) throw std::invalid_argument("Planar matrix count overflow");
    return a*b;
}
double cross(const V2 &a,const V2 &b) {return a[0]*b[1]-a[1]*b[0];}
V2 sub(const V2 &a,const V2 &b) {return {{a[0]-b[0],a[1]-b[1]}};}
//! Separating axes detect positive area overlap, allowing conforming shared edges.
bool plateOverlap(const std::vector<V2> &a,const std::vector<V2> &b)
{
    for(const auto *polygon:{&a,&b}) for(size_t i=0;i<4;++i)
    {
        const auto edge=sub((*polygon)[(i+1)%4],(*polygon)[i]);
        double amin=std::numeric_limits<double>::infinity(),amax=-amin,bmin=amin,bmax=-amin;
        // Use a relative origin to avoid projection cancellation on translated meshes.
        for(const auto &p:a) {const double v=cross(edge,sub(p,a[0])); amin=std::min(amin,v); amax=std::max(amax,v);}
        for(const auto &p:b) {const double v=cross(edge,sub(p,a[0])); bmin=std::min(bmin,v); bmax=std::max(bmax,v);}
        if(std::min(amax,bmax)<=std::max(amin,bmin)) return false;
    }
    return true;
}
bool beamOverlap(const std::vector<V2> &a,const std::vector<V2> &b,double t)
{
    const auto d=sub(a[1],a[0]); const double length=std::hypot(d[0],d[1]);
    const auto b0=sub(b[0],a[0]),b1=sub(b[1],a[0]);
    const double scale=std::max({length,std::hypot(b0[0],b0[1]),std::hypot(b1[0],b1[1])});
    const V2 direction{{d[0]/length,d[1]/length}};
    if(std::abs(cross(direction,b0))>t*scale || std::abs(cross(direction,b1))>t*scale) return false;
    const double s0=(d[0]/length)*b0[0]+(d[1]/length)*b0[1];
    const double s1=(d[0]/length)*b1[0]+(d[1]/length)*b1[1];
    return std::min(length,std::max(s0,s1))>std::max(0.0,std::min(s0,s1));
}
}

RStructuralPiezoLayeredMeshOperators RStructuralPiezoLayeredMesh::assemble(
    RStructuralPiezoLayeredFamily family,const RStructuralPiezoPlanarFrame &frame,
    const std::vector<RStructuralPiezoLayeredNode> &nodes,
    const std::vector<RStructuralPiezoLayeredMember> &members,double tol,double shearCorrection)
{
    if(family==RStructuralPiezoLayeredFamily::Solid)
        return rStructuralPiezoAssembleSolidMesh(frame,nodes,members,tol);
    tolerance(tol);
    if(family!=RStructuralPiezoLayeredFamily::Beam && family!=RStructuralPiezoLayeredFamily::Plate)
        throw std::invalid_argument("Unknown planar mesh family");
    if(nodes.empty()||members.empty()) throw std::invalid_argument("Planar mesh requires nodes and members");
    finite(frame.origin); finite(frame.ex); finite(frame.ey); finite(frame.ez);
    const V3 xy{{frame.ex[1]*frame.ey[2]-frame.ex[2]*frame.ey[1],frame.ex[2]*frame.ey[0]-frame.ex[0]*frame.ey[2],frame.ex[0]*frame.ey[1]-frame.ex[1]*frame.ey[0]}};
    if(std::abs(norm(frame.ex)-1)>tol||std::abs(norm(frame.ey)-1)>tol||std::abs(norm(frame.ez)-1)>tol
       ||std::abs(dot(frame.ex,frame.ey))>tol||std::abs(dot(frame.ex,frame.ez))>tol
       ||std::abs(dot(frame.ey,frame.ez))>tol||norm(difference(xy,frame.ez))>tol)
        throw std::invalid_argument("Planar frame must be orthonormal and right handed");
    const uint nd=family==RStructuralPiezoLayeredFamily::Beam?3:5;
    if(nodes.size()>std::numeric_limits<uint>::max()/nd) throw std::invalid_argument("Planar DOF count overflow");
    RStructuralPiezoLayeredMeshOperators result; result.family=family;
    std::map<uint,V3> catalog;
    for(const auto &node:nodes)
    {
        if(node.physicalId==0 || node.physicalId==std::numeric_limits<uint>::max())
            throw std::invalid_argument("Planar node ID must be nonzero and not the reserved sentinel");
        finite(node.position);
        if(!catalog.emplace(node.physicalId,node.position).second) throw std::invalid_argument("Duplicate planar node ID");
    }
    // The plane's tolerance follows the mesh span, not its distance along the
    // plane from the frame origin. A distant in-plane translation must never
    // turn a real out-of-plane chord into a shorter projected element.
    double meshScale=0.0;
    const auto &anchor=catalog.begin()->second;
    for(const auto &entry:catalog)
        meshScale=std::max(meshScale,norm(difference(entry.second,anchor)));
    if(!std::isfinite(meshScale)) throw std::invalid_argument("Planar mesh span exceeds finite range");
    std::map<uint,size_t> nodeIndex;
    std::vector<V2> projected;
    for(const auto &entry:catalog)
    {
        nodeIndex[entry.first]=result.nodeIds.size(); result.nodeIds.push_back(entry.first);
        const V3 relative=difference(entry.second,frame.origin); finite(relative);
        if(!std::isfinite(norm(relative))) throw std::invalid_argument("Planar geometry scale exceeds finite range");
        const auto &normalAxis=family==RStructuralPiezoLayeredFamily::Beam?frame.ey:frame.ez;
        const double normal=dot(relative,normalAxis);
        double projectionMagnitude=0.0;
        for(size_t axis=0;axis<3;++axis) projectionMagnitude+=std::abs(relative[axis]*normalAxis[axis]);
        const double projectionRoundoff=32.0*std::numeric_limits<double>::epsilon()*projectionMagnitude;
        if(!std::isfinite(normal)||!std::isfinite(projectionRoundoff)||
           std::abs(normal)>tol*meshScale+projectionRoundoff)
            throw std::invalid_argument("Node lies outside the declared mesh plane");
        projected.push_back({{dot(relative,frame.ex),dot(relative,family==RStructuralPiezoLayeredFamily::Beam?frame.ez:frame.ey)}});
        for(double v:projected.back()) if(!std::isfinite(v)) throw std::invalid_argument("Projected geometry overflow");
    }
    std::map<uint,const RStructuralPiezoLayeredMember*> ordered;
    for(const auto &member:members)
    {
        if(member.physicalId==0 || member.physicalId==std::numeric_limits<uint>::max())
            throw std::invalid_argument("Planar member ID must be nonzero and not the reserved sentinel");
        if(!ordered.emplace(member.physicalId,&member).second) throw std::invalid_argument("Duplicate planar member ID");
    }
    struct Element { RStructuralPiezoElementOperators op; std::vector<uint> indices; std::vector<double> transform; std::vector<V2> geometry; };
    std::vector<Element> elements;
    std::set<uint> used; std::set<QString> ports;
    for(const auto &entry:ordered)
    {
        const auto &member=*entry.second; Element element;
        if(member.nodeIds.size()!=(nd==3?2u:4u)) throw std::invalid_argument("Invalid planar member node count");
        std::set<uint> unique;
        for(uint id:member.nodeIds)
        {
            if(!unique.insert(id).second||!nodeIndex.count(id)) throw std::invalid_argument("Invalid planar member connectivity");
            used.insert(id); element.indices.push_back(uint(nodeIndex.at(id))); element.geometry.push_back(projected[nodeIndex.at(id)]);
        }
        if(nd==3)
        {
            const auto d=sub(element.geometry[1],element.geometry[0]); const double l=std::hypot(d[0],d[1]);
            element.op=RStructuralPiezoBeamLayerKernel::assemble(l,member.layers);
            element.transform.assign(36,0.0);
            for(size_t n=0;n<2;++n)
            {
                const size_t i=3*n;
                element.transform[i*6+i]=d[0]/l; element.transform[i*6+i+1]=d[1]/l;
                element.transform[(i+1)*6+i]=-d[1]/l; element.transform[(i+1)*6+i+1]=d[0]/l;
                element.transform[(i+2)*6+i+2]=1;
            }
        }
        else
        {
            std::array<V2,4> xyNodes; std::copy(element.geometry.begin(),element.geometry.end(),xyNodes.begin());
            element.op=RStructuralPiezoPlateMitc4Kernel::assemble(xyNodes,member.layers,shearCorrection);
            element.transform.assign(400,0.0); for(size_t i=0;i<20;++i) element.transform[i*20+i]=1;
        }
        for(const auto &other:elements)
            if(nd==3?beamOverlap(element.geometry,other.geometry,tol):plateOverlap(element.geometry,other.geometry))
                throw std::invalid_argument("Planar host member interiors overlap");
        for(const auto &port:element.op.portIds) ports.insert(port);
        result.memberIds.push_back(entry.first); result.memberNodeIndices.push_back(element.indices);
        elements.push_back(std::move(element));
    }
    if(used.size()!=nodes.size()) throw std::invalid_argument("Planar mesh contains unused nodes");
    auto &op=result.operators; op.dofs=uint(nodes.size()*nd); op.portIds.assign(ports.begin(),ports.end());
    const size_t n=op.dofs,p=ports.size();
    op.mass.assign(count(n,n),0); op.stiffness.assign(count(n,n),0); op.coupling.assign(count(n,p),0); op.capacitance.assign(count(p,p),0);
    for(const auto &e:elements)
    {
        const size_t m=e.op.dofs,ep=e.op.portIds.size(); std::vector<size_t> pm;
        for(const auto &id:e.op.portIds) pm.push_back(size_t(std::lower_bound(op.portIds.begin(),op.portIds.end(),id)-op.portIds.begin()));
        for(size_t i=0;i<m;++i)
        {
            const size_t gi=nd*e.indices[i/nd]+i%nd;
            for(size_t j=0;j<m;++j)
            {
                const size_t gj=nd*e.indices[j/nd]+j%nd;
                // Exact congruence of the full consistent element operators.
                for(size_t a=0;a<m;++a) if(e.transform[a*m+i]!=0)
                    for(size_t b=0;b<m;++b) if(e.transform[b*m+j]!=0)
                    {
                        const double factor=e.transform[a*m+i]*e.transform[b*m+j];
                        op.mass[gi*n+gj]+=factor*e.op.mass[a*m+b];
                        op.stiffness[gi*n+gj]+=factor*e.op.stiffness[a*m+b];
                    }
            }
            for(size_t a=0;a<m;++a) for(size_t k=0;k<ep;++k) op.coupling[gi*p+pm[k]]+=e.transform[a*m+i]*e.op.coupling[a*ep+k];
        }
        for(size_t i=0;i<ep;++i) for(size_t j=0;j<ep;++j) op.capacitance[pm[i]*p+pm[j]]+=e.op.capacitance[i*ep+j];
        op.totalMass+=e.op.totalMass;
    }
    op.validate();
    result.physicalTransform.assign(count(count(nodes.size(),6),n),0);
    for(size_t i=0;i<nodes.size();++i)
    {
        const std::vector<V3> axes=nd==3?std::vector<V3>{frame.ex,frame.ez,frame.ey}:std::vector<V3>{frame.ex,frame.ey,frame.ez,frame.ex,frame.ey};
        for(size_t j=0;j<nd;++j) for(size_t k=0;k<3;++k)
            result.physicalTransform[(6*i+(j<(nd==3?2u:3u)?0:3)+k)*n+nd*i+j]=axes[j][k];
    }
    return result;
}

std::vector<double> RStructuralPiezoLayeredMeshOperators::reconstructPhysical(const std::vector<double> &q) const
{
    const size_t n=operators.dofs,rows=count(nodeIds.size(),6);
    if(q.size()!=n||physicalTransform.size()!=count(rows,n)) throw std::invalid_argument("Invalid planar state or transform count");
    for(double v:q) if(!std::isfinite(v)) throw std::invalid_argument("Planar state must be finite");
    std::vector<double> result(rows,0);
    for(size_t i=0;i<rows;++i) {for(size_t j=0;j<n;++j) result[i]+=physicalTransform[i*n+j]*q[j]; if(!std::isfinite(result[i])) throw std::invalid_argument("Physical state overflow");}
    return result;
}
std::vector<double> RStructuralPiezoLayeredMeshOperators::projectPhysical(const std::vector<double> &x,RStructuralPiezoPhysicalDomain domain,double tol) const
{
    tolerance(tol);
    if(domain!=RStructuralPiezoPhysicalDomain::Displacement&&domain!=RStructuralPiezoPhysicalDomain::Velocity&&domain!=RStructuralPiezoPhysicalDomain::ForceAndMoment) throw std::invalid_argument("Unknown physical projection domain");
    const size_t n=operators.dofs,rows=count(nodeIds.size(),6);
    if(x.size()!=rows||physicalTransform.size()!=count(rows,n)) throw std::invalid_argument("Invalid physical vector count");
    for(double v:x) if(!std::isfinite(v)) throw std::invalid_argument("Physical inputs must be finite");
    std::vector<double> q(n,0);
    for(size_t i=0;i<rows;++i) for(size_t j=0;j<n;++j) q[j]+=physicalTransform[i*n+j]*x[i];
    const auto represented=reconstructPhysical(q);
    for(size_t i=0;i<rows;i+=3)
    {
        const double scale=std::max({std::abs(x[i]),std::abs(x[i+1]),std::abs(x[i+2])});
        V3 block{},residual{};
        for(size_t k=0;k<3;++k)
        {
            block[k]=scale==0?0:x[i+k]/scale;
            residual[k]=scale==0?-represented[i+k]:block[k]-represented[i+k]/scale;
        }
        if(!std::isfinite(norm(residual))||norm(residual)>tol*norm(block)) throw std::invalid_argument(domain==RStructuralPiezoPhysicalDomain::ForceAndMoment?"Unsupported planar force or moment component":"Unsupported planar displacement or velocity component");
    }
    return q;
}
