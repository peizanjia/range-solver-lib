#include <QtTest>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include "rstructural_piezo_electrical_system.h"
#include "rstructural_piezo_layered_mesh.h"

namespace
{
using Vector=std::vector<double>;
using Node=RStructuralPiezoLayeredNode;
using Member=RStructuralPiezoLayeredMember;
using Family=RStructuralPiezoLayeredFamily;
using Domain=RStructuralPiezoPhysicalDomain;
RStructuralPiezoLayer layer(bool active=false,double lower=.1,double upper=.3)
{
    RStructuralPiezoLayer l; l.lowerZ=lower; l.upperZ=upper;
    l.material=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
    if(active) {l.material.piezoelectric=true; l.material.coupling[12]=2; l.material.coupling[13]=2;
        l.material.permittivity={{4,0,0,0,4,0,0,0,4}}; l.portId=QStringLiteral("continuous");}
    return l;
}
bool close(double a,double b) {return std::abs(a-b)<1e-10*(1+std::max(std::abs(a),std::abs(b)));}
Vector product(const Vector &a,const Vector &q)
{
    Vector r(q.size(),0); for(size_t i=0;i<q.size();++i) for(size_t j=0;j<q.size();++j) r[i]+=a[i*q.size()+j]*q[j]; return r;
}
double work(const Vector &a,const Vector &b) {double r=0; for(size_t i=0;i<a.size();++i) r+=a[i]*b[i]; return r;}
std::vector<RStructuralPiezoLayer> cantileverRefinementLayers()
{
    auto host=layer(); host.material=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
    host.lowerZ=-.1; host.upperZ=.1; host.width=.25;
    auto patch=layer(true,.1,.15); patch.material=RStructuralPiezoContinuumMaterial::isotropic(20,.25,7.5);
    patch.material.piezoelectric=true; patch.material.coupling[12]=2; patch.material.coupling[13]=2;
    patch.material.permittivity={{4,0,0,0,4,0,0,0,4}}; patch.width=.18;
    patch.portId=QStringLiteral("shared-full-length-electrode");
    return {host,patch};
}
Vector solveDense(Vector matrix,Vector rhs)
{
    const size_t n=rhs.size();
    if(matrix.size()!=n*n) throw std::invalid_argument("Invalid dense system dimensions");
    for(size_t k=0;k<n;++k)
    {
        size_t pivot=k; double largest=std::abs(matrix[k*n+k]);
        for(size_t i=k+1;i<n;++i) if(std::abs(matrix[i*n+k])>largest)
        { pivot=i; largest=std::abs(matrix[i*n+k]); }
        if(!(largest>0) || !std::isfinite(largest)) throw std::runtime_error("Singular dense refinement stiffness");
        if(pivot!=k)
        {
            for(size_t j=k;j<n;++j) std::swap(matrix[k*n+j],matrix[pivot*n+j]);
            std::swap(rhs[k],rhs[pivot]);
        }
        for(size_t i=k+1;i<n;++i)
        {
            const double factor=matrix[i*n+k]/matrix[k*n+k]; matrix[i*n+k]=0;
            for(size_t j=k+1;j<n;++j) matrix[i*n+j]-=factor*matrix[k*n+j];
            rhs[i]-=factor*rhs[k];
        }
    }
    Vector solution(n,0);
    for(size_t i=n;i-- > 0;)
    {
        double value=rhs[i]; for(size_t j=i+1;j<n;++j) value-=matrix[i*n+j]*solution[j];
        solution[i]=value/matrix[i*n+i];
        if(!std::isfinite(solution[i])) throw std::runtime_error("Non-finite dense refinement solution");
    }
    return solution;
}
class MindlinNavierSeries
{
public:
    MindlinNavierSeries(double length,double width,double thickness,double young,double poisson,
                        double pressure,uint maximumOddMode)
    {
        const double pi=std::acos(-1.0),g=young/(2.0*(1.0+poisson));
        const double d11=young/(1.0-poisson*poisson)*thickness*thickness*thickness/12.0;
        const double d12=poisson*young/(1.0-poisson*poisson)*thickness*thickness*thickness/12.0;
        const double d66=g*thickness*thickness*thickness/12.0;
        const double shear=(5.0/6.0)*g*thickness;
        for(uint m=1;m<=maximumOddMode;m+=2) for(uint n=1;n<=maximumOddMode;n+=2)
        {
            const double alpha=double(m)*pi/length,beta=double(n)*pi/width;
            // Production kinematics are gamma_x=w_x+RY, gamma_y=w_y-RX;
            // these modal amplitudes use w=W sin(ax)sin(by), RX=A sin(ax)cos(by),
            // and RY=B cos(ax)sin(by). The bending cross term is therefore
            // -(D12+D66)*alpha*beta in the RX/RY block.
            const Vector modal{
                shear*(alpha*alpha+beta*beta),-shear*beta,shear*alpha,
                -shear*beta,d11*beta*beta+d66*alpha*alpha+shear,-(d12+d66)*alpha*beta,
                shear*alpha,-(d12+d66)*alpha*beta,d11*alpha*alpha+d66*beta*beta+shear};
            const double pressureMode=16.0*pressure/(pi*pi*double(m)*double(n));
            const auto response=solveDense(modal,{pressureMode,0.0,0.0});
            modes.push_back({alpha,beta,response[0]});
        }
    }

    double displacement(double x,double y) const
    {
        double value=0.0;
        for(const auto &mode:modes)
            value+=mode.transverse*std::sin(mode.alpha*x)*std::sin(mode.beta*y);
        return value;
    }

private:
    struct Mode { double alpha,beta,transverse; };
    std::vector<Mode> modes;
};
std::vector<Node> beamNodes() {return {{30,{{2,0,1}}},{10,{{0,0,0}}},{20,{{1,0,0}}}};}
std::vector<Member> beamMembers() {return {{9,{20,30},{layer(),layer(true,.3,.4)}},{4,{10,20},{layer()}}};}
RStructuralPiezoPlanarFrame rotatedFrame()
{
    // Exact cyclic spatial rotation followed by an oblique rotation about ex.
    const double c=std::cos(.63),s=std::sin(.63);
    RStructuralPiezoPlanarFrame f; f.origin={{3,-4,7}}; f.ex={{0,1,0}}; f.ey={{s,0,c}}; f.ez={{c,0,-s}}; return f;
}
std::vector<Node> rotated(std::vector<Node> nodes,const RStructuralPiezoPlanarFrame &f)
{
    for(auto &node:nodes) {const auto p=node.position; for(size_t k=0;k<3;++k) node.position[k]=f.origin[k]+p[0]*f.ex[k]+p[1]*f.ey[k]+p[2]*f.ez[k];} return nodes;
}
}
class TestStructuralPiezoLayeredMesh:public QObject
{
    Q_OBJECT
private slots:
    void beamProductionAssemblyMassPortsAndOrdering();
    void multiElementBeamUniformLoadConvergesToCompositeClosedForm();
    void spatialRotationAndPhysicalAction();
    void unsupportedLoadsAndInitialState();
    void continuousPlateElectrodeAndDistortedPatch();
    void plateMeshRefinementAgainstMindlinContinuumAndLayeredElectrodeIntegral();
    void malformedInputs();
};
void TestStructuralPiezoLayeredMesh::beamProductionAssemblyMassPortsAndOrdering()
{
    const auto mesh=RStructuralPiezoLayeredMesh::assemble(Family::Beam,{},beamNodes(),beamMembers());
    QCOMPARE(mesh.nodeIds,(std::vector<uint>{10,20,30})); QCOMPARE(mesh.memberIds,(std::vector<uint>{4,9}));
    const auto a=RStructuralPiezoBeamLayerKernel::assemble(1,{layer()});
    const auto b=RStructuralPiezoBeamLayerKernel::assemble(std::sqrt(2.),{layer(),layer(true,.3,.4)});
    QVERIFY(close(mesh.operators.totalMass,7*(.2+std::sqrt(2.)*.3)));
    QVERIFY(close(mesh.operators.capacitance[0],b.capacitance[0]));
    // Independently apply chord kinematics to a common nodal test state.
    Vector q{.4,-.2,.1,.3,.6,-.4,-.7,.2,.5};
    Vector qa(q.begin(),q.begin()+6),qb(6);
    const double c=1/std::sqrt(2.);
    for(size_t i=0;i<2;++i) {qb[3*i]=c*(q[3+3*i]+q[4+3*i]); qb[3*i+1]=c*(-q[3+3*i]+q[4+3*i]); qb[3*i+2]=q[5+3*i];}
    QVERIFY(close(work(q,product(mesh.operators.mass,q)),work(qa,product(a.mass,qa))+work(qb,product(b.mass,qb))));
    QVERIFY(close(work(q,product(mesh.operators.stiffness,q)),work(qa,product(a.stiffness,qa))+work(qb,product(b.stiffness,qb))));
    QVERIFY(close(work(q,mesh.operators.coupling),work(qb,b.coupling)));
    // Eccentric translation/rotation inertia survives the shared node assembly.
    QVERIFY(std::abs(mesh.operators.mass[2])>1e-5);
    auto nodes=beamNodes(); auto members=beamMembers(); std::reverse(nodes.begin(),nodes.end()); std::reverse(members.begin(),members.end());
    const auto reordered=RStructuralPiezoLayeredMesh::assemble(Family::Beam,{},nodes,members);
    QCOMPARE(reordered.operators.mass,mesh.operators.mass); QCOMPARE(reordered.operators.stiffness,mesh.operators.stiffness);
}
void TestStructuralPiezoLayeredMesh::multiElementBeamUniformLoadConvergesToCompositeClosedForm()
{
    constexpr double length=1.2,distributedLoad=.1,sampleX=.37*length;
    const auto layers=cantileverRefinementLayers();
    // This section reference is calculated directly from E*width z-moments,
    // independently of the production beam-layer reduction and assembly.
    const std::array<double,2> young{{30,20}};
    double A=0,B=0,D=0;
    for(size_t i=0;i<layers.size();++i)
    {
        const auto &l=layers[i]; const double lo=l.lowerZ,hi=l.upperZ,ew=young[i]*l.width;
        A+=ew*(hi-lo); B+=ew*(hi*hi-lo*lo)/2; D+=ew*(hi*hi*hi-lo*lo*lo)/3;
    }
    const double effectiveEI=D-B*B/A;
    QVERIFY(A>0); QVERIFY(B!=0); QVERIFY(effectiveEI>0);
    const double exact=distributedLoad*sampleX*sampleX*
        (6*length*length-4*length*sampleX+sampleX*sampleX)/(24*effectiveEI);
    QVERIFY(std::isfinite(exact)); QVERIFY(exact>0);

    const std::array<uint,3> segmentCounts{{2,4,8}};
    std::array<double,3> errors{};
    for(size_t refinement=0;refinement<segmentCounts.size();++refinement)
    {
        const uint segments=segmentCounts[refinement]; const double h=length/segments;
        std::vector<Node> nodes; std::vector<Member> members;
        for(uint i=0;i<=segments;++i) nodes.push_back({i+1,{{i*h,0,0}}});
        for(uint i=0;i<segments;++i) members.push_back({i+1,{i+1,i+2},layers});
        const auto mesh=RStructuralPiezoLayeredMesh::assemble(Family::Beam,{},nodes,members);
        const size_t n=mesh.operators.dofs;
        QCOMPARE(n,size_t(3)*(segments+1)); QCOMPARE(mesh.operators.portIds.size(),size_t(1));
        QCOMPARE(mesh.operators.portIds[0],QStringLiteral("shared-full-length-electrode"));

        // Re-sum each production element's M/K/G/Ce into the shared mesh DOFs.
        // This checks the repeated BEAM1 assembly and the single continuous port.
        Vector summedMass(n*n,0),summedStiffness(n*n,0),summedCoupling(n,0);
        double summedCapacitance=0,summedTotalMass=0;
        for(size_t e=0;e<members.size();++e)
        {
            const auto element=RStructuralPiezoBeamLayerKernel::assemble(h,layers);
            const auto &connectivity=mesh.memberNodeIndices[e];
            for(size_t i=0;i<6;++i)
            {
                const size_t gi=3*connectivity[i/3]+i%3;
                summedCoupling[gi]+=element.coupling[i];
                for(size_t j=0;j<6;++j)
                {
                    const size_t gj=3*connectivity[j/3]+j%3;
                    summedMass[gi*n+gj]+=element.mass[i*6+j];
                    summedStiffness[gi*n+gj]+=element.stiffness[i*6+j];
                }
            }
            summedCapacitance+=element.capacitance[0]; summedTotalMass+=element.totalMass;
        }
        for(size_t i=0;i<n*n;++i)
        { QVERIFY(close(mesh.operators.mass[i],summedMass[i])); QVERIFY(close(mesh.operators.stiffness[i],summedStiffness[i])); }
        for(size_t i=0;i<n;++i) QVERIFY(close(mesh.operators.coupling[i],summedCoupling[i]));
        QVERIFY(close(mesh.operators.capacitance[0],summedCapacitance));
        QVERIFY(close(mesh.operators.totalMass,summedTotalMass));
        double couplingNorm=0;
        for(double value:mesh.operators.coupling) couplingNorm=std::hypot(couplingNorm,value);
        QVERIFY(couplingNorm>0);
        QVERIFY(mesh.operators.capacitance[0]>0);
        QVERIFY(summedCapacitance>0);

        // A short-circuit charge-sensor boundary fixes this shared electrode at
        // V=0; its condensed stiffness must be the assembled c^E stiffness.
        RStructuralPiezoElectricalBoundary shortCircuit;
        shortCircuit.portId=mesh.operators.portIds[0];
        shortCircuit.role=RStructuralPiezoElectricalRole::ShortCircuitChargeSensor;
        RStructuralPiezoElectricalSystem electrical(mesh.operators,{shortCircuit},Vector(n,0),{});
        const auto &shorted=electrical.condensed();
        QCOMPARE(shorted.stiffness,mesh.operators.stiffness);

        // Integrate q against the exact cubic Hermite shape functions on every
        // element. thetaY=-dw/dx, so the consistent end moments have these signs.
        Vector force(n,0);
        for(size_t e=0;e<members.size();++e)
        {
            const size_t left=mesh.memberNodeIndices[e][0],right=mesh.memberNodeIndices[e][1];
            force[3*left+1]+=distributedLoad*h/2;
            force[3*left+2]-=distributedLoad*h*h/12;
            force[3*right+1]+=distributedLoad*h/2;
            force[3*right+2]+=distributedLoad*h*h/12;
        }
        // Clamp u, w and theta at x=0, then solve the actual assembled K_ff.
        const size_t freeCount=n-3; Vector reducedK(freeCount*freeCount),reducedForce(freeCount);
        for(size_t i=0;i<freeCount;++i)
        {
            reducedForce[i]=force[i+3];
            for(size_t j=0;j<freeCount;++j)
                reducedK[i*freeCount+j]=mesh.operators.stiffness[(i+3)*n+j+3];
        }
        const Vector freeSolution=solveDense(std::move(reducedK),std::move(reducedForce));
        Vector displacement(n,0); std::copy(freeSolution.begin(),freeSolution.end(),displacement.begin()+3);
        const auto shortedResponse=electrical.evaluate(displacement,{});
        QCOMPARE(shortedResponse.ports.size(),size_t(1));
        QVERIFY(close(shortedResponse.ports[0].voltage,0));
        double chargeScale=0;
        for(size_t i=0;i<n;++i) chargeScale+=std::abs(mesh.operators.coupling[i]*displacement[i]);
        QVERIFY(chargeScale>0);
        QVERIFY(std::abs(shortedResponse.ports[0].charge)>1e-10*chargeScale);
        const Vector residual=product(mesh.operators.stiffness,displacement);
        double residualNorm=0,loadNorm=0;
        for(size_t i=3;i<n;++i) {residualNorm=std::max(residualNorm,std::abs(residual[i]-force[i])); loadNorm=std::max(loadNorm,std::abs(force[i]));}
        QVERIFY(residualNorm<1e-10*std::max(loadNorm,1e-30));

        const size_t element=std::min<size_t>(segments-1,size_t(std::floor(sampleX/h)));
        const double s=(sampleX-element*h)/h,s2=s*s,s3=s2*s;
        const size_t left=3*element,right=3*(element+1);
        const double interpolated=(1-3*s2+2*s3)*displacement[left+1]
            -h*(s-2*s2+s3)*displacement[left+2]
            +(3*s2-2*s3)*displacement[right+1]
            +h*(s2-s3)*displacement[right+2];
        errors[refinement]=std::abs(interpolated-exact);
        QVERIFY(std::isfinite(interpolated));
        qInfo().noquote() << QStringLiteral("PC-02 A03 n=%1 x/L=0.37 relativeError=%2 chargeRatio=%3 Ce=%4 GNorm=%5")
                .arg(segments)
                .arg(QString::number(errors[refinement]/std::abs(exact),'g',17))
                .arg(QString::number(std::abs(shortedResponse.ports[0].charge)/chargeScale,'g',17))
                .arg(QString::number(mesh.operators.capacitance[0],'g',17))
                .arg(QString::number(couplingNorm,'g',17));
    }
    const double errorFloor=1e-11*std::abs(exact);
    for(double error:errors) QVERIFY(error>errorFloor);
    QVERIFY(errors[1]<.7*errors[0]); QVERIFY(errors[2]<.7*errors[1]);
    // Explicit final qualification tolerance: 0.25% of the independent EB value.
    QVERIFY(errors[2]<2.5e-3*std::abs(exact));
}
void TestStructuralPiezoLayeredMesh::spatialRotationAndPhysicalAction()
{
    for(Family family:{Family::Beam,Family::Plate})
    {
        std::vector<Node> nodes=family==Family::Beam?beamNodes():std::vector<Node>{{1,{{0,0,0}}},{2,{{2,0,0}}},{3,{{2.3,1.8,0}}},{4,{{-.1,2,0}}}};
        std::vector<Member> members=family==Family::Beam?beamMembers():std::vector<Member>{{1,{1,2,3,4},{layer(true)}}};
        const auto a=RStructuralPiezoLayeredMesh::assemble(family,{},nodes,members);
        const auto f=rotatedFrame(); const auto b=RStructuralPiezoLayeredMesh::assemble(family,f,rotated(nodes,f),members);
        for(const auto &pair:{std::make_pair(a.operators.mass,b.operators.mass),std::make_pair(a.operators.stiffness,b.operators.stiffness),std::make_pair(a.operators.coupling,b.operators.coupling)})
            for(size_t i=0;i<pair.first.size();++i) QVERIFY(close(pair.first[i],pair.second[i]));
        Vector q(a.operators.dofs),force(q.size()); for(size_t i=0;i<q.size();++i) {q[i]=.01*(double(i)-4); force[i]=.07*(double(i)+2);}
        const auto global=b.reconstructPhysical(q),load=b.reconstructPhysical(force);
        for(Domain domain:{Domain::Displacement,Domain::Velocity,Domain::ForceAndMoment})
        {
            const auto back=b.projectPhysical(global,domain); for(size_t i=0;i<q.size();++i) QVERIFY(close(back[i],q[i]));
        }
        const auto generalizedLoad=b.projectPhysical(load,Domain::ForceAndMoment);
        QVERIFY(close(work(global,load),work(q,generalizedLoad)));
        // P M P^T physical action, projected back, is the full consistent M action.
        const auto acceleration=b.reconstructPhysical(product(b.operators.mass,q));
        const auto mapped=b.projectPhysical(acceleration,Domain::ForceAndMoment);
        const auto expected=product(a.operators.mass,q); for(size_t i=0;i<q.size();++i) QVERIFY(close(mapped[i],expected[i]));
    }
}
void TestStructuralPiezoLayeredMesh::unsupportedLoadsAndInitialState()
{
    const auto b=RStructuralPiezoLayeredMesh::assemble(Family::Beam,{},beamNodes(),beamMembers());
    for(size_t component:{1u,3u,5u}) for(Domain domain:{Domain::Displacement,Domain::Velocity,Domain::ForceAndMoment})
    {
        Vector x(18,0); x[component]=1e-200;
        QVERIFY_EXCEPTION_THROWN(b.projectPhysical(x,domain),std::invalid_argument);
    }
    const auto p=RStructuralPiezoLayeredMesh::assemble(Family::Plate,{},{{1,{{0,0,0}}},{2,{{1,0,0}}},{3,{{1,1,0}}},{4,{{0,1,0}}}},{{1,{1,2,3,4},{layer()}}});
    Vector x(24,0); x[5]=1;
    QVERIFY_EXCEPTION_THROWN(p.projectPhysical(x,Domain::Displacement),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(p.projectPhysical(x,Domain::ForceAndMoment),std::invalid_argument);
    x[5]=0; x[0]=1; x[5]=1e-12; // Force magnitude must not hide unsupported moment (different units).
    QVERIFY_EXCEPTION_THROWN(p.projectPhysical(x,Domain::ForceAndMoment),std::invalid_argument);
}
void TestStructuralPiezoLayeredMesh::continuousPlateElectrodeAndDistortedPatch()
{
    std::vector<Node> nodes; for(uint y=0;y<3;++y) for(uint x=0;x<3;++x)
        nodes.push_back({3*y+x+1,{{double(x)+(x==1&&y==1?.17:0),double(y)+(x==1&&y==1?-.13:0),0}}});
    std::vector<Member> members;
    for(uint y=0;y<2;++y) for(uint x=0;x<2;++x) members.push_back({2*y+x+1,{3*y+x+1,3*y+x+2,3*(y+1)+x+2,3*(y+1)+x+1},{layer(true)}});
    const auto mesh=RStructuralPiezoLayeredMesh::assemble(Family::Plate,{},nodes,members);
    std::map<uint,size_t> nodeLookup,equationLookup;
    for(size_t i=0;i<nodes.size();++i) nodeLookup.emplace(nodes[i].physicalId,i);
    for(size_t i=0;i<mesh.nodeIds.size();++i) equationLookup.emplace(mesh.nodeIds[i],i);
    QCOMPARE(mesh.operators.portIds.size(),size_t(1));
    Vector summedG(45,0); double summedCe=0;
    for(const auto &member:members)
    {
        std::array<std::array<double,2>,4> xy;
        for(size_t i=0;i<4;++i)
        {
            const auto &position=nodes[nodeLookup.at(member.nodeIds[i])].position;
            xy[i]={{position[0],position[1]}};
        }
        const auto e=RStructuralPiezoPlateMitc4Kernel::assemble(xy,member.layers); summedCe+=e.capacitance[0];
        for(size_t i=0;i<20;++i) summedG[5*equationLookup.at(member.nodeIds[i/5])+i%5]+=e.coupling[i];
    }
    QVERIFY(close(mesh.operators.capacitance[0],summedCe));
    for(size_t i=0;i<45;++i) QVERIFY(close(mesh.operators.coupling[i],summedG[i]));
    QVERIFY(close(mesh.operators.totalMass,4*.2*7));
    for(bool bending:{false,true})
    {
        Vector q(45,0);
        for(size_t i=0;i<9;++i)
        {
            const auto &position=nodes[nodeLookup.at(mesh.nodeIds[i])].position;
            const double x=position[0],y=position[1];
            if(!bending) {q[5*i]=.07*x+.02*y; q[5*i+1]=-.03*y+.02*x;}
            else {q[5*i+2]=.5*(.13*x*x-.11*y*y)+.08*x*y; q[5*i+3]=-.11*y+.08*x; q[5*i+4]=-.13*x-.08*y;}
        }
        const auto force=product(mesh.operators.stiffness,q);
        // Constant membrane stress / bending moment has no interior body load.
        // The quadratic Kirchhoff nodal field has zero MITC tied shear; boundary
        // equivalent tractions and moments remain and must not be zeroed.
        for(size_t k=0;k<5;++k) QVERIFY(std::abs(force[20+k])<2e-11);
        double boundaryNorm=0; for(size_t i=0;i<45;++i) if(i/5!=4) boundaryNorm+=force[i]*force[i];
        QVERIFY(boundaryNorm>1e-6);
    }
}
void TestStructuralPiezoLayeredMesh::plateMeshRefinementAgainstMindlinContinuumAndLayeredElectrodeIntegral()
{
    // The high-order Navier series is an independent Reissner-Mindlin continuum
    // solution for a simply supported square under uniform pressure. A symmetric
    // equal-stiffness laminate keeps the reference homogeneous while retaining a
    // continuous active layer in the assembled electrode operators.
    constexpr double length=1.0,width=1.0,totalThickness=.2,activeThickness=.1;
    constexpr double young=30e9,poisson=.2,pressure=1e4,membraneStrain=1e-5;
    // Compare the independent series truncations at the plate center, then use
    // the same fixed interior sampling grid for every FE refinement. Its 32x32
    // subcell centers never coincide with nodes or element boundaries of the
    // 2x2, 4x4, 8x8, or 16x16 meshes.
    constexpr uint sampleGridSide=32;
    const double pi=std::acos(-1.0),shearModulus=young/(2.0*(1.0+poisson));
    const MindlinNavierSeries coarseSeries(length,width,totalThickness,young,poisson,pressure,301);
    const MindlinNavierSeries fineSeries(length,width,totalThickness,young,poisson,pressure,501);
    const double coarseCenter=coarseSeries.displacement(.5*length,.5*width);
    const double fineCenter=fineSeries.displacement(.5*length,.5*width);
    QVERIFY(std::isfinite(coarseCenter));
    QVERIFY(std::isfinite(fineCenter));
    QVERIFY(fineCenter>0.0);
    const double seriesRelativeDelta=std::abs(fineCenter-coarseCenter)/std::abs(fineCenter);
    QVERIFY(std::isfinite(seriesRelativeDelta));
    QVERIFY(seriesRelativeDelta>0.0);
    qInfo().noquote()<<QStringLiteral("Mindlin Navier reference modes=301:%1 modes=501:%2 relativeDelta:%3")
        .arg(QString::number(coarseCenter,'g',17)).arg(QString::number(fineCenter,'g',17))
        .arg(QString::number(seriesRelativeDelta,'g',17));
    QVERIFY(seriesRelativeDelta<1e-8);
    std::array<double,sampleGridSide*sampleGridSide> referenceSamples{};
    for(uint y=0;y<sampleGridSide;++y) for(uint x=0;x<sampleGridSide;++x)
    {
        const double sampleX=(double(x)+.5)*length/double(sampleGridSide);
        const double sampleY=(double(y)+.5)*width/double(sampleGridSide);
        const double value=fineSeries.displacement(sampleX,sampleY);
        QVERIFY(std::isfinite(value));
        referenceSamples[y*sampleGridSide+x]=value;
    }

    const auto hostMaterial=RStructuralPiezoContinuumMaterial::isotropic(young,poisson,7500.0);
    auto piezoMaterial=hostMaterial;piezoMaterial.piezoelectric=true;
    piezoMaterial.coupling[12]=2.0;piezoMaterial.coupling[13]=2.0;
    piezoMaterial.permittivity={{4e-8,0,0,0,4e-8,0,0,0,4e-8}};
    RStructuralPiezoLayer lower;lower.material=hostMaterial;lower.lowerZ=-.1;lower.upperZ=-.05;
    RStructuralPiezoLayer active;active.material=piezoMaterial;active.lowerZ=-.05;active.upperZ=.05;
    active.portId=QStringLiteral("continuous");
    RStructuralPiezoLayer upper;upper.material=hostMaterial;upper.lowerZ=.05;upper.upperZ=.1;
    const std::vector<RStructuralPiezoLayer> laminate{lower,active,upper};
    const double expectedCapacitance=piezoMaterial.permittivity[8]*length*width/activeThickness;
    const double expectedAffineCharge=piezoMaterial.coupling[12]*membraneStrain*length*width;
    QVERIFY(std::isfinite(expectedCapacitance));
    QVERIFY(std::isfinite(expectedAffineCharge));
    QVERIFY(expectedCapacitance>0.0);QVERIFY(expectedAffineCharge>0.0);

    RStructuralPiezoElectricalBoundary shortCircuit;
    shortCircuit.portId=QStringLiteral("continuous");
    shortCircuit.role=RStructuralPiezoElectricalRole::ShortCircuitChargeSensor;
    const std::array<uint,4> divisions{{2,4,8,16}};
    std::array<double,4> relativeErrors{},affineCharges{};
    for(size_t refinement=0;refinement<divisions.size();++refinement)
    {
        const uint segments=divisions[refinement];
        const double dx=length/segments,dy=width/segments,cellArea=dx*dy;
        std::vector<Node> nodes;
        const auto physicalId=[&](uint x,uint y){return 1+y*(segments+1)+x;};
        for(uint y=0;y<=segments;++y) for(uint x=0;x<=segments;++x)
            nodes.push_back({physicalId(x,y),{{double(x)*dx,double(y)*dy,0.0}}});
        std::map<uint,size_t> nodeIndex;
        for(size_t i=0;i<nodes.size();++i)nodeIndex.emplace(nodes[i].physicalId,i);
        std::vector<Member> members;
        Vector force(5*nodes.size(),0.0);
        for(uint y=0;y<segments;++y) for(uint x=0;x<segments;++x)
        {
            const std::vector<uint> connectivity{
                physicalId(x,y),physicalId(x+1,y),physicalId(x+1,y+1),physicalId(x,y+1)};
            members.push_back({1+y*segments+x+1,connectivity,laminate});
            const double nodalPressure=pressure*cellArea/4.0;
            for(uint id:connectivity)force[5*nodeIndex.at(id)+2]+=nodalPressure;
        }
        const auto mesh=RStructuralPiezoLayeredMesh::assemble(Family::Plate,{},nodes,members);
        QCOMPARE(mesh.operators.dofs,uint(5*nodes.size()));
        QCOMPARE(mesh.operators.portIds,(std::vector<QString>{QStringLiteral("continuous")}));
        QVERIFY(std::isfinite(mesh.operators.capacitance[0]));
        const double capacitanceRelativeError=std::abs(mesh.operators.capacitance[0]-expectedCapacitance)/expectedCapacitance;
        QVERIFY(std::isfinite(capacitanceRelativeError));
        QVERIFY(capacitanceRelativeError<1e-10);

        // Match the Navier basis boundary conditions in the production DOFs:
        // w=0 on all edges, RX=0 on x edges, and RY=0 on y edges. The remaining
        // edge rotations are natural. In-plane motion is fixed to isolate the
        // symmetric plate response; short-circuit condensation leaves K unchanged.
        RStructuralPiezoElectricalSystem electrical(mesh.operators,{shortCircuit},Vector(mesh.operators.dofs,0.0),{});
        const auto &stiffness=electrical.condensed().stiffness;
        const size_t fullDofs=mesh.operators.dofs;
        std::vector<size_t> freeDofs;
        for(size_t ni=0;ni<mesh.nodeIds.size();++ni)
        {
            const auto &position=nodes[nodeIndex.at(mesh.nodeIds[ni])].position;
            const uint x=uint(std::llround(position[0]/dx)),y=uint(std::llround(position[1]/dy));
            if(x>0&&x<segments&&y>0&&y<segments)freeDofs.push_back(5*ni+2);
            if(x>0&&x<segments)freeDofs.push_back(5*ni+3);
            if(y>0&&y<segments)freeDofs.push_back(5*ni+4);
        }
        const size_t freeCount=freeDofs.size();
        Vector reducedStiffness(freeCount*freeCount),reducedForce(freeCount);
        for(size_t i=0;i<freeCount;++i)
        {
            reducedForce[i]=force[freeDofs[i]];
            for(size_t j=0;j<freeCount;++j)
                reducedStiffness[i*freeCount+j]=stiffness[freeDofs[i]*fullDofs+freeDofs[j]];
        }
        const Vector freeSolution=solveDense(std::move(reducedStiffness),std::move(reducedForce));
        Vector q(fullDofs,0.0);
        for(size_t i=0;i<freeCount;++i)q[freeDofs[i]]=freeSolution[i];

        double squaredError=0.0,squaredReference=0.0;
        for(uint sampleYIndex=0;sampleYIndex<sampleGridSide;++sampleYIndex)
            for(uint sampleXIndex=0;sampleXIndex<sampleGridSide;++sampleXIndex)
            {
                const double sampleX=(double(sampleXIndex)+.5)*length/double(sampleGridSide);
                const double sampleY=(double(sampleYIndex)+.5)*width/double(sampleGridSide);
                const uint sampleCellX=uint(std::floor(sampleX/dx));
                const uint sampleCellY=uint(std::floor(sampleY/dy));
                const double localX=(sampleX-double(sampleCellX)*dx)/dx;
                const double localY=(sampleY-double(sampleCellY)*dy)/dy;
                QVERIFY(localX>0.0&&localX<1.0&&localY>0.0&&localY<1.0);
                const std::array<uint,4> sampleNodes{{physicalId(sampleCellX,sampleCellY),
                    physicalId(sampleCellX+1,sampleCellY),physicalId(sampleCellX+1,sampleCellY+1),
                    physicalId(sampleCellX,sampleCellY+1)}};
                const std::array<double,4> shape{{(1-localX)*(1-localY),localX*(1-localY),
                    localX*localY,(1-localX)*localY}};
                double finiteElementValue=0.0;
                for(size_t a=0;a<4;++a)
                    finiteElementValue+=shape[a]*q[5*nodeIndex.at(sampleNodes[a])+2];
                const double continuumValue=referenceSamples[sampleYIndex*sampleGridSide+sampleXIndex];
                QVERIFY(std::isfinite(finiteElementValue));
                const double difference=finiteElementValue-continuumValue;
                QVERIFY(std::isfinite(difference));
                squaredError+=difference*difference;
                squaredReference+=continuumValue*continuumValue;
            }
        QVERIFY(std::isfinite(squaredError));
        QVERIFY(std::isfinite(squaredReference));
        QVERIFY(squaredReference>0.0);
        relativeErrors[refinement]=std::sqrt(squaredError/squaredReference);
        QVERIFY(std::isfinite(relativeErrors[refinement]));
        QVERIFY(relativeErrors[refinement]>0.0);

        Vector affine(fullDofs,0.0);
        for(size_t ni=0;ni<mesh.nodeIds.size();++ni)
            affine[5*ni]=membraneStrain*nodes[nodeIndex.at(mesh.nodeIds[ni])].position[0];
        double generalizedCharge=0.0,couplingScale=0.0;
        for(size_t i=0;i<fullDofs;++i)
        {
            const double term=mesh.operators.coupling[i]*affine[i];
            QVERIFY(std::isfinite(term));
            generalizedCharge+=term;couplingScale+=std::abs(term);
        }
        QVERIFY(std::isfinite(generalizedCharge));
        QVERIFY(std::isfinite(couplingScale));
        QVERIFY(couplingScale>0.0);
        const double generalizedChargeRelativeError=
            std::abs(generalizedCharge-expectedAffineCharge)/std::abs(expectedAffineCharge);
        QVERIFY(std::isfinite(generalizedChargeRelativeError));
        QVERIFY(generalizedChargeRelativeError<1e-10);
        affineCharges[refinement]=generalizedCharge;
        const auto electricalResponse=electrical.evaluate(affine,{});
        QCOMPARE(electricalResponse.ports.size(),size_t(1));
        QCOMPARE(electricalResponse.ports[0].portId,QStringLiteral("continuous"));
        QVERIFY(std::isfinite(electricalResponse.ports[0].charge));
        const double electricalChargeRelativeError=
            std::abs(electricalResponse.ports[0].charge+expectedAffineCharge)/std::abs(expectedAffineCharge);
        QVERIFY(std::isfinite(electricalChargeRelativeError));
        QVERIFY(electricalChargeRelativeError<1e-10);
        qInfo().noquote()<<QStringLiteral("MITC4 continuum plate nx=%1 interiorGridSamples=%2 sampledRelativeL2Error=%3 electrodeCharge=%4 capacitance=%5")
            .arg(segments)
            .arg(sampleGridSide*sampleGridSide)
            .arg(QString::number(relativeErrors[refinement],'g',17))
            .arg(QString::number(electricalResponse.ports[0].charge,'g',17))
            .arg(QString::number(mesh.operators.capacitance[0],'g',17));
    }
    for(double relativeError:relativeErrors)QVERIFY(std::isfinite(relativeError));
    for(double charge:affineCharges)QVERIFY(std::isfinite(charge));
    QVERIFY(relativeErrors[1]<.5*relativeErrors[0]);
    QVERIFY(relativeErrors[2]<.5*relativeErrors[1]);
    QVERIFY(relativeErrors[3]<.5*relativeErrors[2]);
    QVERIFY(relativeErrors[3]<.02); // Explicit 2% cap on the 16x16 fixed-grid relative L2 error.
    const double chargeScale=std::abs(expectedAffineCharge);
    QVERIFY(std::isfinite(chargeScale));
    for(size_t refinement=1;refinement<affineCharges.size();++refinement)
    {
        const double relativeChargeDrift=std::abs(affineCharges[refinement]-affineCharges[0])/chargeScale;
        QVERIFY(std::isfinite(relativeChargeDrift));
        QVERIFY(relativeChargeDrift<1e-10);
    }
}
void TestStructuralPiezoLayeredMesh::malformedInputs()
{
    const auto assemble=[](const std::vector<Node> &n,const std::vector<Member> &m,const RStructuralPiezoPlanarFrame &f=RStructuralPiezoPlanarFrame{}) {return RStructuralPiezoLayeredMesh::assemble(Family::Beam,f,n,m);};
    for(uint invalidId:{0u,std::numeric_limits<uint>::max()})
    {
        auto invalidNodes=beamNodes(); invalidNodes[0].physicalId=invalidId;
        QVERIFY_EXCEPTION_THROWN(assemble(invalidNodes,beamMembers()),std::invalid_argument);
        auto invalidMembers=beamMembers(); invalidMembers[0].physicalId=invalidId;
        QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),invalidMembers),std::invalid_argument);
    }
    auto n=beamNodes(); n.push_back(n.front()); QVERIFY_EXCEPTION_THROWN(assemble(n,beamMembers()),std::invalid_argument);
    n=beamNodes(); n.push_back({99,{{8,0,0}}}); QVERIFY_EXCEPTION_THROWN(assemble(n,beamMembers()),std::invalid_argument);
    n=beamNodes(); n[0].position[1]=.1; QVERIFY_EXCEPTION_THROWN(assemble(n,beamMembers()),std::invalid_argument);
    // A large in-plane origin offset cannot hide a 10 m normal separation.
    n={{10,{{1e12,0,0}}},{20,{{1e12+1,10,0}}}};
    QVERIFY_EXCEPTION_THROWN(assemble(n,{{4,{10,20},{layer()}}}),std::invalid_argument);
    n[1].position[1]=0;
    const auto distant=assemble(n,{{4,{10,20},{layer()}}});
    const auto local=assemble({{10,{{0,0,0}}},{20,{{1,0,0}}}},{{4,{10,20},{layer()}}});
    QCOMPARE(distant.operators.mass,local.operators.mass);
    QCOMPARE(distant.operators.stiffness,local.operators.stiffness);
    n[0].position[1]=10; n[1].position[1]=10;
    QVERIFY_EXCEPTION_THROWN(assemble(n,{{4,{10,20},{layer()}}}),std::invalid_argument);
    n=beamNodes(); n[0].position[0]=std::numeric_limits<double>::quiet_NaN(); QVERIFY_EXCEPTION_THROWN(assemble(n,beamMembers()),std::invalid_argument);
    auto m=beamMembers(); m.push_back(m.front()); QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),m),std::invalid_argument);
    m=beamMembers(); m[0].nodeIds={20,20}; QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),m),std::invalid_argument);
    m=beamMembers(); m[0].nodeIds={20,99}; QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),m),std::invalid_argument);
    m=beamMembers(); m[0].layers.push_back(m[0].layers.front()); QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),m),std::invalid_argument);
    m=beamMembers(); m.push_back({99,{10,20},{layer()}}); QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),m),std::invalid_argument);
    auto f=RStructuralPiezoPlanarFrame{}; f.ez[2]=-1; QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),beamMembers(),f),std::invalid_argument);
    f={}; f.ex[0]=2; QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),beamMembers(),f),std::invalid_argument);
    f={}; f.origin[0]=std::numeric_limits<double>::infinity(); QVERIFY_EXCEPTION_THROWN(assemble(beamNodes(),beamMembers(),f),std::invalid_argument);
    const auto mesh=assemble(beamNodes(),beamMembers());
    QVERIFY_EXCEPTION_THROWN(mesh.reconstructPhysical({}),std::invalid_argument);
    Vector x(18,0); x[0]=std::numeric_limits<double>::infinity(); QVERIFY_EXCEPTION_THROWN(mesh.projectPhysical(x,Domain::Velocity),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(mesh.projectPhysical(Vector(18,0),Domain::Velocity,-1),std::invalid_argument);
}
QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredMesh)
#include "tst_structural_piezo_layered_mesh.moc"
