#include "rstructural_piezo_element_operators.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
constexpr unsigned Dofs = 20;
using Nodes = std::array<std::array<double,2>,4>;
using Row = std::array<double,Dofs>;
using Strains = std::array<Row,5>;

struct Shape
{
    std::array<double,4> n{},xi{},eta{},dx{},dy{};
    double xXi=0.0,yXi=0.0,xEta=0.0,yEta=0.0,det=0.0;
};

Shape shape(const Nodes &nodes,double xi,double eta)
{
    constexpr double sx[4]={-1.0,1.0,1.0,-1.0};
    constexpr double sy[4]={-1.0,-1.0,1.0,1.0};
    Shape s;
    for (unsigned i=0;i<4;++i)
    {
        s.n[i]=0.25*(1.0+sx[i]*xi)*(1.0+sy[i]*eta);
        s.xi[i]=0.25*sx[i]*(1.0+sy[i]*eta);
        s.eta[i]=0.25*sy[i]*(1.0+sx[i]*xi);
        // Relative coordinates avoid cancellation from a large translated origin.
        const double x=nodes[i][0]-nodes[0][0],y=nodes[i][1]-nodes[0][1];
        s.xXi+=s.xi[i]*x; s.yXi+=s.xi[i]*y;
        s.xEta+=s.eta[i]*x; s.yEta+=s.eta[i]*y;
    }
    s.det=s.xXi*s.yEta-s.yXi*s.xEta;
    if (!std::isfinite(s.det) || s.det<=0.0)
        throw std::invalid_argument("MITC4 plate requires a finite positive Jacobian");
    for (unsigned i=0;i<4;++i)
    {
        s.dx[i]=(s.yEta*s.xi[i]-s.yXi*s.eta[i])/s.det;
        s.dy[i]=(-s.xEta*s.xi[i]+s.xXi*s.eta[i])/s.det;
        if (!std::isfinite(s.dx[i]) || !std::isfinite(s.dy[i]))
            throw std::invalid_argument("MITC4 plate has non-finite physical derivatives");
    }
    return s;
}

Row covariantShear(const Shape &s,bool xi)
{
    Row row{};
    for (unsigned i=0;i<4;++i)
    {
        row[5*i+2]=xi?s.xi[i]:s.eta[i];
        row[5*i+3]=-(xi?s.yXi:s.yEta)*s.n[i];
        row[5*i+4]=(xi?s.xXi:s.xEta)*s.n[i];
    }
    return row;
}

void validateGeometry(const Nodes &nodes)
{
    for (const auto &node:nodes)
        for (double value:node)
            if (!std::isfinite(value))
                throw std::invalid_argument("MITC4 plate coordinates must be finite");
    for (unsigned i=0;i<4;++i)
    {
        const auto &a=nodes[i],&b=nodes[(i+1)%4],&c=nodes[(i+2)%4];
        const double cross=(b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]);
        if (!std::isfinite(cross) || cross<=0.0)
            throw std::invalid_argument("MITC4 plate nodes must form a strictly convex counter-clockwise Q4");
    }
    // For a bilinear planar map det(J) is affine; positive corners cover the
    // entire parent square, including every layer subcell and tying point.
    for (double xi:{-1.0,1.0})
        for (double eta:{-1.0,1.0})
            shape(nodes,xi,eta);
}

bool overlaps(const RStructuralPiezoLayer &a,const RStructuralPiezoLayer &b)
{
    return std::max(a.lowerZ,b.lowerZ)<std::min(a.upperZ,b.upperZ)
        && std::max(a.coverage[0],b.coverage[0])<std::min(a.coverage[1],b.coverage[1])
        && std::max(a.coverage[2],b.coverage[2])<std::min(a.coverage[3],b.coverage[3]);
}
}

RStructuralPiezoElementOperators RStructuralPiezoPlateMitc4Kernel::assemble(
    const Nodes &nodes,const std::vector<RStructuralPiezoLayer> &layers,double shearCorrection)
{
    validateGeometry(nodes);
    if (!std::isfinite(shearCorrection) || shearCorrection<=0.0 || shearCorrection>1.0)
        throw std::invalid_argument("MITC4 shear correction must be in (0,1]");
    if (layers.empty())
        throw std::invalid_argument("MITC4 plate needs at least one layer");

    RStructuralPiezoElementOperators result;
    result.dofs=Dofs;
    for (std::size_t i=0;i<layers.size();++i)
    {
        const auto &layer=layers[i];
        layer.validate();
        if (layer.width!=1.0)
            throw std::invalid_argument("MITC4 layer width must be 1; the plate geometry defines area");
        for (std::size_t j=0;j<i;++j)
            if (overlaps(layer,layers[j]))
                throw std::invalid_argument("MITC4 plate layers overlap in both coverage and thickness");
        if (!layer.portId.isEmpty()
            && std::find(result.portIds.begin(),result.portIds.end(),layer.portId)==result.portIds.end())
            result.portIds.push_back(layer.portId);
    }
    const std::size_t ports=result.portIds.size();
    result.mass.assign(Dofs*Dofs,0.0);
    result.stiffness.assign(Dofs*Dofs,0.0);
    result.coupling.assign(Dofs*ports,0.0);
    result.capacitance.assign(ports*ports,0.0);

    // Independent small-strain planar specialization of covariant MITC4 tying:
    // Dvorkin & Bathe (1984), Engineering Computations 1, pp. 77-88, eqs. 3/6.
    // xi shear is tied at (0,+/-1); eta shear at (+/-1,0). Kinematics
    // are always evaluated on the entire host, even for partial electrodes.
    const Row a=covariantShear(shape(nodes,0.0,-1.0),true);
    const Row c=covariantShear(shape(nodes,0.0,1.0),true);
    const Row d=covariantShear(shape(nodes,-1.0,0.0),false);
    const Row b=covariantShear(shape(nodes,1.0,0.0),false);
    const double shearScale=std::sqrt(shearCorrection);
    const double gauss=1.0/std::sqrt(3.0);
    for (const auto &layer:layers)
    {
        // sigma_zz=0: retain xx, yy, xy, xz, yz in engineering notation.
        const auto material=rStructuralPiezoReduceFreeStress(layer.material,{0,1,5,4,3});
        const bool active=!layer.portId.isEmpty();
        const std::size_t port=active?std::size_t(std::find(result.portIds.begin(),result.portIds.end(),layer.portId)
                                                 -result.portIds.begin()):0;
        const double thickness=layer.upperZ-layer.lowerZ;
        const double middle=layer.lowerZ+0.5*thickness;
        const double i0=layer.material.density*thickness;
        const double i1=i0*middle;
        const double i2=i0*(middle*middle+thickness*thickness/12.0);
        const double hx=0.5*(layer.coverage[1]-layer.coverage[0]);
        const double hy=0.5*(layer.coverage[3]-layer.coverage[2]);
        const double cx=layer.coverage[0]+hx,cy=layer.coverage[2]+hy;
        if (!std::isfinite(thickness) || thickness<=0.0 || !std::isfinite(middle)
            || !std::isfinite(i0) || i0<=0.0 || !std::isfinite(i1)
            || !std::isfinite(i2) || i2<=0.0 || hx<=0.0 || hy<=0.0)
            throw std::invalid_argument("MITC4 layer has unrepresentable thickness or inertia");
        for (double gx:{-gauss,gauss})
            for (double gy:{-gauss,gauss})
            {
                const double xi=cx+hx*gx,eta=cy+hy*gy;
                const Shape s=shape(nodes,xi,eta);
                const double areaWeight=s.det*hx*hy;
                if (!std::isfinite(areaWeight) || areaWeight<=0.0
                    || !std::isfinite(areaWeight*thickness) || areaWeight*thickness<=0.0)
                    throw std::invalid_argument("MITC4 layer has unrepresentable integration volume");
                result.totalMass+=areaWeight*i0;
                // Exact section moments and consistent bilinear shape mass.
                // No drilling coordinate is present or supplied with artificial inertia.
                for (unsigned i=0;i<4;++i)
                    for (unsigned j=0;j<4;++j)
                    {
                        const double nn=areaWeight*s.n[i]*s.n[j];
                        for (unsigned k=0;k<3;++k)
                            result.mass[(5*i+k)*Dofs+5*j+k]+=nn*i0;
                        result.mass[(5*i+3)*Dofs+5*j+3]+=nn*i2;
                        result.mass[(5*i+4)*Dofs+5*j+4]+=nn*i2;
                        result.mass[(5*i)*Dofs+5*j+4]+=nn*i1;
                        result.mass[(5*i+4)*Dofs+5*j]+=nn*i1;
                        result.mass[(5*i+1)*Dofs+5*j+3]-=nn*i1;
                        result.mass[(5*i+3)*Dofs+5*j+1]-=nn*i1;
                    }
                Strains strain{};
                for (unsigned j=0;j<Dofs;++j)
                {
                    const double covXi=0.5*((1.0-eta)*a[j]+(1.0+eta)*c[j]);
                    const double covEta=0.5*((1.0-xi)*d[j]+(1.0+xi)*b[j]);
                    // J has rows (x_,xi,y_,xi), (x_,eta,y_,eta).
                    // The same inverse maps both covariant shear and shape gradients.
                    strain[3][j]=shearScale*(s.yEta*covXi-s.yXi*covEta)/s.det;
                    strain[4][j]=shearScale*(-s.xEta*covXi+s.xXi*covEta)/s.det;
                }
                if (active)
                    result.capacitance[port*ports+port]+=areaWeight*(material.permittivity/thickness);
                // Two through-thickness Gauss points integrate every z^0/z^1/z^2
                // term, including anisotropic normal/shear constitutive couplings.
                for (double gz:{-gauss,gauss})
                {
                    const double z=middle+0.5*thickness*gz;
                    const double weight=areaWeight*0.5*thickness;
                    if (!std::isfinite(weight) || weight<=0.0 || !std::isfinite(z))
                        throw std::invalid_argument("MITC4 layer has unrepresentable thickness quadrature");
                    for (unsigned i=0;i<4;++i)
                    {
                        strain[0][5*i]=s.dx[i];
                        strain[0][5*i+4]=z*s.dx[i];
                        strain[1][5*i+1]=s.dy[i];
                        strain[1][5*i+3]=-z*s.dy[i];
                        strain[2][5*i]=s.dy[i];
                        strain[2][5*i+1]=s.dx[i];
                        strain[2][5*i+3]=-z*s.dx[i];
                        strain[2][5*i+4]=z*s.dy[i];
                    }
                    for (unsigned i=0;i<Dofs;++i)
                    {
                        if (active)
                            for (unsigned k=0;k<5;++k)
                                // dz/t=1/2 per point; cancel before dividing to
                                // avoid overflow for a thin, finite electrode.
                                result.coupling[i*ports+port]+=areaWeight*0.5*strain[k][i]
                                    *material.coupling[k]*layer.polarity;
                        for (unsigned j=i;j<Dofs;++j)
                        {
                            double value=0.0;
                            for (unsigned k=0;k<5;++k)
                                for (unsigned l=0;l<5;++l)
                                    value+=strain[k][i]*material.elasticity[k*5+l]*strain[l][j];
                            result.stiffness[i*Dofs+j]+=weight*value;
                        }
                    }
                }
            }
    }
    for (unsigned i=0;i<Dofs;++i)
        for (unsigned j=i+1;j<Dofs;++j)
            result.stiffness[j*Dofs+i]=result.stiffness[i*Dofs+j];
    result.validate();
    return result;
}
