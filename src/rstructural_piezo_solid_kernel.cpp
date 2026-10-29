#include "rstructural_piezo_solid_kernel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
using V3=std::array<double,3>;
void require(bool ok,const char *s) {if(!ok) throw std::invalid_argument(s);}
V3 sub(const V3 &a,const V3 &b) {return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
double norm(const V3 &a) {return std::hypot(std::hypot(a[0],a[1]),a[2]);}
struct Shape {std::vector<double> n; std::vector<V3> gradient; double determinant=0;};
Shape shape(const std::vector<V3> &x,const V3 &parent) {
    const size_t nn=x.size(); Shape s; s.n.resize(nn); s.gradient.resize(nn); std::vector<V3> dn(nn);
    if(nn==4) {s.n={1-parent[0]-parent[1]-parent[2],parent[0],parent[1],parent[2]}; dn={{{-1,-1,-1}},{{1,0,0}},{{0,1,0}},{{0,0,1}}};}
    else {
        constexpr int signs[8][3]={{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
        for(size_t a=0;a<8;++a) {
            s.n[a]=.125; for(size_t k=0;k<3;++k) s.n[a]*=1+signs[a][k]*parent[k];
            for(size_t k=0;k<3;++k) {dn[a][k]=.125*signs[a][k]; for(size_t j=0;j<3;++j) if(j!=k) dn[a][k]*=1+signs[a][j]*parent[j];}
        }
    }
    // J physical rows / parent columns. Relative origin retains translation invariance.
    double j[3][3]{};
    for(size_t a=0;a<nn;++a) for(size_t i=0;i<3;++i) for(size_t k=0;k<3;++k) j[i][k]+=(x[a][i]-x[0][i])*dn[a][k];
    double co[3][3];
    for(size_t i=0;i<3;++i) for(size_t k=0;k<3;++k)
        co[i][k]=j[(i+1)%3][(k+1)%3]*j[(i+2)%3][(k+2)%3]-j[(i+1)%3][(k+2)%3]*j[(i+2)%3][(k+1)%3];
    s.determinant=j[0][0]*co[0][0]+j[0][1]*co[0][1]+j[0][2]*co[0][2];
    double scale=1; for(size_t k=0;k<3;++k) scale*=std::hypot(std::hypot(j[0][k],j[1][k]),j[2][k]);
    require(std::isfinite(scale) && scale>0 && std::isfinite(s.determinant) && s.determinant>1e-12*scale,
            "Solid requires a positive, resolved Jacobian (inverted or degenerate element)");
    for(size_t a=0;a<nn;++a) for(size_t i=0;i<3;++i) {
        for(size_t k=0;k<3;++k) s.gradient[a][i]+=co[i][k]/s.determinant*dn[a][k];
        require(std::isfinite(s.gradient[a][i]),"Solid derivative overflow");
    }
    return s;
}
}

RStructuralPiezoSolidNodalOperators RStructuralPiezoSolidKernel::assemble(
    const std::vector<V3> &nodes,const RStructuralPiezoContinuumMaterial &material) {
    require(nodes.size()==4 || nodes.size()==8,"Solid supports only TET4 and affine HEX8"); material.validate();
    for(const auto &x:nodes) for(double v:x) require(std::isfinite(v),"Non-finite solid geometry");
    const size_t nn=nodes.size(),nd=3*nn;
    if(nn==8) {
        const V3 a=sub(nodes[1],nodes[0]),b=sub(nodes[3],nodes[0]),c=sub(nodes[4],nodes[0]);
        const double span=std::max({norm(a),norm(b),norm(c)});
        constexpr int bits[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for(size_t n=0;n<8;++n) {
            V3 residual=sub(nodes[n],nodes[0]);
            for(size_t k=0;k<3;++k) residual[k]-=bits[n][0]*a[k]+bits[n][1]*b[k]+bits[n][2]*c[k];
            require(norm(residual)<=1e-10*span,"HEX8 requires an affine parallelepiped; warped HEX8 is unsupported");
        }
    }
    RStructuralPiezoSolidNodalOperators result; auto &op=result.mechanical;
    op.dofs=uint(nd); op.mass.assign(nd*nd,0);op.stiffness.assign(nd*nd,0);
    result.coupling.assign(nd*nn,0);result.dielectric.assign(nn*nn,0);
    std::vector<std::pair<V3,double>> quadrature;
    if(nn==4) {
        const double a=(5+3*std::sqrt(5.))/20,b=(5-std::sqrt(5.))/20;
        quadrature={{{{b,b,b}},1./24},{{{a,b,b}},1./24},{{{b,a,b}},1./24},{{{b,b,a}},1./24}};
    } else {const double g=1/std::sqrt(3.);for(double x:{-g,g})for(double y:{-g,g})for(double z:{-g,g})quadrature.push_back({{{x,y,z}},1});}
    for(const auto &point:quadrature) {
        const auto s=shape(nodes,point.first); const double w=s.determinant*point.second;
        std::vector<double> B(6*nd);
        for(size_t a=0;a<nn;++a) {
            const auto &d=s.gradient[a];const size_t i=3*a;
            B[i]=d[0];B[nd+i+1]=d[1];B[2*nd+i+2]=d[2];
            B[3*nd+i+1]=d[2];B[3*nd+i+2]=d[1];
            B[4*nd+i]=d[2];B[4*nd+i+2]=d[0];B[5*nd+i]=d[1];B[5*nd+i+1]=d[0];
            for(size_t b=0;b<nn;++b) for(size_t k=0;k<3;++k) op.mass[(3*a+k)*nd+3*b+k]+=w*material.density*s.n[a]*s.n[b];
        }
        for(size_t i=0;i<nd;++i) {
            for(size_t j=0;j<nd;++j) for(size_t a=0;a<6;++a) for(size_t b=0;b<6;++b)
                op.stiffness[i*nd+j]+=w*B[a*nd+i]*material.elasticity[a*6+b]*B[b*nd+j];
            if(material.piezoelectric) for(size_t j=0;j<nn;++j) for(size_t a=0;a<6;++a) for(size_t k=0;k<3;++k)
                result.coupling[i*nn+j]+=w*B[a*nd+i]*material.coupling[k*6+a]*s.gradient[j][k];
        }
        if(material.piezoelectric) for(size_t i=0;i<nn;++i) for(size_t j=0;j<nn;++j) for(size_t a=0;a<3;++a) for(size_t b=0;b<3;++b)
            result.dielectric[i*nn+j]+=w*s.gradient[i][a]*material.permittivity[a*3+b]*s.gradient[j][b];
        op.totalMass+=w*material.density;
    }
    op.validate();
    for(const auto *v:{&result.coupling,&result.dielectric})for(double x:*v)require(std::isfinite(x),"Solid electrical operator overflow");
    return result;
}
