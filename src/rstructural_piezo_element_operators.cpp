#include "rstructural_piezo_element_operators.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
constexpr double symmetryTolerance = 1.0e-11;

void require(bool condition,const char *message)
{
    if (!condition) throw std::invalid_argument(message);
}

template<class Values>
void requireFinite(const Values &values)
{
    for (double value : values)
        require(std::isfinite(value),"Piezo operator contains a non-finite entry");
}

void requireSymmetric(const std::vector<double> &matrix,size_t n)
{
    require(matrix.size() == n*n,"Piezo matrix size mismatch");
    requireFinite(matrix);
    for (size_t i=0;i<n;++i)
        for (size_t j=0;j<i;++j)
        {
            const double a = matrix[i*n+j], b = matrix[j*n+i];
            const double diagonalScale = std::sqrt(std::abs(matrix[i*n+i]))*
                                         std::sqrt(std::abs(matrix[j*n+j]));
            const double scale = std::max({std::abs(a),std::abs(b),diagonalScale});
            require(std::abs(a-b) <= symmetryTolerance*scale,
                    "Piezo matrix must be symmetric");
        }
}

// Normalize each coordinate independently: dielectric, elastic and rotational
// inertia entries have very different SI magnitudes. There is no absolute floor.
struct PositiveFactor
{
    size_t n;
    std::vector<double> scale;
    std::vector<double> lower;

    explicit PositiveFactor(const std::vector<double> &matrix,size_t size)
        : n(size),scale(n),lower(n*n,0.0)
    {
        requireSymmetric(matrix,n);
        for (size_t i=0;i<n;++i)
        {
            require(matrix[i*n+i] > 0.0,"Piezo matrix must be positive definite");
            scale[i] = std::sqrt(matrix[i*n+i]);
        }
        for (size_t i=0;i<n;++i)
            for (size_t j=0;j<=i;++j)
            {
                double value = matrix[i*n+j]/scale[i]/scale[j];
                for (size_t k=0;k<j;++k) value -= lower[i*n+k]*lower[j*n+k];
                require(std::isfinite(value),"Piezo factorization is non-finite");
                if (i == j)
                {
                    require(value > 0.0,"Piezo matrix must be positive definite");
                    lower[i*n+j] = std::sqrt(value);
                }
                else lower[i*n+j] = value/lower[j*n+j];
            }
    }

    std::vector<double> solve(const std::vector<double> &rhs) const
    {
        require(rhs.size() == n,"Piezo linear solve size mismatch");
        std::vector<double> solution(n);
        for (size_t i=0;i<n;++i)
        {
            double value = rhs[i]/scale[i];
            for (size_t j=0;j<i;++j) value -= lower[i*n+j]*solution[j];
            solution[i] = value/lower[i*n+i];
        }
        for (size_t i=n;i-- > 0;)
        {
            double value = solution[i];
            for (size_t j=i+1;j<n;++j) value -= lower[j*n+i]*solution[j];
            solution[i] = value/lower[i*n+i];
        }
        for (size_t i=0;i<n;++i) solution[i] /= scale[i];
        requireFinite(solution);
        return solution;
    }
};

bool validId(const QString &id)
{
    return !id.trimmed().isEmpty() && !id.contains(QChar(0));
}

//! Sum products at a shared binary exponent, avoiding intermediate overflow
//! even on MSVC where long double has the same range as double.
double productSum(const std::vector<double> &a,const std::vector<double> &b)
{
    require(a.size()==b.size(),"Piezo product sum count mismatch");
    std::vector<double> mantissa(a.size()); std::vector<int> exponent(a.size());
    int maximum=std::numeric_limits<int>::min();
    for(size_t i=0;i<a.size();++i)
    {
        if(a[i]==0 || b[i]==0) continue;
        int ea=0,eb=0; mantissa[i]=std::frexp(a[i],&ea)*std::frexp(b[i],&eb); exponent[i]=ea+eb;
        maximum=std::max(maximum,exponent[i]);
    }
    if(maximum==std::numeric_limits<int>::min()) return 0;
    double value=0;
    for(size_t i=0;i<a.size();++i) if(mantissa[i]!=0) value+=std::scalbn(mantissa[i],exponent[i]-maximum);
    const double result=std::scalbn(value,maximum);
    require(std::isfinite(result),"Piezo constitutive conversion exceeds finite range");
    return result;
}

double physicalProduct(double a,double b,double c,double divisor=1.0)
{
    if(a==0 || b==0 || c==0) return 0;
    int ea=0,eb=0,ec=0,ed=0;
    const double value=std::frexp(a,&ea)*std::frexp(b,&eb)*std::frexp(c,&ec)/std::frexp(divisor,&ed);
    const double result=std::scalbn(value,ea+eb+ec-ed);
    require(std::isfinite(result)&&result!=0,"Axial piezo property overflows or underflows finite range");
    return result;
}
}

RStructuralPiezoContinuumMaterial RStructuralPiezoContinuumMaterial::fromStrainCharge(
    const std::array<double,36> &compliance,const std::array<double,18> &strainCoupling,
    const std::array<double,9> &stressFreePermittivity,double density)
{
    requireFinite(compliance); requireFinite(strainCoupling); requireFinite(stressFreePermittivity);
    const PositiveFactor complianceFactor(std::vector<double>(compliance.begin(),compliance.end()),6);
    const PositiveFactor dielectricFactor(std::vector<double>(stressFreePermittivity.begin(),stressFreePermittivity.end()),3);
    RStructuralPiezoContinuumMaterial result; result.density=density; result.piezoelectric=true;
    for(size_t j=0;j<6;++j)
    {
        std::vector<double> unit(6,0); unit[j]=1;
        const auto column=complianceFactor.solve(unit);
        for(size_t i=0;i<6;++i) result.elasticity[6*i+j]=column[i];
    }
    for(size_t i=0;i<6;++i) for(size_t j=0;j<i;++j)
        result.elasticity[6*i+j]=result.elasticity[6*j+i]=result.elasticity[6*i+j]/2+result.elasticity[6*j+i]/2;
    // Formula reference (independent implementation): Abaqus piezoelectric
    // material definition, c^E=(s^E)^-1, e=d*c^E,
    // epsilon^S=epsilon^T-d*c^E*d^T.
    // https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-piezoelect.htm
    for(size_t i=0;i<3;++i) for(size_t j=0;j<6;++j)
    {
        std::vector<double> a(6),b(6);
        for(size_t k=0;k<6;++k) {a[k]=strainCoupling[6*i+k]; b[k]=result.elasticity[6*k+j];}
        result.coupling[6*i+j]=productSum(a,b);
    }
    for(size_t i=0;i<3;++i) for(size_t j=0;j<3;++j)
    {
        std::vector<double> a(7),b(7); a[6]=stressFreePermittivity[3*i+j]; b[6]=1;
        for(size_t k=0;k<6;++k) {a[k]=-result.coupling[6*i+k]; b[k]=strainCoupling[6*j+k];}
        result.permittivity[3*i+j]=productSum(a,b);
    }
    for(size_t i=0;i<3;++i) for(size_t j=0;j<i;++j)
        result.permittivity[3*i+j]=result.permittivity[3*j+i]=result.permittivity[3*i+j]/2+result.permittivity[3*j+i]/2;
    result.validate(); return result;
}

RStructuralPiezoContinuumMaterial RStructuralPiezoContinuumMaterial::isotropic(
    double young,double poisson,double density)
{
    require(std::isfinite(young) && young > 0.0 && std::isfinite(poisson) &&
            poisson > -1.0 && poisson < 0.5,"Invalid isotropic elastic constants");
    RStructuralPiezoContinuumMaterial material;
    material.density = density;
    const double mu = young/(2.0*(1.0+poisson));
    const double lambda = young*poisson/((1.0+poisson)*(1.0-2.0*poisson));
    for (size_t i=0;i<3;++i)
        for (size_t j=0;j<3;++j)
            material.elasticity[i*6+j] = lambda+(i == j ? 2.0*mu : 0.0);
    for (size_t i=3;i<6;++i) material.elasticity[i*6+i] = mu;
    material.validate();
    return material;
}

void RStructuralPiezoContinuumMaterial::validate() const
{
    require(std::isfinite(density) && density > 0.0,"Piezo material density must be positive");
    requireFinite(coupling);
    requireFinite(permittivity);
    PositiveFactor elasticFactor(std::vector<double>(elasticity.begin(),elasticity.end()),6);
    if (!piezoelectric)
        for (double value : coupling)
            require(value == 0.0,"Passive material cannot have piezoelectric coupling");
    const bool zeroDielectric = std::all_of(permittivity.begin(),permittivity.end(),
                                         [](double value) { return value == 0.0; });
    if (piezoelectric || !zeroDielectric)
        PositiveFactor dielectricFactor(std::vector<double>(permittivity.begin(),permittivity.end()),3);
}

void RStructuralPiezoLayer::validate() const
{
    material.validate();
    require(std::isfinite(lowerZ) && std::isfinite(upperZ) && lowerZ < upperZ,
            "Piezo layer thickness must be positive and finite");
    require(std::isfinite(width) && width > 0.0,"Piezo layer width must be positive and finite");
    require(polarity == 1 || polarity == -1,"Piezo layer polarity must be +1 or -1");
    requireFinite(coverage);
    require(coverage[0] >= -1.0 && coverage[1] <= 1.0 && coverage[0] < coverage[1] &&
            coverage[2] >= -1.0 && coverage[3] <= 1.0 && coverage[2] < coverage[3],
            "Piezo layer coverage must be an ordered parent-coordinate rectangle");
    require(material.piezoelectric ? validId(portId) : portId.isEmpty(),
            "Piezo layer requires an electrode port; passive layer must have none");
}

RStructuralPiezoReducedConstitutive rStructuralPiezoReduceFreeStress(
    const RStructuralPiezoContinuumMaterial &material,const std::vector<uint> &retained,uint electricAxis)
{
    require(electricAxis<3,"Invalid piezoelectric field axis");
    material.validate();
    require(!retained.empty() && retained.size() <= 6,"Invalid retained stress components");
    std::array<bool,6> selected{};
    for (uint component : retained)
    {
        require(component < 6 && !selected[component],"Invalid or duplicate retained component");
        selected[component] = true;
    }
    std::vector<uint> eliminated;
    for (uint i=0;i<6;++i) if (!selected[i]) eliminated.push_back(i);
    const size_t a = retained.size(), b = eliminated.size();
    RStructuralPiezoReducedConstitutive result;
    result.components = retained;
    result.elasticity.resize(a*a);
    result.coupling.resize(a);
    result.permittivity = material.permittivity[4*electricAxis];
    for (size_t i=0;i<a;++i)
    {
        result.coupling[i] = material.coupling[6*electricAxis+retained[i]];
        for (size_t j=0;j<a;++j)
            result.elasticity[i*a+j] = material.elasticity[retained[i]*6+retained[j]];
    }
    if (b != 0)
    {
        std::vector<double> cbb(b*b),eb(b);
        for (size_t i=0;i<b;++i)
        {
            eb[i] = material.coupling[6*electricAxis+eliminated[i]];
            for (size_t j=0;j<b;++j)
                cbb[i*b+j] = material.elasticity[eliminated[i]*6+eliminated[j]];
        }
        const PositiveFactor factor(cbb,b);
        for (size_t j=0;j<a;++j)
        {
            std::vector<double> column(b);
            for (size_t k=0;k<b;++k)
                column[k] = material.elasticity[eliminated[k]*6+retained[j]];
            const auto solution = factor.solve(column);
            for (size_t k=0;k<b;++k)
            {
                result.coupling[j] -= eb[k]*solution[k];
                for (size_t i=0;i<a;++i)
                    result.elasticity[i*a+j] -=
                        material.elasticity[retained[i]*6+eliminated[k]]*solution[k];
            }
        }
        const auto electricSolution = factor.solve(eb);
        for (size_t k=0;k<b;++k) result.permittivity += eb[k]*electricSolution[k];
    }
    requireFinite(result.coupling);
    require(std::isfinite(result.permittivity) &&
            (material.piezoelectric ? result.permittivity > 0.0 : result.permittivity >= 0.0),
            "Reduced piezo permittivity must be finite and physical");
    PositiveFactor reducedFactor(result.elasticity,a);
    return result;
}

RStructuralPiezoAxialMemberProperties rStructuralPiezoCompileAxialMember(
    const RStructuralPiezoContinuumMaterial &material,double area,double length)
{
    require(std::isfinite(area)&&area>0&&std::isfinite(length)&&length>0,
            "Axial piezo area and length must be positive and finite");
    require(material.piezoelectric,"Axial piezo replacement requires a piezoelectric material");
    const auto axial=rStructuralPiezoReduceFreeStress(material,{0},0);
    RStructuralPiezoAxialMemberProperties result;
    result.k=physicalProduct(area,axial.elasticity[0],1,length);
    result.g=physicalProduct(area,axial.coupling[0],1,length);
    result.Ce=physicalProduct(area,axial.permittivity,1,length);
    result.fixedMass=physicalProduct(material.density,area,length);
    require(result.k>0&&result.Ce>0&&result.fixedMass>0,"Unusable axial piezo replacement properties");
    return result;
}

void RStructuralPiezoElementOperators::validate() const
{
    require(dofs > 0,"Piezo element must have mechanical DOFs");
    const size_t n = dofs, p = portIds.size();
    require(n <= std::numeric_limits<size_t>::max()/n &&
            (p == 0 || (p <= std::numeric_limits<size_t>::max()/p &&
                        n <= std::numeric_limits<size_t>::max()/p)),
            "Piezo element matrix dimensions overflow");
    require(mass.size() == n*n && stiffness.size() == n*n &&
            coupling.size() == n*p && capacitance.size() == p*p,
            "Piezo element matrix dimensions mismatch");
    require(std::isfinite(totalMass) && totalMass > 0.0,"Piezo element mass must be positive");
    for (size_t i=0;i<p;++i)
    {
        require(validId(portIds[i]),"Invalid piezo electrode port ID");
        for (size_t j=0;j<i;++j)
            require(portIds[i] != portIds[j],"Duplicate piezo electrode port ID");
    }
    requireFinite(coupling);
    PositiveFactor massFactor(mass,n);
    requireSymmetric(stiffness,n);
    // Rigid motions make unconstrained element stiffness singular. Full PSD,
    // rank and patch-test evidence belong to each element kernel's tests.
    for (size_t i=0;i<n;++i)
        require(stiffness[i*n+i] >= 0.0,"Piezo element has negative diagonal stiffness");
    if (p != 0) PositiveFactor capacitanceFactor(capacitance,p);
}
