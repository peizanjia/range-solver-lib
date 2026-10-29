#ifndef RSTRUCTURAL_PIEZO_ELEMENT_OPERATORS_H
#define RSTRUCTURAL_PIEZO_ELEMENT_OPERATORS_H

#include <array>
#include <QString>
#include <vector>

//! Engineering Voigt order xx, yy, zz, yz, xz, xy; matrices are row-major SI.
//! Material axes are already aligned to the element's local x/y/z frame.
struct RStructuralPiezoContinuumMaterial
{
    std::array<double,36> elasticity{};
    std::array<double,18> coupling{}; //!< e, electric rows x/y/z, C/m^2.
    std::array<double,9> permittivity{}; //!< epsilon^S, F/m.
    double density = 0.0;
    bool piezoelectric = false;
    static RStructuralPiezoContinuumMaterial isotropic(double young,double poisson,double density);
    //! Engineering strain-charge data s^E,d,epsilon^T -> c^E,e,epsilon^S (SI).
    static RStructuralPiezoContinuumMaterial fromStrainCharge(
        const std::array<double,36> &compliance,const std::array<double,18> &strainCoupling,
        const std::array<double,9> &stressFreePermittivity,double density);
    void validate() const;
};

struct RStructuralPiezoLayer
{
    RStructuralPiezoContinuumMaterial material;
    double lowerZ = 0.0;
    double upperZ = 0.0;
    double width = 1.0; //!< Beam width; plate requires the default (geometry owns area).
    QString portId; //!< Continuous isolated electrode pair; passive layer has none.
    int polarity = 1; //!< +1 or -1, positive electrode top minus bottom.
    std::array<double,4> coverage{{-1.0,1.0,-1.0,1.0}}; //!< Parent xi/eta rectangle.
    void validate() const;
};

struct RStructuralPiezoReducedConstitutive
{
    std::vector<uint> components;
    std::vector<double> elasticity;
    std::vector<double> coupling; //!< Selected electric row after free-stress elimination.
    double permittivity = 0.0;
};

//! Eliminate all unretained engineering stresses using a constitutive Schur complement.
RStructuralPiezoReducedConstitutive rStructuralPiezoReduceFreeStress(
    const RStructuralPiezoContinuumMaterial &material,const std::vector<uint> &retained,
    uint electricAxis=2);

//! Scalar two-node axial replacement; material local x is member/electric axis.
//! Native assembly owns incidence signs. fixedMass replaces the original member
//! mass once, rather than adding piezo material on top of a retained host.
struct RStructuralPiezoAxialMemberProperties
{
    double k=0.0,g=0.0,Ce=0.0,fixedMass=0.0;
};
RStructuralPiezoAxialMemberProperties rStructuralPiezoCompileAxialMember(
    const RStructuralPiezoContinuumMaterial &material,double area,double length);

//! Same operators feed actual native FE dynamics, modes and controller/export paths.
//! M/K are dofs*dofs; G is dofs*ports; Ce is ports*ports. No electrical inertia.
struct RStructuralPiezoElementOperators
{
    uint dofs = 0;
    std::vector<QString> portIds;
    std::vector<double> mass;
    std::vector<double> stiffness;
    std::vector<double> coupling;
    std::vector<double> capacitance;
    double totalMass = 0.0;
    void validate() const;
};

class RStructuralPiezoBeamLayerKernel
{
public:
    //! Planar local DOFs: u1,w1,thetaY1,u2,w2,thetaY2; thetaY=-dw/dx.
    //! Exact supported layer inertia includes u(z)=u0+z*thetaY and w=w0.
    static RStructuralPiezoElementOperators assemble(
        double length,const std::vector<RStructuralPiezoLayer> &layers);
};

class RStructuralPiezoPlateMitc4Kernel
{
public:
    //! Counter-clockwise convex planar Q4, local xy coordinates.
    //! Per-node DOFs u,v,w,thetaX,thetaY; u(z)=u+z*thetaY,v(z)=v-z*thetaX.
    //! MITC4 covariant shear tying, complete membrane/bending/shear and I0/I1/I2.
    static RStructuralPiezoElementOperators assemble(
        const std::array<std::array<double,2>,4> &nodes,
        const std::vector<RStructuralPiezoLayer> &layers,
        double shearCorrection=5.0/6.0);
};

#endif // RSTRUCTURAL_PIEZO_ELEMENT_OPERATORS_H
