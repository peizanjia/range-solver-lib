#ifndef RSTRUCTURAL_PIEZO_KERNEL_H
#define RSTRUCTURAL_PIEZO_KERNEL_H

#include <array>

using RStructuralPiezoMatrix2 = std::array<std::array<double,2>,2>;

//! Pure, terminal 1-D piezoelectric rod result under voltage control.
struct RStructuralPiezoKernelResult
{
    double enthalpy = 0.0;
    double mechanicalForce = 0.0;
    //! A*D at the dielectric right face with outward normal +x.
    //! This is not the conductor electrode free charge, whose sign is opposite.
    double dielectricDisplacementFlux = 0.0;
    double strain = 0.0;
    double electricField = 0.0;
    double stress = 0.0;
    double electricDisplacement = 0.0;
    double voltageFixedStiffness = 0.0;
    double openCircuitStiffness = 0.0;
    //! Hessian of enthalpy with respect to [elongation, voltage difference].
    RStructuralPiezoMatrix2 hessian{};
};

//! Terminal 1-D piezo rod electric-enthalpy kernel.
//!
//! Uses u = u_right-u_left, V = phi_right-phi_left, S = u/L, and E = -V/L.
//! The Hessian has a negative electrical diagonal because this is electric
//! enthalpy; it is not a positive-definite mechanical stiffness matrix.
class RStructuralPiezoKernel
{
    public:
        //! cE and epsS are positive; e is signed and encodes poling direction.
        static RStructuralPiezoKernelResult evaluate(
            double elongation,
            double voltageDifference,
            double cE,
            double e,
            double epsS,
            double area,
            double length);
};

#endif // RSTRUCTURAL_PIEZO_KERNEL_H
