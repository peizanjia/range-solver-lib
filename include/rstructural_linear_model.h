#ifndef RSTRUCTURAL_LINEAR_MODEL_H
#define RSTRUCTURAL_LINEAR_MODEL_H

#include <QJsonObject>
class RModel;

//! Bounded, side-effect-free SI operating-point export for Stage 9.
//! Fixed-support TRUSS1/BEAM1, positive mass/stiffness, nodal force inputs
//! and translation q/v/a sensors. Cable exports use the converged taut
//! tangent and the same reference-geometry lumped mass as cable dynamics.
class RStructuralLinearModel
{
public:
    static constexpr unsigned MaxEquations = 256;
    static QJsonObject build(const RModel &model);
};
#endif
