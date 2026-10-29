#ifndef RSTRUCTURAL_PIEZO_ELECTRICAL_SYSTEM_H
#define RSTRUCTURAL_PIEZO_ELECTRICAL_SYSTEM_H

#include "rstructural_piezo_element_operators.h"

enum class RStructuralPiezoElectricalRole
{
    VoltageActuator,
    FixedChargeOpenCircuit,
    ShortCircuitChargeSensor,
    IdealChargeAmplifier
};

struct RStructuralPiezoElectricalBoundary
{
    QString portId;
    RStructuralPiezoElectricalRole role = RStructuralPiezoElectricalRole::VoltageActuator;
    double fixedCharge = 0.0; //!< Absolute positive-electrode charge, C; open circuit only.
    double feedbackCapacitance = 0.0; //!< Cf > 0, F; ideal amplifier only.
    double referenceCharge = 0.0; //!< Explicit Qref, C; ideal amplifier only.
};

struct RStructuralPiezoElectricalPortResponse
{
    QString portId;
    RStructuralPiezoElectricalRole role = RStructuralPiezoElectricalRole::VoltageActuator;
    double voltage = 0.0; //!< phi(+) - phi(-), V; actuator voltage is an action record.
    double charge = 0.0; //!< Free positive-electrode charge, C.
    double deltaVoltage = 0.0;
    double deltaCharge = 0.0;
    double amplifierVoltage = 0.0; //!< -(Q-Qref)/Cf, V; zero for other roles.
    double deltaAmplifierVoltage = 0.0;
};

struct RStructuralPiezoElectricalResponse
{
    std::vector<double> internalForce; //!< K^E q + G V, N (or N m for rotations).
    std::vector<double> deltaInternalForce;
    std::vector<RStructuralPiezoElectricalPortResponse> ports; //!< Original operator port order.
};

struct RStructuralPiezoElectricalCondensedOperators
{
    uint dofs = 0;
    std::vector<double> mass;
    std::vector<double> stiffness; //!< K^E + Gs Css^-1 Gs^T; ONLY open sensors condensed.
    std::vector<QString> actuatorPortIds; //!< Original operator order, filtered to actuators.
    std::vector<double> coupling; //!< dofs*actuators, Ga; isolated-port first release.
    std::vector<double> fixedChargeForce; //!< Gs Css^-1 Qfixed, absolute force offset.
};

//! Validated, immutable electrical boundary for shared native FE operators.
//! Matrices are row-major SI. Equations: force=K^E q+G V, Q=Ce V-G^T q.
//! Every exact port ID needs one role. Nonzero inter-port capacitance is rejected.
//! Amplifier electrode is virtual ground (V=0); its output uses inverting polarity
//! -(Q-Qref)/Cf. This ideal model has no leakage, input impedance or bandwidth.
//! q0 and Va0 freeze a reference for increments, not a certified force equilibrium:
//! the native plant must certify its external-force residual separately.
//! All invalid/nonrepresentable inputs throw std::invalid_argument; evaluateInto
//! commits its workspace only after the complete response has been validated.
class RStructuralPiezoElectricalSystem
{
public:
    RStructuralPiezoElectricalSystem(
        const RStructuralPiezoElementOperators &operators,
        const std::vector<RStructuralPiezoElectricalBoundary> &boundaries,
        const std::vector<double> &q0,
        const std::vector<double> &actuatorVoltages0);

    const RStructuralPiezoElectricalCondensedOperators &condensed() const;
    const std::vector<double> &referenceDisplacement() const;
    const std::vector<double> &referenceActuatorVoltages() const;
    const RStructuralPiezoElectricalResponse &referenceResponse() const;

    //! q and Va are ABSOLUTE; Va follows condensed().actuatorPortIds order.
    RStructuralPiezoElectricalResponse evaluate(
        const std::vector<double> &q,const std::vector<double> &actuatorVoltages) const;
    void evaluateInto(const std::vector<double> &q,
                      const std::vector<double> &actuatorVoltages,
                      RStructuralPiezoElectricalResponse &workspace) const;

private:
    RStructuralPiezoElementOperators operators;
    std::vector<RStructuralPiezoElectricalBoundary> boundaries;
    std::vector<double> q0;
    std::vector<double> actuatorVoltages0;
    RStructuralPiezoElectricalCondensedOperators condensedOperators;
    RStructuralPiezoElectricalResponse initialResponse;
    RStructuralPiezoElectricalResponse absoluteResponse(
        const std::vector<double> &q,const std::vector<double> &actuatorVoltages) const;
};

#endif // RSTRUCTURAL_PIEZO_ELECTRICAL_SYSTEM_H
