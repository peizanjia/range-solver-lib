#include "rstructural_truss_dynamics_restart_adapter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <QCryptographicHash>
#include <rbl_error.h>

namespace
{
[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}

bool derivedEqual(double a,double b)
{
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    if (a == b) return true;
    return std::fabs(a-b) <= 128.0*std::numeric_limits<double>::epsilon()*
        std::max({std::numeric_limits<double>::min(),std::fabs(a),std::fabs(b)});
}

double inputValue(const RStructuralTrussInput &input)
{
    if (const auto *v = std::get_if<RStructuralTrussVoltage>(&input)) return v->value;
    if (const auto *q = std::get_if<RStructuralTrussCharge>(&input)) return q->value;
    if (const auto *l = std::get_if<RStructuralTrussCableLength>(&input)) return l->value;
    return 0.0;
}

void verifyPhase(const RStructuralTrussDynamicsModelAdapter &adapter,
    const RStructuralTrussDynamicsState &state,uint completed,uint cursor)
{
    if (completed > adapter.targetSteps() || cursor != adapter.commandCursor(completed) ||
        state.time != adapter.targetTime(completed) || state.systemIdentity != adapter.driver().systemIdentity())
        fail("Native mixed accepted state has a different system, canonical physical clock or command phase.");
    const auto replay = adapter.replayInputs(completed);
    const auto &expected = replay.completeInputs;
    const auto &actual = state.appliedInputs;
    if (actual.externalForces != expected.externalForces || actual.memberInputs.size() != expected.memberInputs.size())
        fail("Native mixed accepted applied inputs differ from complete frozen command replay.");
    for (std::size_t i=0;i<actual.memberInputs.size();++i)
        if (actual.memberInputs[i].index() != expected.memberInputs[i].index() ||
            inputValue(actual.memberInputs[i]) != inputValue(expected.memberInputs[i]))
            fail("Native mixed accepted member input kind/value differs from complete frozen command replay.");
    const auto &system=adapter.system();
    const auto n=system.prescribedMask.size();
    if (state.displacements.size()!=n || state.velocities.size()!=n || state.accelerations.size()!=n)
        fail("Native mixed accepted kinematics have invalid full-DOF dimensions.");
    auto prescribedQ=system.fixedPrescribedDisplacements;
    std::vector<double> prescribedV(n,0.0),prescribedA(n,0.0);
    for (const auto &motion:system.prescribedMotions)
    {
        const auto sample=motion.evaluate(state.time);
        const auto d=std::size_t(3)*motion.target.nodeId+uint(motion.target.dof);
        prescribedQ[d]=sample.displacement;prescribedV[d]=sample.velocity;prescribedA[d]=sample.acceleration;
    }
    // Native accepted supports are produced by this exact frozen evaluator,
    // not by approximate display fields. A tiny unsupported shift is still a
    // corrupted restart even when it changes no force on a slack member.
    for (std::size_t d=0;d<n;++d) if (system.prescribedMask[d])
        if (state.displacements[d]!=prescribedQ[d] || state.velocities[d]!=prescribedV[d] ||
            state.accelerations[d]!=prescribedA[d])
            fail("Native mixed accepted support kinematics differ from the exact frozen prescription.");
    // refresh validates the left state's support kinematics, inertial balance,
    // observations and ledgers. Its re-equilibrated acceleration is never used
    // as the original accepted Newmark acceleration or as a restored ledger.
    adapter.driver().refresh(state,state.appliedInputs);
    const auto initial = adapter.initialize();
    if (state.initialTotalEnergy != initial.initialTotalEnergy)
        fail("Native mixed accepted initial total energy differs from frozen source initialization.");
    if (completed == 0 && (state.displacements != initial.displacements || state.velocities != initial.velocities ||
        state.accelerations != initial.accelerations || state.reactions != initial.reactions ||
        state.cumulativeExternalWork != 0.0 || state.cumulativeSupportWork != 0.0 ||
        state.cumulativeElectricalWork != 0.0 || state.cumulativeLengthWork != 0.0))
        fail("Native mixed zero-step checkpoint differs from frozen source initialization.");
}

RStructuralTrussAssemblyResult observe(const RStructuralTrussDynamicsModelAdapter &adapter,
                                      const RStructuralTrussDynamicsState &state)
{
    const auto &system = adapter.system();
    auto assembly = RStructuralTrussAssembly::evaluate(system.referenceNodes,state.displacements,
        system.members,state.appliedInputs.memberInputs);
    if (!derivedEqual(assembly.storedEnergy,state.storedEnergy))
        fail("Native mixed accepted stored energy differs from fresh assembly.");
    double kinetic = 0.0;
    for (std::size_t node=0;node<system.nodalLumpedMass.size();++node)
        for (std::size_t axis=0;axis<3;++axis)
        {
            const auto d = 3*node+axis;
            // Same component accumulation order and arithmetic as the driver.
            const double root = std::sqrt(system.nodalLumpedMass[node])*state.velocities[d];
            kinetic += (0.5*root)*root;
            const double reaction = system.nodalLumpedMass[node]*state.accelerations[d]+
                (assembly.internalForces[d]-state.appliedInputs.externalForces[d]);
            if (!derivedEqual(reaction,state.reactions[d]))
                fail("Native mixed accepted reactions differ from fresh assembly and accepted acceleration.");
        }
    if (!derivedEqual(kinetic,state.kineticEnergy))
        fail("Native mixed accepted kinetic energy differs from frozen mass and accepted velocity.");
    return assembly;
}

void portObservations(const RStructuralTrussDynamicsModelAdapter &adapter,
    const RStructuralTrussAssemblyResult &assembly,
    std::vector<double> &voltages,std::vector<double> &charges,std::vector<double> &enthalpies)
{
    for (const auto &id : adapter.maps().orderedPortIds)
    {
        const auto &response = assembly.memberResponses.at(adapter.maps().memberByPortId.at(id));
        if (!response.electrical) fail("Native mixed independent port has no electrical assembly response.");
        voltages.push_back(response.electrical->voltage);
        charges.push_back(response.electrical->conductorCharge);
        enthalpies.push_back(response.electrical->enthalpy);
    }
}

void verifyDerivedVector(const std::vector<double> &record,const std::vector<double> &fresh)
{
    if (record.size() != fresh.size()) fail("Native mixed restart independent port observation count differs.");
    for (std::size_t i=0;i<record.size();++i)
        if (!derivedEqual(record[i],fresh[i]))
            fail("Native mixed restart independent port observations differ from fresh assembly.");
}

class StateWriter
{
    QCryptographicHash hash{QCryptographicHash::Sha256};
public:
    void integer(std::uint64_t value)
    {
        char bytes[8];
        for (unsigned i=0;i<8;++i) bytes[i] = char((value >> (8*i)) & 0xffu);
        hash.addData(bytes,8);
    }
    void real(double value)
    {
        std::uint64_t bits = 0;
        static_assert(sizeof(value) == sizeof(bits),"Native mixed state fingerprint requires binary64 doubles.");
        std::memcpy(&bits,&value,sizeof(bits));
        integer(bits); // Preserve every bit, including signed zero.
    }
    void bytes(const QByteArray &value)
    {
        integer(std::uint64_t(value.size())); hash.addData(value);
    }
    void text(const QString &value) { bytes(value.toUtf8()); }
    void reals(const std::vector<double> &values)
    {
        integer(values.size()); for (double value : values) real(value);
    }
    void integers(const std::vector<uint> &values)
    {
        integer(values.size()); for (uint value : values) integer(value);
    }
    void texts(const std::vector<QString> &values)
    {
        integer(values.size()); for (const auto &value : values) text(value);
    }
    QString digest() { return QString::fromLatin1(hash.result().toHex()); }
};
}

RStructuralTrussDynamicsRestartData RStructuralTrussDynamicsRestartAdapter::build(
    const RStructuralTrussDynamicsModelAdapter &adapter,const RStructuralTrussDynamicsState &state,
    uint completed,uint cursor)
{
    verifyPhase(adapter,state,completed,cursor);
    const auto assembly = observe(adapter,state);
    RStructuralTrussDynamicsRestartData record;
    record.available = true;
    record.completedSteps = completed; record.runTargetSteps = adapter.targetSteps();
    record.initialTime = adapter.initialTime(); record.time = state.time;
    record.fixedTimeStep = adapter.nominalTimeStep(); record.nextCommandIndex = cursor;
    record.orderedPhysicalNodeIds = adapter.maps().orderedPhysicalNodeIds;
    record.orderedMemberPhysicalIds = adapter.maps().orderedMemberPhysicalIds;
    record.orderedPortIds = adapter.maps().orderedPortIds;
    record.sourceFingerprint = adapter.sourceFingerprint(); record.systemIdentity = state.systemIdentity;
    record.q = state.displacements; record.v = state.velocities; record.a = state.accelerations;
    record.reactions = state.reactions; record.externalForce = state.appliedInputs.externalForces;
    for (const auto &input : state.appliedInputs.memberInputs)
    {
        record.memberInputKinds.push_back(uint(input.index())); record.memberInputValues.push_back(inputValue(input));
    }
    record.storedEnergy = state.storedEnergy; record.kineticEnergy = state.kineticEnergy;
    record.initialTotalEnergy = state.initialTotalEnergy;
    record.cumulativeExternalWork = state.cumulativeExternalWork;
    record.cumulativeSupportWork = state.cumulativeSupportWork;
    record.cumulativeElectricalWork = state.cumulativeElectricalWork;
    record.cumulativeRestLengthWork = state.cumulativeLengthWork;
    portObservations(adapter,assembly,record.portVoltages,record.portCharges,record.portEnthalpies);
    record.validate();
    return record;
}

RStructuralTrussDynamicsRestoredState RStructuralTrussDynamicsRestartAdapter::restore(
    const RStructuralTrussDynamicsModelAdapter &adapter,const RStructuralTrussDynamicsRestartData &record)
{
    record.validate();
    if (!record.available || record.sourceFingerprint != adapter.sourceFingerprint() ||
        record.systemIdentity != adapter.driver().systemIdentity() ||
        record.orderedPhysicalNodeIds != adapter.maps().orderedPhysicalNodeIds ||
        record.orderedMemberPhysicalIds != adapter.maps().orderedMemberPhysicalIds ||
        record.orderedPortIds != adapter.maps().orderedPortIds ||
        record.initialTime != adapter.initialTime() || record.fixedTimeStep != adapter.nominalTimeStep() ||
        record.runTargetSteps != adapter.targetSteps())
        fail("Native mixed restart does not match frozen source/system, ordered physical maps or original run time grid.");
    RStructuralTrussDynamicsRestoredState restored;
    auto &state = restored.state;
    state.time = record.time; state.systemIdentity = record.systemIdentity;
    state.displacements = record.q; state.velocities = record.v; state.accelerations = record.a;
    state.reactions = record.reactions; state.appliedInputs.externalForces = record.externalForce;
    for (std::size_t i=0;i<record.memberInputKinds.size();++i)
    {
        const auto value = record.memberInputValues[i];
        switch (record.memberInputKinds[i])
        {
            case 0: state.appliedInputs.memberInputs.emplace_back(RStructuralTrussNoInput{}); break;
            case 1: state.appliedInputs.memberInputs.emplace_back(RStructuralTrussVoltage{value}); break;
            case 2: state.appliedInputs.memberInputs.emplace_back(RStructuralTrussCharge{value}); break;
            case 3: state.appliedInputs.memberInputs.emplace_back(RStructuralTrussCableLength{value}); break;
            default: fail("Native mixed restart member input kind is unsupported.");
        }
    }
    state.storedEnergy = record.storedEnergy; state.kineticEnergy = record.kineticEnergy;
    state.initialTotalEnergy = record.initialTotalEnergy;
    state.cumulativeExternalWork = record.cumulativeExternalWork;
    state.cumulativeSupportWork = record.cumulativeSupportWork;
    state.cumulativeElectricalWork = record.cumulativeElectricalWork;
    state.cumulativeLengthWork = record.cumulativeRestLengthWork;
    verifyPhase(adapter,state,record.completedSteps,record.nextCommandIndex);
    const auto assembly = observe(adapter,state);
    std::vector<double> voltages,charges,enthalpies;
    portObservations(adapter,assembly,voltages,charges,enthalpies);
    verifyDerivedVector(record.portVoltages,voltages);
    verifyDerivedVector(record.portCharges,charges);
    verifyDerivedVector(record.portEnthalpies,enthalpies);
    restored.completedSteps = record.completedSteps; restored.nextCommandIndex = record.nextCommandIndex;
    return restored;
}

QString RStructuralTrussDynamicsRestartAdapter::stateFingerprint(const RStructuralTrussDynamicsRestartData &record)
{
    record.validate();
    StateWriter writer;
    writer.text(QStringLiteral("range-fea/native-mixed-truss/accepted-state/v1/binary64"));
    writer.integer(record.available); writer.integer(record.completedSteps); writer.integer(record.runTargetSteps);
    writer.real(record.initialTime); writer.real(record.time); writer.real(record.fixedTimeStep);
    writer.integer(record.nextCommandIndex);
    writer.integers(record.orderedPhysicalNodeIds); writer.integers(record.orderedMemberPhysicalIds);
    writer.texts(record.orderedPortIds);
    writer.reals(record.portVoltages); writer.reals(record.portCharges); writer.reals(record.portEnthalpies);
    writer.text(record.sourceFingerprint); writer.bytes(record.systemIdentity);
    writer.reals(record.q); writer.reals(record.v); writer.reals(record.a); writer.reals(record.reactions);
    writer.reals(record.externalForce); writer.integers(record.memberInputKinds); writer.reals(record.memberInputValues);
    writer.real(record.storedEnergy); writer.real(record.kineticEnergy); writer.real(record.initialTotalEnergy);
    writer.real(record.cumulativeExternalWork); writer.real(record.cumulativeSupportWork);
    writer.real(record.cumulativeElectricalWork); writer.real(record.cumulativeRestLengthWork);
    return writer.digest();
}
