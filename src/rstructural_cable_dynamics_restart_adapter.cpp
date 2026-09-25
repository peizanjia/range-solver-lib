#include "rstructural_cable_dynamics_restart_adapter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>

#include <QCryptographicHash>

#include <rbl_error.h>

#include <rml_dof.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

namespace
{

[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}

void appendUInt64(QByteArray &bytes, std::uint64_t value)
{
    for (unsigned int shift=0;shift<64;shift+=8)
    {
        bytes.append(char((value >> shift) & 0xffu));
    }
}

void appendDouble(QByteArray &bytes, double value)
{
    if (value == 0.0) value = 0.0;
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),"Cable restart fingerprint requires 64-bit doubles.");
    std::memcpy(&bits,&value,sizeof(bits));
    appendUInt64(bytes,bits);
}

void appendVector(QByteArray &bytes, const std::vector<double> &values)
{
    appendUInt64(bytes,std::uint64_t(values.size()));
    for (double value : values) appendDouble(bytes,value);
}

std::vector<RStructuralDofRef> translationalDofs(uint nodeCount)
{
    std::vector<RStructuralDofRef> dofs;
    dofs.reserve(std::size_t(3)*nodeCount);
    for (uint node=0;node<nodeCount;node++)
    {
        dofs.push_back({node,RDofType::UX});
        dofs.push_back({node,RDofType::UY});
        dofs.push_back({node,RDofType::UZ});
    }
    return dofs;
}

double vectorScale(const std::vector<double> &values)
{
    double scale = 0.0;
    for (double value : values) scale = std::max(scale,std::fabs(value));
    return std::max(scale,1.0e-14);
}

void requireCloseVectors(const std::vector<double> &stored,
                         const std::vector<double> &checkpoint,
                         const char *description)
{
    if (stored.size() != checkpoint.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Cable transient restart %s length does not match its checkpoint.",description);
    }
    const double scale = std::max(vectorScale(stored),vectorScale(checkpoint));
    const double tolerance = 2.0e-6*scale +
        64.0*std::numeric_limits<double>::epsilon()*scale;
    for (std::size_t i=0;i<stored.size();i++)
    {
        if (!std::isfinite(stored[i]) || !std::isfinite(checkpoint[i]) ||
            std::fabs(stored[i]-checkpoint[i]) > tolerance)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable transient restart %s differs from its accepted checkpoint.",description);
        }
    }
}

void validateCheckpoint(const RModel &model,
                        const QString &systemFingerprint,
                        const std::vector<double> &externalForces)
{
    const RStructuralDynamicsRestartData &checkpoint =
        model.getProblemSetup().getStructuralDynamicsRestartData();
    checkpoint.validate();
    if (!checkpoint.available)
    {
        fail("Cable transient restart has no accepted structural checkpoint.");
    }
    const std::vector<RStructuralDofRef> dofs = translationalDofs(model.getNNodes());
    if (!checkpoint.matches(model,dofs) ||
        checkpoint.systemFingerprint != systemFingerprint ||
        checkpoint.q.size() != std::size_t(3)*model.getNNodes() ||
        checkpoint.externalForce.size() != externalForces.size())
    {
        fail("Cable transient restart checkpoint does not match its mesh, DOF order or cable system fingerprint.");
    }
    requireCloseVectors(checkpoint.externalForce,externalForces,"external force");

    const double configuredDt =
        model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    if (!std::isfinite(configuredDt) || configuredDt <= 0.0 ||
        std::fabs(checkpoint.fixedTimeStep-configuredDt) >
            1.0e-12*std::max(checkpoint.fixedTimeStep,configuredDt))
    {
        fail("Cable transient restart fixed time step differs from its accepted checkpoint.");
    }
}

void validateAcceptedIndexAndTime(const RModel &model,
                                  const RStructuralDynamicsRestartData &checkpoint)
{
    const RTimeSolver &time = model.getTimeSolver();
    if (checkpoint.completedSteps != time.getCurrentTimeStep())
    {
        fail("Cable transient trial does not follow the last accepted checkpoint step.");
    }
    const double previousTime = time.getPreviousTime();
    if (!std::isfinite(previousTime) ||
        std::fabs(checkpoint.time-previousTime) >
            1.0e-10*std::max({1.0,std::fabs(checkpoint.time),std::fabs(previousTime)}))
    {
        fail("Cable transient trial time does not follow the last accepted checkpoint time.");
    }
}

}

QString RStructuralCableDynamicsRestartAdapter::computeSystemFingerprint(
    const RModel &model,
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<double> &elementAreas,
    const std::vector<double> &nodalLumpedMass,
    const std::vector<bool> &prescribedMask,
    const std::vector<double> &prescribedDisplacements,
    const std::vector<double> &externalForces,
    const RStructuralCableDynamicsOptions &options)
{
    const std::size_t expectedDofs = std::size_t(3)*model.getNNodes();
    if (referenceNodes.size() != model.getNNodes() ||
        members.size() != elementAreas.size() ||
        nodalLumpedMass.size() != model.getNNodes() ||
        prescribedMask.size() != expectedDofs ||
        prescribedDisplacements.size() != expectedDofs ||
        externalForces.size() != expectedDofs)
    {
        fail("Cable transient fingerprint inputs do not match the current model dimensions.");
    }

    QByteArray bytes("range-fea/cable-dynamics/restart/v1");
    appendUInt64(bytes,model.getNNodes());
    appendUInt64(bytes,model.getNElements());
    appendUInt64(bytes,std::uint64_t(members.size()));
    for (const RStructuralCableVector3 &node : referenceNodes)
    {
        for (double coordinate : node) appendDouble(bytes,coordinate);
    }
    for (std::size_t index=0;index<members.size();index++)
    {
        const RStructuralCableEquilibriumMember &member = members[index];
        appendUInt64(bytes,std::uint64_t(member.node1));
        appendUInt64(bytes,std::uint64_t(member.node2));
        appendDouble(bytes,member.referenceLength);
        appendDouble(bytes,member.initialTension);
        appendDouble(bytes,member.axialRigidity);
        appendDouble(bytes,elementAreas[index]);
    }
    appendVector(bytes,nodalLumpedMass);
    appendUInt64(bytes,std::uint64_t(prescribedMask.size()));
    for (std::size_t dof=0;dof<prescribedMask.size();dof++)
    {
        bytes.append(prescribedMask[dof] ? char(1) : char(0));
        appendDouble(bytes,prescribedDisplacements[dof]);
    }
    appendVector(bytes,externalForces);
    appendDouble(bytes,options.dt);
    appendDouble(bytes,options.beta);
    appendDouble(bytes,options.gamma);
    appendUInt64(bytes,std::uint64_t(options.maxIterations));
    appendUInt64(bytes,std::uint64_t(options.maxLineSearchIterations));
    appendDouble(bytes,options.absForceTol);
    appendDouble(bytes,options.relForceTol);

    return QString::fromLatin1(QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex());
}

void RStructuralCableDynamicsRestartAdapter::validateAcceptedState(
    const RModel &model,
    const QString &systemFingerprint,
    const std::vector<double> &externalForces)
{
    if (!model.getProblemSetup().getRestart()) return;
    validateCheckpoint(model,systemFingerprint,externalForces);
}

void RStructuralCableDynamicsRestartAdapter::restoreAcceptedState(
    const RModel &model,
    const QString &systemFingerprint,
    const std::vector<double> &externalForces,
    std::vector<double> &displacements,
    std::vector<double> &velocities,
    std::vector<double> &accelerations)
{
    const RStructuralDynamicsRestartData &checkpoint =
        model.getProblemSetup().getStructuralDynamicsRestartData();
    if (!checkpoint.available)
    {
        if (model.getProblemSetup().getRestart())
        {
            fail("Cable transient restart has no accepted structural checkpoint.");
        }
        return;
    }

    validateCheckpoint(model,systemFingerprint,externalForces);
    validateAcceptedIndexAndTime(model,checkpoint);
    requireCloseVectors(displacements,checkpoint.q,"displacement");
    requireCloseVectors(velocities,checkpoint.v,"velocity");
    requireCloseVectors(accelerations,checkpoint.a,"acceleration");

    // ASCII model results may be rounded for display. The checkpoint stores
    // 17-digit accepted state and is the authoritative continuation state.
    displacements = checkpoint.q;
    velocities = checkpoint.v;
    accelerations = checkpoint.a;
}

RStructuralDynamicsRestartData RStructuralCableDynamicsRestartAdapter::stageAcceptedState(
    const RModel &model,
    const RStructuralCableDynamicsResult &result,
    const QString &systemFingerprint,
    const std::vector<double> &externalForces)
{
    const std::size_t expectedDofs = std::size_t(3)*model.getNNodes();
    if (result.displacements.size() != expectedDofs ||
        result.velocities.size() != expectedDofs ||
        result.accelerations.size() != expectedDofs ||
        externalForces.size() != expectedDofs)
    {
        fail("Cable transient accepted checkpoint vectors do not match the model DOF count.");
    }

    const RStructuralDynamicsRestartData &previous =
        model.getProblemSetup().getStructuralDynamicsRestartData();
    if (model.getProblemSetup().getRestart() && !previous.available)
    {
        fail("Cable transient restart has no accepted structural checkpoint.");
    }
    RStructuralDynamicsRestartData candidate;
    candidate.meshNodeCount = model.getNNodes();
    candidate.meshElementCount = model.getNElements();
    candidate.physicalDofCount = uint(expectedDofs);
    candidate.meshFingerprint = RStructuralDynamicsRestartData::computeMeshFingerprint(model);
    candidate.dofFingerprint = RStructuralDynamicsRestartData::computeDofFingerprint(
        translationalDofs(model.getNNodes()));
    candidate.systemFingerprint = systemFingerprint;
    candidate.completedSteps = model.getTimeSolver().getCurrentTimeStep()+1;
    candidate.time = model.getTimeSolver().getCurrentTime();
    candidate.fixedTimeStep =
        model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    candidate.q = result.displacements;
    candidate.v = result.velocities;
    candidate.a = result.accelerations;
    candidate.externalForce = externalForces;

    if (previous.available)
    {
        validateCheckpoint(model,systemFingerprint,externalForces);
        validateAcceptedIndexAndTime(model,previous);
        candidate.cumulativeExternalWork = previous.cumulativeExternalWork+
            result.externalWorkIncrement+result.supportWorkIncrement;
        candidate.cumulativeDampingDissipation = previous.cumulativeDampingDissipation;
        candidate.initialMechanicalEnergy = previous.initialMechanicalEnergy;
    }
    else
    {
        candidate.cumulativeExternalWork = result.externalWorkIncrement+
            result.supportWorkIncrement;
        candidate.cumulativeDampingDissipation = 0.0;
        candidate.initialMechanicalEnergy = result.previousStrainEnergy+
            result.previousKineticEnergy;
    }

    candidate.available = true;
    candidate.validate();
    return candidate;
}

void RStructuralCableDynamicsRestartAdapter::commitAcceptedState(
    RModel &model,
    RStructuralDynamicsRestartData candidate)
{
    candidate.validate();
    model.getProblemSetup().getStructuralDynamicsRestartData() = std::move(candidate);
}
