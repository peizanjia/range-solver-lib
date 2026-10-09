#ifndef RSTRUCTURAL_CONTROL_RUNTIME_PARAMETERS_H
#define RSTRUCTURAL_CONTROL_RUNTIME_PARAMETERS_H

#include <QByteArray>
#include <QMap>
#include <QString>
#include <QtGlobal>

#include <cstddef>
#include <vector>

//! Pure preparation and checkpoint codec for accepted WF09 runtime parameters.
//!
//! The codec stores only a current parameter set and its accepted change
//! history. It does not represent proposals, approvals, solver rollback state,
//! or an RCCE composite envelope.
namespace RStructuralControlRuntimeParameters
{
struct ActuatorLimit
{
    QString actuatorId;
    double minN = 0.0;
    double maxN = 0.0;
    double rateNPerS = 0.0; // Zero means unlimited.
};

//! parameterSha256 uses a portable binary encoding: the ASCII marker
//! "range.structural-runtime-parameters/hash/v1" plus NUL, 32 raw base-SHA
//! bytes, u64 revision BE, f64 gain BE, u32 actuator count BE, then each sorted
//! actuator's u32 UTF-8 ID length BE, ID bytes, and f64 minN/maxN/rateNPerS BE.
//! Signed zero is encoded as positive zero. runtimeControllerRevision is
//! SHA-256 over ASCII base-revision and parameter-SHA strings, each prefixed
//! by its four-byte big-endian byte length.
struct Parameters
{
    quint64 revision = 0;
    double velocityGainNsm = 1.0;
    std::vector<ActuatorLimit> actuatorLimits;
    QString fixedBaseControllerRevision;
    QString parameterSha256;
    QString runtimeControllerRevision;
};

struct Change
{
    QString proposalId;
    quint64 previousRevision = 0;
    quint64 revision = 0;
    quint64 appliedSeq = 0;
    double acceptedTimeS = 0.0;
    QString parameterSha256;
};

struct State
{
    Parameters parameters;
    std::vector<Change> history;
};

//! Feasible scalar-force bounds after intersecting amplitude and rate limits.
struct AppliedForceInterval
{
    QString actuatorId;
    double minN = 0.0;
    double maxN = 0.0;
};

constexpr std::size_t maximumActuatorCount = 256;
constexpr std::size_t maximumHistoryCount = 128;
constexpr qsizetype maximumStateBytes = 64 * 1024;

//! Construct revision zero. Limits are validated and sorted by actuator ID.
State fresh(const QString &fixedBaseControllerRevision,
            double velocityGainNsm,
            std::vector<ActuatorLimit> actuatorLimits);

//! Prepare revision current+1 with the same base revision and actuator IDs.
//! This creates a value object only; it does not add it to accepted history.
Parameters candidate(const State &accepted,
                    double velocityGainNsm,
                    std::vector<ActuatorLimit> actuatorLimits);

//! Validate a prepared candidate against immutable base limits and the latest
//! applied force. Returns the non-empty per-actuator force interval for this
//! step. Call this before acceptedChange().
std::vector<AppliedForceInterval> validateCandidateAgainstBase(
    const Parameters &candidate,
    const std::vector<ActuatorLimit> &immutableLimits,
    const QMap<QString, double> &lastAppliedById,
    double dt);

//! Return a new accepted state. The input state is never mutated. The caller
//! must first validate the candidate against its immutable limits and latest
//! applied forces. Proposal IDs are canonical lowercase UUID strings.
State acceptedChange(const State &accepted,
                     const QString &proposalId,
                     const Parameters &candidate,
                     quint64 appliedSeq,
                     double acceptedTimeS);

//! Encode the exact range.structural-runtime-parameters/v1 canonical JSON
//! byte form. Decoding accepts only bytes identical to compact Qt JSON emitted
//! by this codec; duplicate keys, whitespace variants and trailing data fail.
//! Parsing also recomputes both SHA-256 values from the decoded state.
//! The encoded document is bounded to 65536 bytes.
QByteArray serialize(const State &state);
State deserialize(const QByteArray &bytes);
}

#endif // RSTRUCTURAL_CONTROL_RUNTIME_PARAMETERS_H
