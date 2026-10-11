#include "rstructural_truss_dynamics.h"
#include "rmatrixsolver.h"

#include <QCryptographicHash>
#include <QDataStream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace
{
double finite(double x)
{
    if (!std::isfinite(x)) throw std::invalid_argument("Mixed truss dynamics non-finite or overflowing value");
    return x;
}
void require(bool condition, const char *message)
{
    if (!condition) throw std::invalid_argument(message);
}
double norm(const std::vector<double> &v)
{
    double n = 0.0;
    for (double x : v) n = finite(std::hypot(n,finite(x)));
    return n;
}
bool close(double a, double b, double relative = 128.0*std::numeric_limits<double>::epsilon())
{
    finite(a); finite(b);
    return std::fabs(finite(a-b)) <= finite(relative*std::max({1.0,std::fabs(a),std::fabs(b)}));
}
std::size_t motionDof(const RStructuralPrescribedMotion &motion)
{
    require(motion.target.dof == RDofType::UX || motion.target.dof == RDofType::UY ||
            motion.target.dof == RDofType::UZ,"Mixed truss motion requires UX/UY/UZ");
    const std::size_t axis = motion.target.dof == RDofType::UX ? 0 : motion.target.dof == RDofType::UY ? 1 : 2;
    return std::size_t(motion.target.nodeId)*3+axis;
}
void optionsValid(const RStructuralTrussDynamicsOptions &o)
{
    require(o.maxIterations > 0 && o.maxLineSearchIterations > 0 &&
            finite(o.absForceTol) >= 0.0 && finite(o.relForceTol) >= 0.0 &&
            o.absForceTol <= 1.0e-8 && o.relForceTol <= 1.0e-8 &&
            (o.absForceTol > 0.0 || o.relForceTol > 0.0) &&
            finite(o.linearRelativeTolerance) > 0.0 && o.linearRelativeTolerance < 1.0,
            "Invalid mixed truss dynamics convergence options");
}
double diagnostic(const RStructuralTrussDynamicsState &s)
{
    double error = finite(finite(s.storedEnergy+s.kineticEnergy)-s.initialTotalEnergy);
    for (double w : {s.cumulativeExternalWork,s.cumulativeSupportWork,
                     s.cumulativeElectricalWork,s.cumulativeLengthWork}) error = finite(error-w);
    return error;
}
void accumulate(RStructuralTrussDynamicsState &s, const RStructuralTrussDynamicsWork &w)
{
    s.cumulativeExternalWork = finite(s.cumulativeExternalWork+w.externalWork);
    s.cumulativeSupportWork = finite(s.cumulativeSupportWork+w.supportWork);
    s.cumulativeElectricalWork = finite(s.cumulativeElectricalWork+w.electricalWork);
    s.cumulativeLengthWork = finite(s.cumulativeLengthWork+w.lengthWork);
}
struct Evaluation
{
    RStructuralTrussAssemblyResult assembly;
    std::vector<double> residual;
    double residualNorm = 0.0, scale = 0.0;
};
Evaluation evaluate(const RStructuralTrussDynamicsSystem &system,
                    const RStructuralTrussDynamicsState &state,
                    const std::vector<std::size_t> &freeDofs)
{
    Evaluation e;
    e.assembly = RStructuralTrussAssembly::evaluate(system.referenceNodes,state.displacements,
                                                   system.members,state.appliedInputs.memberInputs);
    double internal = 0.0, external = 0.0, inertia = 0.0;
    for (std::size_t d : freeDofs)
    {
        const double ma = finite(system.nodalLumpedMass[d/3]*state.accelerations[d]);
        const double r = finite(finite(ma+e.assembly.internalForces[d])-state.appliedInputs.externalForces[d]);
        e.residual.push_back(r);
        internal = finite(std::hypot(internal,e.assembly.internalForces[d]));
        external = finite(std::hypot(external,state.appliedInputs.externalForces[d]));
        inertia = finite(std::hypot(inertia,ma));
    }
    e.residualNorm = norm(e.residual);
    e.scale = std::max({internal,external,inertia});
    return e;
}
std::vector<double> linearSolve(const RSparseMatrix &matrix, const std::vector<double> &rhs,
                                double relativeTolerance)
{
    const uint n = uint(rhs.size());
    double matrixScale = 0.0;
    for (uint i=0;i<n;++i)
        for (uint j=0;j<matrix.getNColumns(i);++j)
            matrixScale = std::max(matrixScale,std::fabs(finite(matrix.getValue(i,j))));
    require(matrixScale > 0.0,"Mixed truss effective tangent is zero");
    // Keep the legacy GMRES squared-norm implementation away from overflowing
    // matrix magnitudes; certification below uses the original matrix/rhs.
    RSparseMatrix scaled;
    scaled.setNRows(n);
    RRVector b(n,0.0), x(n,0.0);
    double matrixNorm = 0.0;
    for (uint i=0;i<n;++i)
    {
        b[i] = finite(rhs[i]/matrixScale);
        const auto columns = matrix.getRowIndexes(i);
        for (uint j=0;j<columns.size();++j)
        {
            const double entry = finite(matrix.getValue(i,j)/matrixScale);
            scaled.addValue(i,columns[j],entry);
            matrixNorm = finite(std::hypot(matrixNorm,entry));
        }
    }
    const double scaledNorm = norm(std::vector<double>(b.begin(),b.end()));
    require(scaledNorm > 0.0,"Mixed truss linear rhs lost representability");
    require(finite(scaledNorm*scaledNorm) > 0.0,"Mixed truss native linear norm underflow");
    const double equationScale = finite(1.0e9/std::fabs(scaledNorm/matrixNorm));
    require(equationScale > 0.0,"Mixed truss native equation scale underflow");
    finite(finite(scaledNorm*equationScale)*finite(scaledNorm*equationScale));
    RMatrixSolverConf conf(RMatrixSolverConf::GMRES);
    conf.setNInnerIterations(std::min<uint>(n,100));
    conf.setNOuterIterations(std::max<uint>(100,n));
    require(relativeTolerance*0.1 > 0.0,"Mixed truss linear tolerance underflow");
    conf.setSolverCvgValue(relativeTolerance*0.1);
    conf.setOutputFrequency(0);
    RMatrixSolver solver(conf);
    solver.disableConvergenceLogFile();
    solver.solve(scaled,b,x,R_MATRIX_PRECONDITIONER_NONE);
    std::vector<double> result(x.begin(),x.end()), residual(n,0.0);
    for (uint i=0;i<n;++i)
    {
        finite(result[i]);
        double ax = 0.0;
        const auto columns = matrix.getRowIndexes(i);
        for (uint j=0;j<columns.size();++j)
            ax = finite(ax+finite(matrix.getValue(i,j)*result[columns[j]]));
        residual[i] = finite(ax-rhs[i]);
    }
    if (norm(residual) > finite(relativeTolerance*norm(rhs)))
        throw std::runtime_error("Mixed truss GMRES failed fresh linear residual certification");
    return result;
}
}

RStructuralTrussDynamics::RStructuralTrussDynamics(const RStructuralTrussDynamicsSystem &input)
    : system(input)
{
    require(!system.referenceNodes.empty() && system.referenceNodes.size() <= std::numeric_limits<uint>::max()/3,
            "Mixed truss dynamics requires a nonempty representable node array");
    const std::size_t n = 3*system.referenceNodes.size();
    require(system.nodalLumpedMass.size() == system.referenceNodes.size() &&
            system.prescribedMask.size() == n && system.fixedPrescribedDisplacements.size() == n,
            "Mixed truss system dimensions are invalid");
    for (const auto &node : system.referenceNodes) for (double x : node) finite(x);
    for (double m : system.nodalLumpedMass) require(finite(m) >= 0.0,"Negative mixed truss mass");
    std::set<QString> motionIds;
    std::set<std::size_t> motionDofs;
    for (const auto &motion : system.prescribedMotions)
    {
        motion.validate();
        const auto d = motionDof(motion);
        require(d < n && system.prescribedMask[d] && motionIds.insert(motion.id).second &&
                motionDofs.insert(d).second,"Invalid or duplicate mixed truss prescribed motion ownership");
    }
    for (std::size_t d=0;d<n;++d)
    {
        finite(system.fixedPrescribedDisplacements[d]);
        if (!system.prescribedMask[d])
        {
            require(system.nodalLumpedMass[d/3] > 0.0,"Every free mixed truss DOF requires positive mass");
            freeDofs.push_back(d);
        }
        if (!system.prescribedMask[d] || motionDofs.count(d))
            require(system.fixedPrescribedDisplacements[d] == 0.0,"Ambiguous fixed prescribed displacement ownership");
    }
    // Validate frozen data alone. The initial physical geometry and real
    // electrical/length boundary are evaluated only by initialize.
    std::set<std::uint64_t> memberIds;
    for (const auto &member : system.members)
    {
        require(member.physicalElementId != 0 && memberIds.insert(member.physicalElementId).second &&
                member.node1 < system.referenceNodes.size() && member.node2 < system.referenceNodes.size() &&
                member.node1 != member.node2,"Invalid mixed truss frozen member identity/connectivity");
        RStructuralTrussVector3 chord{};
        for (std::size_t axis=0;axis<3;++axis)
            chord[axis] = finite(system.referenceNodes[member.node2][axis]-system.referenceNodes[member.node1][axis]);
        require(finite(std::hypot(chord[0],chord[1],chord[2])) > 0.0,"Zero mixed truss reference chord");
        if (const auto *passive = std::get_if<RStructuralTrussPassiveLaw>(&member.law))
        {
            require(finite(passive->referenceLength) > 0.0 && finite(passive->axialRigidity) > 0.0,
                    "Invalid frozen passive rod law");
            require(finite(passive->axialRigidity/passive->referenceLength) > 0.0,
                    "Unrepresentable frozen passive rod stiffness");
        }
        else if (const auto *piezo = std::get_if<RStructuralTrussPiezoLaw>(&member.law))
        {
            require(finite(piezo->referenceLength) > 0.0 && finite(piezo->axialStiffness) > 0.0 &&
                    finite(piezo->capacitance) > 0.0,"Invalid frozen piezo rod law");
            finite(piezo->coupling);
        }
        else
        {
            const auto &cable = std::get<RStructuralTrussCableLaw>(member.law);
            require(finite(cable.initialTension) >= 0.0 && finite(cable.axialRigidity) > 0.0,
                    "Invalid frozen cable law");
        }
    }
    QByteArray bytes;
    QDataStream stream(&bytes,QIODevice::WriteOnly);
    stream.setVersion(QDataStream::Qt_5_0);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    stream << QStringLiteral("mixed-truss-dynamics-v1") << quint64(system.referenceNodes.size());
    for (const auto &node : system.referenceNodes) for (double x : node) stream << x;
    stream << quint64(system.members.size());
    for (const auto &m : system.members)
    {
        stream << quint64(m.physicalElementId) << quint64(m.node1) << quint64(m.node2) << quint64(m.law.index());
        if (const auto *law = std::get_if<RStructuralTrussPassiveLaw>(&m.law)) stream << law->referenceLength << law->axialRigidity;
        else if (const auto *law = std::get_if<RStructuralTrussPiezoLaw>(&m.law))
            stream << law->referenceLength << law->axialStiffness << law->coupling << law->capacitance;
        else { const auto &cableLaw = std::get<RStructuralTrussCableLaw>(m.law); stream << cableLaw.initialTension << cableLaw.axialRigidity; }
    }
    for (double m : system.nodalLumpedMass) stream << m;
    for (std::size_t d=0;d<n;++d) stream << bool(system.prescribedMask[d]) << system.fixedPrescribedDisplacements[d];
    stream << quint64(system.prescribedMotions.size());
    for (const auto &motion : system.prescribedMotions)
    {
        stream << motion.id << quint64(motion.target.nodeId) << quint64(motion.target.dof) << quint64(motion.knots.size());
        for (const auto &k : motion.knots) stream << k.time << k.displacement << k.velocity << k.acceleration;
    }
    require(stream.status() == QDataStream::Ok,"Cannot encode mixed truss system identity");
    identity = QCryptographicHash::hash(bytes,QCryptographicHash::Sha256);
}

void RStructuralTrussDynamics::prescribed(double time, std::vector<double> &q,
                                        std::vector<double> &v, std::vector<double> &a) const
{
    finite(time);
    for (std::size_t d=0;d<q.size();++d) if (system.prescribedMask[d])
    { q[d] = system.fixedPrescribedDisplacements[d]; v[d] = a[d] = 0.0; }
    for (const auto &motion : system.prescribedMotions)
    {
        // The shared analytic evaluator accepts bounded endpoint roundoff,
        // while rejecting genuine extrapolation. A stricter gate here would
        // reject an ordinary decimal-step schedule at its final knot.
        const auto value = motion.evaluate(time);
        const auto d = motionDof(motion);
        q[d] = finite(value.displacement); v[d] = finite(value.velocity); a[d] = finite(value.acceleration);
    }
}

void RStructuralTrussDynamics::validateInputs(const RStructuralTrussDynamicsInputs &inputs) const
{
    require(inputs.memberInputs.size() == system.members.size() &&
            inputs.externalForces.size() == system.prescribedMask.size(),"Mixed truss input dimensions are invalid");
    for (double f : inputs.externalForces) finite(f);
    // Assembly validates the law/input pairing at the actual evaluation geometry.
    for (std::size_t i=0;i<inputs.memberInputs.size();++i)
    {
        const auto &law = system.members[i].law;
        const auto &command = inputs.memberInputs[i];
        if (std::holds_alternative<RStructuralTrussPassiveLaw>(law))
            require(std::holds_alternative<RStructuralTrussNoInput>(command),"Passive rod input mismatch");
        else if (std::holds_alternative<RStructuralTrussPiezoLaw>(law))
        {
            require(std::holds_alternative<RStructuralTrussVoltage>(command) ||
                    std::holds_alternative<RStructuralTrussCharge>(command),"Piezo input mismatch");
            if (const auto *c = std::get_if<RStructuralTrussVoltage>(&command)) finite(c->value);
            else finite(std::get<RStructuralTrussCharge>(command).value);
        }
        else
        {
            require(std::holds_alternative<RStructuralTrussCableLength>(command),"Cable input mismatch");
            require(finite(std::get<RStructuralTrussCableLength>(command).value) > 0.0,"Invalid cable held length");
        }
    }
}

void RStructuralTrussDynamics::observe(RStructuralTrussDynamicsState &s,
                                     const RStructuralTrussAssemblyResult &assembly) const
{
    s.storedEnergy = assembly.storedEnergy;
    s.kineticEnergy = 0.0;
    s.reactions.resize(s.displacements.size());
    for (std::size_t d=0;d<s.displacements.size();++d)
    {
        s.reactions[d] = finite(finite(system.nodalLumpedMass[d/3]*s.accelerations[d])+
                               finite(assembly.internalForces[d]-s.appliedInputs.externalForces[d]));
        const double root = finite(std::sqrt(system.nodalLumpedMass[d/3])*s.velocities[d]);
        s.kineticEnergy = finite(s.kineticEnergy+finite((0.5*root)*root));
    }
}

RStructuralTrussAssemblyResult RStructuralTrussDynamics::validateState(const RStructuralTrussDynamicsState &s) const
{
    const auto n = system.prescribedMask.size();
    require(s.systemIdentity == identity && s.displacements.size() == n && s.velocities.size() == n &&
            s.accelerations.size() == n && s.reactions.size() == n,"Mixed truss accepted state identity or dimensions mismatch");
    finite(s.time);
    for (std::size_t d=0;d<n;++d) { finite(s.displacements[d]); finite(s.velocities[d]); finite(s.accelerations[d]); finite(s.reactions[d]); }
    validateInputs(s.appliedInputs);
    auto q = s.displacements, v = s.velocities, a = s.accelerations;
    prescribed(s.time,q,v,a);
    for (std::size_t d=0;d<n;++d) if (system.prescribedMask[d])
        require(close(q[d],s.displacements[d]) && close(v[d],s.velocities[d]) && close(a[d],s.accelerations[d]),
                "Accepted support kinematics disagree with frozen motion");
    const auto e = evaluate(system,s,freeDofs);
    require(e.residualNorm <= finite(1.0e-8+1.0e-8*e.scale),"Accepted acceleration is inconsistent with old applied inputs");
    auto observed = s;
    observe(observed,e.assembly);
    require(close(observed.storedEnergy,s.storedEnergy) && close(observed.kineticEnergy,s.kineticEnergy) &&
            finite(s.storedEnergy) >= 0.0 && finite(s.kineticEnergy) >= 0.0 && finite(s.initialTotalEnergy) >= 0.0,
            "Accepted energy observation is inconsistent");
    for (std::size_t d=0;d<n;++d) require(close(observed.reactions[d],s.reactions[d]),"Accepted reaction observation is inconsistent");
    diagnostic(s); // Certifies all cumulative work fields and their arithmetic.
    return e.assembly;
}

RStructuralTrussDynamicsState RStructuralTrussDynamics::initialize(double time,
    const std::vector<double> &q, const std::vector<double> &v, const RStructuralTrussDynamicsInputs &inputs) const
{
    const auto n = system.prescribedMask.size();
    require(q.size() == n && v.size() == n,"Mixed truss initial state dimensions are invalid");
    validateInputs(inputs);
    RStructuralTrussDynamicsState s;
    s.time = finite(time); s.displacements = q; s.velocities = v;
    s.accelerations.assign(n,0.0); s.appliedInputs = inputs; s.systemIdentity = identity;
    for (std::size_t d=0;d<n;++d) { finite(q[d]); finite(v[d]); }
    prescribed(time,s.displacements,s.velocities,s.accelerations);
    for (std::size_t d=0;d<n;++d) if (system.prescribedMask[d])
        require(close(q[d],s.displacements[d]) && close(v[d],s.velocities[d]),"Initial prescribed q/v mismatch");
    const auto assembly = RStructuralTrussAssembly::evaluate(system.referenceNodes,s.displacements,system.members,inputs.memberInputs);
    for (std::size_t d : freeDofs)
        s.accelerations[d] = finite(finite(inputs.externalForces[d]-assembly.internalForces[d])/system.nodalLumpedMass[d/3]);
    observe(s,assembly);
    s.initialTotalEnergy = finite(s.storedEnergy+s.kineticEnergy);
    validateState(s);
    return s;
}

RStructuralTrussDynamicsResult RStructuralTrussDynamics::refresh(const RStructuralTrussDynamicsState &accepted,
    const RStructuralTrussDynamicsInputs &inputs) const
{
    const auto old = validateState(accepted);
    validateInputs(inputs);
    for (std::size_t i=0;i<system.members.size();++i)
        require(accepted.appliedInputs.memberInputs[i].index() == inputs.memberInputs[i].index(),
                "Mixed truss voltage/charge boundary kind switching is unavailable");
    RStructuralTrussDynamicsResult result;
    auto &s = result.candidate;
    s = accepted; s.appliedInputs = inputs;
    const auto right = RStructuralTrussAssembly::evaluate(system.referenceNodes,s.displacements,system.members,inputs.memberInputs);
    for (std::size_t i=0;i<system.members.size();++i)
    {
        const double jump = finite(right.memberResponses[i].storedEnergy-old.memberResponses[i].storedEnergy);
        if (std::holds_alternative<RStructuralTrussPiezoLaw>(system.members[i].law))
            result.jumpWork.electricalWork = finite(result.jumpWork.electricalWork+jump);
        else if (std::holds_alternative<RStructuralTrussCableLaw>(system.members[i].law))
            result.jumpWork.lengthWork = finite(result.jumpWork.lengthWork+jump);
    }
    for (std::size_t d : freeDofs)
        s.accelerations[d] = finite(finite(inputs.externalForces[d]-right.internalForces[d])/system.nodalLumpedMass[d/3]);
    observe(s,right); accumulate(s,result.jumpWork);
    result.residualNorm = evaluate(system,s,freeDofs).residualNorm;
    result.energyBalanceError = diagnostic(s);
    validateState(s);
    return result;
}

RStructuralTrussDynamicsResult RStructuralTrussDynamics::trialStep(const RStructuralTrussDynamicsState &accepted,
    const RStructuralTrussDynamicsInputs &inputs, double dt, const RStructuralTrussDynamicsOptions &options) const
{
    optionsValid(options);
    require(finite(dt) > 0.0,"Mixed truss dt must be positive");
    const double endTime = finite(accepted.time+dt);
    require(endTime > accepted.time,"Mixed truss dt does not advance representable time");
    const double dt2 = finite(dt*dt);
    const double coefficient = finite(4.0/dt2);
    require(dt2 > 0.0 && coefficient > 0.0,"Mixed truss Newmark time scale is not representable");
    auto result = refresh(accepted,inputs);
    const auto start = result.candidate;
    const auto rightAssembly = RStructuralTrussAssembly::evaluate(system.referenceNodes,start.displacements,system.members,inputs.memberInputs);
    auto &s = result.candidate;
    s.time = endTime;
    const auto n = s.displacements.size();
    std::vector<double> predictor(n,0.0), velocityPredictor(n,0.0);
    for (std::size_t d : freeDofs)
    {
        predictor[d] = finite(finite(start.displacements[d]+finite(dt*start.velocities[d]))+
                              finite((0.25*dt2)*start.accelerations[d]));
        velocityPredictor[d] = finite(start.velocities[d]+finite((0.5*dt)*start.accelerations[d]));
        s.displacements[d] = predictor[d];
        s.accelerations[d] = 0.0;
    }
    prescribed(endTime,s.displacements,s.velocities,s.accelerations);
    auto kinematics = [&](RStructuralTrussDynamicsState &trial) {
        for (std::size_t d : freeDofs)
        {
            // Acceleration is the Newmark correction variable. Recovering it
            // by subtracting rounded O(q) displacements loses the small O(dt²)
            // correction and can make a valid force tolerance unattainable.
            trial.displacements[d] = finite(predictor[d]+finite((0.25*dt2)*trial.accelerations[d]));
            trial.velocities[d] = finite(velocityPredictor[d]+finite((0.5*dt)*trial.accelerations[d]));
        }
    };
    auto tolerance = [&](const Evaluation &e) { return finite(options.absForceTol+finite(options.relForceTol*e.scale)); };
    for (;;)
    {
        kinematics(s);
        const auto e = evaluate(system,s,freeDofs);
        if (e.residualNorm <= tolerance(e)) break;
        if (result.iterations >= options.maxIterations)
            throw std::runtime_error("Mixed truss Newton iteration limit exceeded");
        RSparseMatrix tangent;
        tangent.setNRows(uint(freeDofs.size()));
        std::vector<std::size_t> map(n,n);
        for (std::size_t i=0;i<freeDofs.size();++i) map[freeDofs[i]] = i;
        for (std::size_t i=0;i<freeDofs.size();++i)
        {
            const uint d = uint(freeDofs[i]);
            const auto columns = e.assembly.tangent.getRowIndexes(d);
            for (uint j=0;j<columns.size();++j) if (map[columns[j]] != n)
                tangent.addValue(uint(i),uint(map[columns[j]]),e.assembly.tangent.getValue(d,j));
            finite(tangent.findValue(uint(i),uint(i))+finite(coefficient*system.nodalLumpedMass[d/3]));
            tangent.addValue(uint(i),uint(i),coefficient*system.nodalLumpedMass[d/3]);
        }
        auto rhs = e.residual;
        for (double &r : rhs) r = -r;
        const auto increment = linearSolve(tangent,rhs,options.linearRelativeTolerance);
        bool reduced = false;
        double alpha = 1.0;
        for (std::size_t line=0;line<options.maxLineSearchIterations;++line)
        {
            auto trial = s;
            for (std::size_t i=0;i<freeDofs.size();++i)
                trial.accelerations[freeDofs[i]] = finite(s.accelerations[freeDofs[i]]+
                    finite(coefficient*finite(alpha*increment[i])));
            kinematics(trial);
            // Malformed/overflowing candidates propagate; they are never turned
            // into convergence or silently hidden by a catch-all line search.
            const auto trialEvaluation = evaluate(system,trial,freeDofs);
            if (trialEvaluation.residualNorm <= (1.0-1.0e-4*alpha)*e.residualNorm ||
                trialEvaluation.residualNorm <= tolerance(trialEvaluation))
            { s = std::move(trial); reduced = true; break; }
            alpha *= 0.5;
        }
        if (!reduced) throw std::runtime_error("Mixed truss Newton free residual merit line search failed");
        ++result.iterations;
    }
    kinematics(s);
    const auto fresh = evaluate(system,s,freeDofs);
    if (fresh.residualNorm > tolerance(fresh))
        throw std::runtime_error("Mixed truss failed fresh physical residual certification");
    observe(s,fresh.assembly);
    result.residualNorm = fresh.residualNorm;
    for (std::size_t d=0;d<n;++d)
    {
        result.stepWork.externalWork = finite(result.stepWork.externalWork+
            finite(inputs.externalForces[d]*finite(s.displacements[d]-start.displacements[d])));
        if (system.prescribedMask[d])
            result.stepWork.supportWork = finite(result.stepWork.supportWork+finite((0.5*dt)*
                finite(finite(start.reactions[d]*start.velocities[d])+finite(s.reactions[d]*s.velocities[d]))));
    }
    for (std::size_t i=0;i<system.members.size();++i)
        if (const auto *voltage = std::get_if<RStructuralTrussVoltage>(&inputs.memberInputs[i]))
            result.stepWork.electricalWork = finite(result.stepWork.electricalWork+finite(voltage->value*
                finite(fresh.assembly.memberResponses[i].electrical->conductorCharge-
                       rightAssembly.memberResponses[i].electrical->conductorCharge)));
    accumulate(s,result.stepWork);
    result.energyBalanceError = diagnostic(s);
    validateState(s);
    return result;
}
