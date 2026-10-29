#include "rstructural_truss_dynamics_model_adapter.h"
#include "rstructural_truss_dynamics_restart_adapter.h"
#include "rstructural_static_truss_model.h"
#include "rstructural_native_control.h"
#include "rstructural_piezo_native_readout.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <rbl_error.h>
#include <rml_model.h>
#include <rml_line.h>
#include <rml_variable.h>

namespace
{
[[noreturn]] void fail(const QString &message)
{
    const auto bytes = message.toUtf8();
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",bytes.constData());
}

bool close(double a,double b)
{
    return std::isfinite(a) && std::isfinite(b) &&
        std::fabs(a-b) <= 128.0*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::fabs(a),std::fabs(b)});
}

bool analyticMotionEqual(double a,double b)
{
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    if (a == b) return true;
    // Only interpolation roundoff at the actual kinematic magnitude is
    // allowed. No one-metre / one-metre-per-second absolute floor exists.
    const double magnitude = std::max(std::fabs(a),std::fabs(b));
    const double ulp = magnitude-std::nextafter(magnitude,0.0);
    return std::fabs(a-b) <= 32.0*ulp;
}

double positive(double value,const QString &context)
{
    if (!std::isfinite(value) || value <= 0.0) fail(context);
    return value;
}

double positiveProduct(double a,double b,double c,const QString &context)
{
    positive(a,context); positive(b,context); positive(c,context);
    int ea=0,eb=0,ec=0;
    const double mantissa = std::frexp(a,&ea)*std::frexp(b,&eb)*std::frexp(c,&ec);
    return positive(std::ldexp(mantissa,ea+eb+ec),context);
}

void preflight(const RModel &model)
{
    const auto &problem = model.getProblemSetup();
    const auto &dynamics = problem.getStructuralDynamicsSetup();
    const auto &mixed = problem.getStructuralTrussDynamicsSetup();
    mixed.validate(model); dynamics.validate();
    const auto &nativeControl = problem.getStructuralNativeControlSpec();
    nativeControl.validate(model);
    if (nativeControl.source == RStructuralNativeControlSource::LocalV3 && !mixed.heldCommands.empty())
        fail("Native LocalV3 cannot share action ownership with frozen held commands.");
    if (!mixed.enabled || model.getProblemTaskTree().getProblemTypeMask() != R_PROBLEM_STRESS)
        fail("Native mixed dynamics requires enabled native setup and a STRESS-only task.");
    if (problem.getMassMatrixType() != RMassMatrixType::Lumped)
        fail("Native mixed dynamics requires explicitly selected Lumped mass; Consistent mass is not supported by this increment.");
    if (!dynamics.getEnabled() || dynamics.getMethod() != RStructuralTimeIntegrationMethod::Newmark ||
        dynamics.getNewmarkBeta() != 0.25 || dynamics.getNewmarkGamma() != 0.5 ||
        dynamics.getRayleighAlpha() != 0.0 || dynamics.getRayleighBeta() != 0.0 ||
        !dynamics.getModalDampingRatios().empty())
        fail("Native mixed dynamics supports exactly undamped Newmark beta=.25, gamma=.5 without modal damping records.");
    if (problem.getStructuralCableSetup().enabled || problem.getStructuralPiezoSetup().enabled)
        fail("Native mixed dynamics cannot run with legacy cable or piezo modes enabled.");
    if (problem.getStructuralDynamicsRestartData().available)
        fail("Native mixed dynamics cannot consume legacy structural restart data.");
    const auto &control = problem.getStructuralControlSetup();
    if (control.enabled || !control.controllerId.isEmpty() || !control.channelConditioning.empty() ||
        control.deterministicSeed != 0 || control.controlProtocolVersion != 1)
        fail("Native mixed dynamics held commands do not support legacy control or callback protocols; v3 integration is pending.");
    const auto &thermal = problem.getStructuralThermalControlSetup();
    if (thermal.enabled || !thermal.channels.empty() || !thermal.sparseObservationNodeIds.empty() ||
        thermal.thermalSubstepSeconds != 0.0 || thermal.thermalSubstepsPerMechanicalMacrostep != 1 ||
        thermal.mechanicalMacrostepsPerControllerSample != 1)
        fail("Native mixed dynamics does not support thermal scheduling.");
    const auto &interaction = problem.getStructuralInteractionSetup();
    if (!interaction.getLinearConstraints().empty() || !interaction.getCouplings().empty() ||
        !interaction.getSensors().empty() || !interaction.getActuators().empty())
        fail("Native mixed dynamics does not support MPC, coupling or legacy sensor/actuator catalogs; v3 integration is pending.");
    if (model.findVariable(R_VARIABLE_TEMPERATURE) != RConstants::eod)
        fail("Native mixed dynamics is isothermal at 293.15 K and rejects a temperature field.");
    if (problem.getRestart() && !problem.getStructuralTrussDynamicsRestartData().available)
        fail("Native mixed restart requires a complete native accepted-state record.");
}

std::vector<uint> groupNodes(const RModel &model,const RElementGroup &group)
{
    std::set<uint> nodes;
    for (uint i=0;i<group.size();++i)
    {
        if (group.get(i) >= model.getNElements()) fail("Native mixed initial condition has an invalid element target.");
        const auto &element = model.getElement(group.get(i));
        for (uint j=0;j<element.size();++j)
        {
            if (element.getNodeId(j) >= model.getNNodes()) fail("Native mixed initial condition has an invalid node target.");
            nodes.insert(element.getNodeId(j));
        }
    }
    return {nodes.begin(),nodes.end()};
}

void sourceInitialConditions(const RModel &model,const RStructuralTrussDynamicsSystem &system,
                             double time,std::vector<double> &q,std::vector<double> &v)
{
    const auto n = system.prescribedMask.size();
    std::vector<bool> explicitQ(n,false),explicitV(n,false);
    q.assign(n,0.0); v.assign(n,0.0);
    for (uint gi=0;gi<model.getNElementGroups();++gi)
    {
        const auto *group = model.getElementGroupPtr(gi);
        if (!group || group->getNInitialConditions() == 0) continue;
        const auto nodes = groupNodes(model,*group);
        if (nodes.empty()) fail(QString("Native mixed initial conditions address empty physical group %1.").arg(group->getPhysicalId()));
        for (uint ci=0;ci<group->getNInitialConditions();++ci)
        {
            const auto &condition = group->getInitialCondition(ci);
            const bool displacement = condition.getType() == R_INITIAL_CONDITION_DISPLACEMENT;
            if (!displacement && condition.getType() != R_INITIAL_CONDITION_VELOCITY)
                fail(QString("Unsupported native mixed initial condition on physical group %1.").arg(group->getPhysicalId()));
            for (uint componentIndex=0;componentIndex<condition.size();++componentIndex)
            {
                const auto &component = condition.getComponent(componentIndex);
                if (!component.getEnabled()) continue;
                uint axis = 3;
                const auto type = component.getType();
                if (displacement)
                {
                    if (type == R_VARIABLE_DISPLACEMENT_X) axis=0;
                    else if (type == R_VARIABLE_DISPLACEMENT_Y) axis=1;
                    else if (type == R_VARIABLE_DISPLACEMENT_Z) axis=2;
                }
                else
                {
                    if (type == R_VARIABLE_VELOCITY_X) axis=0;
                    else if (type == R_VARIABLE_VELOCITY_Y) axis=1;
                    else if (type == R_VARIABLE_VELOCITY_Z) axis=2;
                }
                if (axis == 3 || component.size() != 1 || !std::isfinite(component.getValue(0)))
                    fail(QString("Native mixed initial condition on physical group %1 requires a constant finite global translation component.").arg(group->getPhysicalId()));
                auto &values = displacement ? q : v;
                auto &specified = displacement ? explicitQ : explicitV;
                const double value = component.getValue(0);
                for (uint node : nodes)
                {
                    const auto d = std::size_t(3)*node+axis;
                    if (specified[d] && values[d] != value)
                        fail(QString("Conflicting native mixed initial conditions at physical node %1, axis %2.").arg(model.getNodePhysicalId(node)).arg(axis));
                    specified[d] = true; values[d] = value;
                }
            }
        }
    }
    std::vector<double> prescribedQ = system.fixedPrescribedDisplacements,prescribedV(n,0.0);
    std::vector<bool> curveOwned(n,false);
    for (const auto &motion : system.prescribedMotions)
    {
        const auto state = motion.evaluate(time);
        const auto d = std::size_t(3)*motion.target.nodeId+uint(motion.target.dof);
        prescribedQ[d] = state.displacement; prescribedV[d] = state.velocity;
        curveOwned[d] = true;
    }
    for (std::size_t d=0;d<n;++d) if (system.prescribedMask[d])
    {
        const bool qMatches = curveOwned[d] ? analyticMotionEqual(q[d],prescribedQ[d]) : q[d] == prescribedQ[d];
        const bool vMatches = curveOwned[d] ? analyticMotionEqual(v[d],prescribedV[d]) : v[d] == prescribedV[d];
        if ((explicitQ[d] && !qMatches) || (explicitV[d] && !vMatches))
            fail(QString("Native mixed explicit initial condition conflicts with support/motion at physical node %1, axis %2.").arg(model.getNodePhysicalId(uint(d/3))).arg(uint(d%3)));
        q[d] = prescribedQ[d]; v[d] = prescribedV[d];
    }
}

void storeVector(RModel &model,RVariableType type,const std::vector<double> &values)
{
    uint position = model.findVariable(type);
    if (position == RConstants::eod) position = model.addVariable(type);
    auto &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(3,model.getNNodes());
    double maximum = 0.0;
    for (uint node=0;node<model.getNNodes();++node)
    {
        for (uint axis=0;axis<3;++axis) variable.setValue(axis,node,values[std::size_t(3)*node+axis]);
        maximum = std::max(maximum,std::hypot(std::hypot(values[3*node],values[3*node+1]),values[3*node+2]));
    }
    if (!std::isfinite(maximum)) fail("Native mixed kinematic display magnitude is nonfinite.");
    variable.getVariableData().setMinMaxDisplayValue(0.0,maximum);
}

bool sameInputs(const RStructuralTrussDynamicsInputs &a,const RStructuralTrussDynamicsInputs &b)
{
    if (a.externalForces != b.externalForces || a.memberInputs.size() != b.memberInputs.size()) return false;
    for (std::size_t i=0;i<a.memberInputs.size();++i)
    {
        const auto &left = a.memberInputs[i]; const auto &right = b.memberInputs[i];
        if (left.index() != right.index()) return false;
        if (const auto *value = std::get_if<RStructuralTrussVoltage>(&left))
        { if (value->value != std::get<RStructuralTrussVoltage>(right).value) return false; }
        else if (const auto *value = std::get_if<RStructuralTrussCharge>(&left))
        { if (value->value != std::get<RStructuralTrussCharge>(right).value) return false; }
        else if (const auto *value = std::get_if<RStructuralTrussCableLength>(&left))
        { if (value->value != std::get<RStructuralTrussCableLength>(right).value) return false; }
    }
    return true;
}
}

RStructuralTrussDynamicsModelAdapter::RStructuralTrussDynamicsModelAdapter(const RModel &model)
{
    preflight(model);
    frozenSetup = model.getProblemSetup().getStructuralTrussDynamicsSetup();
    frozenControlSpec = model.getProblemSetup().getStructuralNativeControlSpec();
    frozenReadoutSetup = model.getProblemSetup().getStructuralPiezoReadoutSetup();
    frozenFingerprint = frozenSetup.computeSourceFingerprint(model);
    if (frozenSetup.sourceFingerprint != frozenFingerprint)
        fail("Native mixed model must be bound to its exact current source before freezing.");
    const auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    frozenDt = dynamics.getTimeStep(); frozenTargetSteps = dynamics.getNSteps();
    frozenInitialTime = model.getTimeSolver().getInputStartTime()-frozenDt;
    frozenOptions.maxIterations = frozenSetup.maxIterations;
    frozenOptions.maxLineSearchIterations = frozenSetup.maxLineSearchIterations;
    frozenOptions.absForceTol = frozenSetup.absoluteForceTolerance;
    frozenOptions.relForceTol = frozenSetup.relativeForceTolerance;
    frozenOptions.linearRelativeTolerance = frozenSetup.linearRelativeTolerance;

    // Explicit mixed extraction policy; no unsupported setting is erased in a copy.
    const auto common = RStructuralStaticTrussModel::buildMixedDynamics(model);
    frozenSystem.referenceNodes = common.nodes;
    frozenSystem.prescribedMask = common.prescribedMask;
    frozenSystem.fixedPrescribedDisplacements = common.prescribedDisplacements;
    initialInputs.externalForces = common.deadNodalLoads;
    frozenSystem.nodalLumpedMass.assign(common.nodes.size(),0.0);
    for (uint node=0;node<model.getNNodes();++node)
    {
        const uint id = model.getNodePhysicalId(node);
        frozenMaps.orderedPhysicalNodeIds.push_back(id);
        frozenMaps.nodeByPhysicalId.emplace(id,node);
    }
    std::map<uint,const RStructuralTrussPassiveMember*> passive;
    std::map<uint,const RStructuralTrussActiveMember*> active;
    std::map<uint,const RStructuralTrussCableMember*> cables;
    for (const auto &member : frozenSetup.passiveMembers) passive.emplace(member.physicalElementId,&member);
    for (const auto &member : frozenSetup.activeMembers) active.emplace(member.physicalElementId,&member);
    for (const auto &member : frozenSetup.cableMembers) cables.emplace(member.physicalElementId,&member);
    std::map<uint,const RLine*> lineByElement;
    for (uint line=0;line<model.getNLines();++line)
        for (uint i=0;i<model.getLine(line).size();++i) lineByElement.emplace(model.getLine(line).get(i),&model.getLine(line));
    for (const auto &geometry : common.members)
    {
        const uint id = model.getElementPhysicalId(geometry.elementId);
        const auto index = frozenSystem.members.size();
        RStructuralTrussMember member;
        member.physicalElementId = id; member.node1 = geometry.node1; member.node2 = geometry.node2;
        RStructuralTrussMassSource massSource = RStructuralTrussMassSource::MaterialDensity;
        double mass = 0.0;
        if (passive.count(id))
        {
            member.law = RStructuralTrussPassiveLaw{passive.at(id)->referenceLength,geometry.youngModulus*geometry.area};
            initialInputs.memberInputs.emplace_back(RStructuralTrussNoInput{});
        }
        else if (active.count(id))
        {
            const auto &record = *active.at(id);
            member.law = RStructuralTrussPiezoLaw{record.referenceLength,record.axialStiffness,record.coupling,record.capacitance};
            if (record.kind == RStructuralTrussElectricalKind::Voltage)
                initialInputs.memberInputs.emplace_back(RStructuralTrussVoltage{record.initialValue});
            else initialInputs.memberInputs.emplace_back(RStructuralTrussCharge{record.initialValue});
            frozenMaps.orderedPortIds.push_back(record.portId);
            frozenMaps.memberByPortId.emplace(record.portId,index);
            massSource = record.massSource; mass = record.fixedMass;
        }
        else if (cables.count(id))
        {
            const auto &record = *cables.at(id);
            member.law = RStructuralTrussCableLaw{record.initialTension,geometry.youngModulus*geometry.area};
            initialInputs.memberInputs.emplace_back(RStructuralTrussCableLength{record.initialReferenceLength});
            massSource = record.massSource; mass = record.fixedMass;
        }
        else fail(QString("Native mixed physical member %1 lacks classification.").arg(id));
        if (massSource == RStructuralTrussMassSource::MaterialDensity)
        {
            const auto &material = lineByElement.at(geometry.elementId)->getMaterial();
            const uint densityIndex = material.findPosition(RMaterialProperty::Density);
            if (densityIndex >= material.size() || material.get(densityIndex).empty())
                fail(QString("Native mixed physical member %1 lacks density.").arg(id));
            mass = positiveProduct(material.get(densityIndex).get(293.15),geometry.area,geometry.referenceLength,
                QString("Native mixed density mass for physical member %1 is not finite positive.").arg(id));
        }
        const double half = positive(0.5*mass,QString("Native mixed half mass for physical member %1 is not finite positive.").arg(id));
        for (auto node : {member.node1,member.node2})
            frozenSystem.nodalLumpedMass[node] = positive(frozenSystem.nodalLumpedMass[node]+half,
                QString("Native mixed mass accumulation at physical node %1 is invalid.").arg(frozenMaps.orderedPhysicalNodeIds[node]));
        frozenMaps.orderedMemberPhysicalIds.push_back(id);
        frozenMaps.memberByPhysicalId.emplace(id,index);
        frozenSystem.members.push_back(member);
    }
    std::set<std::size_t> movingDofs;
    for (const auto &sourceMotion : model.getProblemSetup().getStructuralInteractionSetup().getPrescribedMotions())
    {
        sourceMotion.validate();
        auto motion = sourceMotion;
        const uint axis = uint(motion.target.dof);
        if (axis > 2 || motion.target.nodeId >= model.getNNodes())
            fail(QString("Native mixed motion '%1' requires a valid translational node target.").arg(motion.id));
        const auto d = std::size_t(3)*motion.target.nodeId+axis;
        if (!frozenSystem.prescribedMask[d] || !movingDofs.insert(d).second)
            fail(QString("Native mixed motion '%1' must uniquely own a held DOF at physical node %2.").arg(motion.id).arg(model.getNodePhysicalId(motion.target.nodeId)));
        motion.evaluate(frozenInitialTime); motion.evaluate(targetTime(frozenTargetSteps));
        frozenSystem.fixedPrescribedDisplacements[d] = 0.0;
        frozenSystem.prescribedMotions.push_back(std::move(motion));
    }
    sourceInitialConditions(model,frozenSystem,frozenInitialTime,initialQ,initialV);
    for (std::size_t d=0;d<frozenSystem.prescribedMask.size();++d)
        if (!frozenSystem.prescribedMask[d] && frozenSystem.nodalLumpedMass[d/3] <= 0.0)
            fail(QString("Native mixed free DOF at physical node %1, axis %2 has no positive frozen mass.")
                 .arg(frozenMaps.orderedPhysicalNodeIds[d/3]).arg(uint(d%3)));
    frozenDriver = std::make_unique<RStructuralTrussDynamics>(frozenSystem);
    if (frozenReadoutSetup.hasNonDefaultData())
        frozenReadout = std::make_unique<RStructuralPiezoNativeReadout>(*this,frozenReadoutSetup);
    if (frozenControlSpec.source == RStructuralNativeControlSource::LocalV3)
    {
        const RStructuralNativeControl provider(*this); // Resolve the actual compiled factory during preflight.
    }
}

RStructuralTrussDynamicsModelAdapter::~RStructuralTrussDynamicsModelAdapter() = default;

RStructuralTrussDynamicsState RStructuralTrussDynamicsModelAdapter::sourceInitialState() const
{
    return frozenDriver->initialize(frozenInitialTime,initialQ,initialV,initialInputs);
}

const RStructuralPiezoNativeReadout &RStructuralTrussDynamicsModelAdapter::piezoReadout() const
{
    if (!frozenReadout) fail("Native piezo readout requires a classified source catalog.");
    return *frozenReadout;
}

double RStructuralTrussDynamicsModelAdapter::targetTime(uint steps) const
{
    if (steps > frozenTargetSteps) fail("Native mixed completed step exceeds the frozen run target.");
    const double time = frozenInitialTime+double(steps)*frozenDt;
    if (!std::isfinite(time)) fail("Native mixed canonical physical time is nonfinite.");
    return time;
}

uint RStructuralTrussDynamicsModelAdapter::commandCursor(uint completedSteps) const
{
    if (completedSteps > frozenTargetSteps) fail("Native mixed command phase exceeds the frozen run target.");
    const auto position = std::lower_bound(frozenSetup.heldCommands.begin(),frozenSetup.heldCommands.end(),completedSteps+1,
        [](const RStructuralTrussHeldCommand &row,uint nextStep) { return row.stepIndex < nextStep; });
    return uint(position-frozenSetup.heldCommands.begin());
}

RStructuralTrussDynamicsState RStructuralTrussDynamicsModelAdapter::initialize() const
{
    if (frozenControlSpec.source == RStructuralNativeControlSource::LocalV3)
    {
        const RStructuralNativeControl control(*this);
        const auto inputs = control.completeInputs(control.initialize().lastAppliedValues);
        return driver().initialize(frozenInitialTime,initialQ,initialV,inputs);
    }
    return driver().initialize(frozenInitialTime,initialQ,initialV,initialInputs);
}

void RStructuralTrussDynamicsModelAdapter::validateCurrentSource(const RModel &model) const
{
    preflight(model);
    const auto &setup = model.getProblemSetup().getStructuralTrussDynamicsSetup();
    if (setup.computeSourceFingerprint(model) != frozenFingerprint || setup.sourceFingerprint != frozenFingerprint)
        fail("Native mixed current source differs from the frozen model/command identity.");
}

RStructuralTrussDynamicsAction RStructuralTrussDynamicsModelAdapter::actionForStep(
    const RStructuralTrussDynamicsInputs &accepted,uint cursor,uint step) const
{
    if (frozenControlSpec.source != RStructuralNativeControlSource::FrozenTable)
        fail("Native LocalV3 requires a complete provider proposal, not frozen command replay.");
    if (step == 0 || step > frozenTargetSteps || cursor != commandCursor(step-1))
        fail("Native mixed command cursor does not match the next one-based FE step.");
    if (!sameInputs(accepted,replayInputs(step-1).completeInputs))
        fail("Native mixed accepted applied inputs disagree with frozen command replay.");
    RStructuralTrussDynamicsAction action{accepted,cursor};
    if (cursor < frozenSetup.heldCommands.size() && frozenSetup.heldCommands[cursor].stepIndex == step)
    {
        const auto &row = frozenSetup.heldCommands[cursor];
        for (const auto &input : row.electrical)
        {
            auto &target = action.completeInputs.memberInputs.at(frozenMaps.memberByPortId.at(input.portId));
            if (std::holds_alternative<RStructuralTrussVoltage>(target)) target = RStructuralTrussVoltage{input.value};
            else target = RStructuralTrussCharge{input.value};
        }
        for (const auto &input : row.cables)
            action.completeInputs.memberInputs.at(frozenMaps.memberByPhysicalId.at(input.physicalElementId)) = RStructuralTrussCableLength{input.referenceLength};
        for (const auto &input : row.nodalForces)
        {
            const auto node = frozenMaps.nodeByPhysicalId.at(input.physicalNodeId);
            for (uint axis=0;axis<3;++axis) action.completeInputs.externalForces[3*node+axis] = input.force[axis];
        }
        ++action.nextCommandIndex;
    }
    return action;
}

RStructuralTrussDynamicsAction RStructuralTrussDynamicsModelAdapter::replayInputs(uint steps) const
{
    if (frozenControlSpec.source != RStructuralNativeControlSource::FrozenTable)
        fail("Native LocalV3 requires complete control state and cannot replay a frozen command table.");
    const uint cursor = commandCursor(steps);
    RStructuralTrussDynamicsAction replay{initialInputs,cursor};
    for (uint i=0;i<cursor;++i)
    {
        const auto &row = frozenSetup.heldCommands[i];
        for (const auto &input : row.electrical)
        {
            auto &target = replay.completeInputs.memberInputs.at(frozenMaps.memberByPortId.at(input.portId));
            if (std::holds_alternative<RStructuralTrussVoltage>(target)) target = RStructuralTrussVoltage{input.value};
            else target = RStructuralTrussCharge{input.value};
        }
        for (const auto &input : row.cables)
            replay.completeInputs.memberInputs.at(frozenMaps.memberByPhysicalId.at(input.physicalElementId)) = RStructuralTrussCableLength{input.referenceLength};
        for (const auto &input : row.nodalForces)
        {
            const auto node = frozenMaps.nodeByPhysicalId.at(input.physicalNodeId);
            for (uint axis=0;axis<3;++axis) replay.completeInputs.externalForces[3*node+axis] = input.force[axis];
        }
    }
    return replay;
}

RStructuralTrussDynamicsModelTrial RStructuralTrussDynamicsModelAdapter::trial(
    const RStructuralTrussDynamicsState &state,uint completed,uint cursor) const
{
    if (completed >= frozenTargetSteps || state.time != targetTime(completed))
        fail("Native mixed trial must follow the canonical accepted physical clock before the frozen target.");
    const auto action = actionForStep(state.appliedInputs,cursor,completed+1);
    const double target = targetTime(completed+1);
    const double dt = positive(target-state.time,"Native mixed canonical step interval is not representable and positive.");
    // Only binary64 grid rounding is allowed; the configured nominal dt stays frozen.
    if (!close(dt,frozenDt)) fail("Native mixed physical clock requires a step different from the frozen fixed time step.");
    RStructuralTrussDynamicsModelTrial candidate;
    candidate.result = driver().trialStep(state,action.completeInputs,dt,frozenOptions);
    if (!close(candidate.result.candidate.time,target)) fail("Native mixed driver returned a different canonical target time.");
    candidate.result.candidate.time = target;
    driver().refresh(candidate.result.candidate,action.completeInputs); // certificate only; preserve accepted Newmark a/ledgers.
    candidate.state = candidate.result.candidate;
    candidate.completeInputs = action.completeInputs;
    candidate.nextCommandIndex = action.nextCommandIndex;
    return candidate;
}

RModel RStructuralTrussDynamicsModelAdapter::candidateModel(const RModel &source,
    const RStructuralTrussDynamicsState &state,uint completed,uint cursor) const
{
    if (frozenControlSpec.source != RStructuralNativeControlSource::FrozenTable)
        fail("Native LocalV3 candidate requires its complete accepted control state.");
    return candidateModel(source,state,completed,cursor,RStructuralNativeControlState{});
}

RStructuralTrussDynamicsModelTrial RStructuralTrussDynamicsModelAdapter::trialWithCompleteInputs(
    const RStructuralTrussDynamicsState &state,uint completed,const RStructuralTrussDynamicsInputs &inputs) const
{
    if (frozenControlSpec.source != RStructuralNativeControlSource::LocalV3 || completed >= frozenTargetSteps ||
        state.time != targetTime(completed))
        fail("Native complete-input trial requires LocalV3 at the canonical accepted clock before its target.");
    const double target = targetTime(completed+1);
    const double dt = positive(target-state.time,"Native canonical physical interval is not positive and representable.");
    if (!close(dt,frozenDt)) fail("Native physical interval differs from the frozen nominal time step.");
    RStructuralTrussDynamicsModelTrial candidate;
    candidate.result = driver().trialStep(state,inputs,dt,frozenOptions);
    if (!close(candidate.result.candidate.time,target)) fail("Native LocalV3 trial returned a different target clock.");
    candidate.result.candidate.time = target;
    driver().refresh(candidate.result.candidate,inputs);
    candidate.state = candidate.result.candidate;
    candidate.completeInputs = inputs;
    candidate.nextCommandIndex = 0;
    return candidate;
}

RModel RStructuralTrussDynamicsModelAdapter::candidateModel(const RModel &source,
    const RStructuralTrussDynamicsState &state,uint completed,uint cursor,
    const RStructuralNativeControlState &controlState) const
{
    validateCurrentSource(source);
    const auto checkpoint = RStructuralTrussDynamicsRestartAdapter::build(*this,state,completed,cursor,controlState);
    const auto assembly = RStructuralTrussAssembly::evaluate(frozenSystem.referenceNodes,state.displacements,
        frozenSystem.members,state.appliedInputs.memberInputs);
    std::vector<double> axial;
    axial.reserve(assembly.memberResponses.size());
    for (const auto &member : assembly.memberResponses) axial.push_back(member.axialForce);
    RModel candidate(source);
    // build() has already certified analytic q/v/a at checkpoint.time, plus
    // the full accepted state and observations. Install its authoritative
    // native clock/state before publishing into this independent candidate.
    candidate.getProblemSetup().setStructuralTrussDynamicsRestartData(checkpoint);
    candidate.getProblemSetup().setStructuralNativeControlState(controlState);
    auto &time = candidate.getTimeSolver();
    // Use the const getter: the mutable getTimes() resets the accepted index.
    const auto &existingTimes = static_cast<const RTimeSolver&>(time).getTimes();
    bool canonicalTimes = existingTimes.size() == frozenTargetSteps;
    for (uint i=0;canonicalTimes && i<frozenTargetSteps;++i)
        canonicalTimes = existingTimes[i] == targetTime(i+1);
    if (!canonicalTimes)
    {
        std::vector<double> fullTimes;
        fullTimes.reserve(frozenTargetSteps);
        for (uint i=0;i<frozenTargetSteps;++i) fullTimes.push_back(targetTime(i+1));
        time.setTimes(fullTimes);
    }
    time.setEnabled(true);
    time.setInputNTimeSteps(frozenTargetSteps);
    time.setInputTimeStepSize(frozenDt);
    time.setOutputFrequency(source.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency());
    // Preserve hashed inputStartTime exactly, including binary64 rounding.
    // The zero-step baseline uses index zero but its computed clock is t0.
    time.setCurrentTimeStep(completed == 0 ? 0 : completed-1);
    time.setComputedTime(state.time);
    // The dedicated native mode takes all physical inputs from frozen source
    // records and ICs, never from display fields. Publish only its freshly
    // certified outputs; stale results from other physics are unavailable.
    // Temperature-field inputs have already been rejected by preflight.
    candidate.removeAllVariables();
    RStructuralStaticTrussModel::publishMechanicalMixed(candidate,state.displacements,assembly.internalForces,axial);
    storeVector(candidate,R_VARIABLE_VELOCITY,state.velocities);
    storeVector(candidate,R_VARIABLE_ACCELERATION,state.accelerations);
    return candidate;
}
