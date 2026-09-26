#include "rstructural_typed_fe_control.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

#include <QCryptographicHash>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTemporaryDir>

#include <rml_model.h>
#include <rml_point.h>

#include "rstructural_cable_model_adapter.h"

namespace
{
using Quantity = RStructuralTypedControlQuantity;

void require(bool condition,const char *message)
{
    if (!condition) throw std::invalid_argument(message);
}

RStructuralCableMember &member(RModel &model,uint elementId)
{
    auto &members = model.getProblemSetup().getStructuralCableSetup().members;
    auto found = std::find_if(members.begin(),members.end(),[elementId](const auto &value)
        { return value.elementId==elementId; });
    require(found!=members.end(),"Typed cable channel references an unconfigured member.");
    return *found;
}

RConditionComponent &anchor(RModel &model,uint nodeId,uint axis)
{
    require(nodeId<model.getNNodes() && axis<3,"Typed cable anchor references an invalid node/axis.");
    const RVariableType types[] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    RConditionComponent *result = nullptr;
    for (uint pointId=0;pointId<model.getNPoints();pointId++)
    {
        RPoint &point = model.getPoint(pointId);
        bool contains = false;
        for (uint i=0;i<point.size();i++)
        {
            const RElement &element = model.getElement(point.get(i));
            if (element.getType()==R_ELEMENT_POINT && element.getNodeId(0)==nodeId) contains=true;
        }
        if (!contains) continue;
        for (uint i=0;i<point.getNBoundaryConditions();i++)
        {
            RBoundaryCondition &condition = point.getBoundaryCondition(i);
            if (condition.getType()!=R_BOUNDARY_CONDITION_DISPLACEMENT) continue;
            const uint position = condition.findComponentPosition(types[axis]);
            if (position==RConstants::eod || !condition.getComponent(position).getEnabled()) continue;
            require(point.size()==1 && result==nullptr,
                    "Typed cable anchor requires one unique prescribed point component.");
            result = &condition.getComponent(position);
        }
    }
    require(result!=nullptr,"Typed cable anchor must already be a prescribed point component.");
    return *result;
}

uint patchPosition(const RStructuralPiezoBeamFEInput &input,uint elementId)
{
    for (uint i=0;i<input.patches.size();i++)
        if (input.patches[i].elementId==elementId) return i;
    throw std::invalid_argument("Typed piezo channel references an unconfigured beam patch.");
}

double currentCommand(RModel &model,const RStructuralPiezoBeamFEInput &input,
                      const RStructuralTypedFEControlBinding &binding,Quantity quantity)
{
    if (quantity==Quantity::CableReferenceLength) return member(model,binding.elementId).referenceLength;
    if (quantity==Quantity::CableInitialTension) return member(model,binding.elementId).initialTension;
    if (quantity==Quantity::CableAnchorDisplacement)
        return anchor(model,binding.nodeId,binding.axis).get(model.getTimeSolver().getCurrentTime());
    return input.patches[patchPosition(input,binding.elementId)].prescribedValue;
}

void applyCommand(RModel &model,RStructuralPiezoBeamFEInput &input,
                  const RStructuralTypedFEControlBinding &binding,Quantity quantity,double value)
{
    require(std::isfinite(value),"Typed FE applied command must be finite.");
    if (quantity==Quantity::CableReferenceLength)
    {
        require(value>0.0,"Controlled cable reference length must be positive.");
        member(model,binding.elementId).referenceLength=value;
    }
    else if (quantity==Quantity::CableInitialTension)
    {
        require(value>=0.0,"Controlled cable initial tension must be nonnegative.");
        member(model,binding.elementId).initialTension=value;
    }
    else if (quantity==Quantity::CableAnchorDisplacement)
    {
        RConditionComponent &component = anchor(model,binding.nodeId,binding.axis);
        component.clear();
        component.add(0.0,value);
    }
    else input.patches[patchPosition(input,binding.elementId)].prescribedValue=value;
}
}

struct RStructuralTypedFEControlDriver::Impl
{
    enum class Mode { Cable, PiezoBeam };
    struct Snapshot
    {
        RModel model;
        RStructuralPiezoBeamFEInput piezoInput;
        RStructuralPiezoBeamFEResult piezoResult;
        RStructuralCableEquilibriumResult cableResult;
        RStructuralTypedControlSession::State state;
        std::vector<RStructuralTypedControlChannelSample> channels;
    };

    Mode mode;
    RStructuralTypedControlSession session;
    std::vector<RStructuralTypedFEControlBinding> bindings;
    std::unique_ptr<Snapshot> accepted;
    bool advancing=false;

    Impl(Mode mode,const RModel &model,const RStructuralPiezoBeamFEInput &input,
         const RStructuralTypedControlSession &session,
         const std::vector<RStructuralTypedFEControlBinding> &bindings)
        : mode(mode),session(session),bindings(bindings),accepted(std::make_unique<Snapshot>())
    {
        accepted->model=model;
        accepted->piezoInput=input;
        accepted->state=session.initialState();
        require(!model.getTimeSolver().getEnabled() &&
                !model.getProblemSetup().getStructuralDynamicsSetup().getEnabled() &&
                !model.getProblemSetup().getStructuralControlSetup().enabled &&
                !model.getProblemSetup().getRestart() &&
                !model.getProblemSetup().getStructuralDynamicsRestartData().available &&
                !model.getProblemSetup().getStructuralThermalControlSetup().enabled,
                "Typed quasistatic FE driver requires an unscheduled static model.");
        // Only controller commands vary during these macrosteps. The native
        // static adapters evaluate their other conditions at the frozen model
        // time; temporal load/support curves must never be silently sampled at t0.
        for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
        {
            const auto *group=model.getElementGroupPtr(groupId);
            if (!group) continue;
            const auto check=[](const auto &condition)
            {
                for (uint i=0;i<condition.size();i++)
                {
                    const auto &component=condition.getComponent(i);
                    if (component.getEnabled()) require(component.size()==1,
                        "Typed quasistatic FE requires constant non-control loads/support/environment values.");
                }
            };
            for (uint i=0;i<group->getNBoundaryConditions();i++) check(group->getBoundaryCondition(i));
            for (uint i=0;i<group->getNEnvironmentConditions();i++) check(group->getEnvironmentCondition(i));
        }
        require(bindings.size()==session.catalog().size() && !bindings.empty(),
                "Typed FE channel bindings must match the ordered catalog.");
        std::set<QString> targets;
        for (uint i=0;i<bindings.size();i++)
        {
            const auto &channel=session.catalog()[i];
            const auto &binding=bindings[i];
            require(channel.channelId==binding.channelId,"Typed FE binding channel ID differs from catalog order.");
            const bool cableQuantity=channel.quantity==Quantity::CableReferenceLength ||
                channel.quantity==Quantity::CableInitialTension ||
                channel.quantity==Quantity::CableAnchorDisplacement;
            require(mode==Mode::Cable ? cableQuantity :
                channel.quantity==(input.electricalMode==RStructuralPiezoBeamElectricalMode::PrescribedVoltage
                    ? Quantity::Voltage : Quantity::Charge),
                "Typed FE quantity is incompatible with the cable/piezo formulation.");
            const QString target=QString::number(uint(channel.quantity))+":"+
                (channel.quantity==Quantity::CableAnchorDisplacement
                    ? QString::number(binding.nodeId)+":"+QString::number(binding.axis)
                    : QString::number(binding.elementId));
            require(targets.insert(target).second,"Typed FE catalog contains duplicate physical targets.");
            require(currentCommand(accepted->model,input,binding,channel.quantity)==channel.initialCommand,
                    "Typed FE initial command differs from the configured physical input.");
        }
        if (mode==Mode::Cable)
        {
            accepted->cableResult=RStructuralCableModelAdapter::solve(accepted->model);
            RStructuralCableModelAdapter::publish(accepted->model,accepted->cableResult);
        }
        else
        {
            accepted->piezoResult=RStructuralPiezoBeamFEAdapter::solve(accepted->model,input);
            RStructuralPiezoBeamFEAdapter::publish(accepted->model,input,accepted->piezoResult);
        }
    }

    std::vector<RStructuralControlSensorSample> sensors() const
    {
        std::vector<RStructuralControlSensorSample> values;
        const auto &displacements=mode==Mode::Cable ? accepted->cableResult.displacements :
            accepted->piezoResult.displacements;
        for (uint node=0;node<accepted->model.getNNodes();node++)
            for (uint axis=0;axis<3;axis++)
                values.push_back({QString("node-%1-u%2").arg(accepted->model.getNodePhysicalId(node)).arg(axis),
                                  displacements[3*node+axis],"m"});
        if (mode==Mode::PiezoBeam)
            for (uint i=0;i<accepted->piezoInput.patches.size();i++)
            {
                const QString id=QString("patch-%1").arg(accepted->model.getElementPhysicalId(
                    accepted->piezoInput.patches[i].elementId));
                values.push_back({id+"-voltage",accepted->piezoResult.patchVoltages[i],"V"});
                values.push_back({id+"-charge",accepted->piezoResult.patchCharges[i],"C"});
            }
        else
        {
            auto members=accepted->model.getProblemSetup().getStructuralCableSetup().members;
            std::sort(members.begin(),members.end(),[](const auto &a,const auto &b)
                { return a.elementId<b.elementId; });
            for (uint i=0;i<members.size();i++)
                values.push_back({QString("member-%1-tension").arg(
                    accepted->model.getElementPhysicalId(members[i].elementId)),
                    accepted->cableResult.elementTensions[i],"N"});
        }
        return values;
    }

    void advance(const Callback &callback)
    {
        require(!advancing,"Typed FE controller callbacks cannot reenter advance.");
        struct Guard
        {
            bool &flag;
            explicit Guard(bool &flag) : flag(flag) { flag=true; }
            ~Guard() { flag=false; }
        } guard(advancing);
        require(bool(callback),"Typed FE controller callback is required.");
        const auto observation=session.makeObservation(accepted->state,sensors());
        const auto decision=callback(observation);
        auto control=session.advance(accepted->state,observation,decision);
        auto candidate=std::make_unique<Snapshot>(*accepted);
        std::vector<double> work(bindings.size(),0.0);
        if (mode==Mode::Cable)
        {
            // A definite, catalog-ordered reversible path assigns the change
            // of stationary total potential to each parameter/support port.
            for (uint i=0;i<bindings.size();i++)
            {
                const double priorPotential=candidate->cableResult.energy;
                applyCommand(candidate->model,candidate->piezoInput,bindings[i],
                             session.catalog()[i].quantity,control.channels[i].appliedValue);
                candidate->cableResult=RStructuralCableModelAdapter::solve(candidate->model);
                work[i]=candidate->cableResult.energy-priorPotential;
            }
            RStructuralCableModelAdapter::publish(candidate->model,candidate->cableResult);
        }
        else
        {
            for (uint i=0;i<bindings.size();i++)
                applyCommand(candidate->model,candidate->piezoInput,bindings[i],
                             session.catalog()[i].quantity,control.channels[i].appliedValue);
            candidate->piezoResult=RStructuralPiezoBeamFEAdapter::solve(candidate->model,candidate->piezoInput);
            for (uint i=0;i<bindings.size();i++)
            {
                const uint patch=patchPosition(candidate->piezoInput,bindings[i].elementId);
                work[i]=0.5*(accepted->piezoResult.patchVoltages[patch]+candidate->piezoResult.patchVoltages[patch])*
                    (candidate->piezoResult.patchCharges[patch]-accepted->piezoResult.patchCharges[patch]);
            }
            RStructuralPiezoBeamFEAdapter::publish(candidate->model,candidate->piezoInput,candidate->piezoResult);
        }
        control=session.withWork(control,work);
        candidate->state=std::move(control.state);
        candidate->channels=std::move(control.channels);
        accepted.swap(candidate); // All throwing work precedes this noexcept commit.
    }
};

RStructuralTypedFEControlDriver::RStructuralTypedFEControlDriver(std::unique_ptr<Impl> implementation)
    : impl(std::move(implementation)) {}
RStructuralTypedFEControlDriver::~RStructuralTypedFEControlDriver() = default;

std::unique_ptr<RStructuralTypedFEControlDriver> RStructuralTypedFEControlDriver::cable(
    const RModel &model,const RStructuralTypedControlSession &session,
    const std::vector<RStructuralTypedFEControlBinding> &bindings)
{
    return std::unique_ptr<RStructuralTypedFEControlDriver>(new RStructuralTypedFEControlDriver(
        std::make_unique<Impl>(Impl::Mode::Cable,model,RStructuralPiezoBeamFEInput{},session,bindings)));
}

std::unique_ptr<RStructuralTypedFEControlDriver> RStructuralTypedFEControlDriver::piezoBeam(
    const RModel &model,const RStructuralPiezoBeamFEInput &input,
    const RStructuralTypedControlSession &session,
    const std::vector<RStructuralTypedFEControlBinding> &bindings)
{
    return std::unique_ptr<RStructuralTypedFEControlDriver>(new RStructuralTypedFEControlDriver(
        std::make_unique<Impl>(Impl::Mode::PiezoBeam,model,input,session,bindings)));
}

const RModel &RStructuralTypedFEControlDriver::model() const { return impl->accepted->model; }
const RStructuralTypedControlSession::State &RStructuralTypedFEControlDriver::state() const
{ return impl->accepted->state; }
const RStructuralTypedControlSession &RStructuralTypedFEControlDriver::session() const
{ return impl->session; }
const std::vector<RStructuralTypedControlChannelSample> &RStructuralTypedFEControlDriver::lastChannels() const
{ return impl->accepted->channels; }
RStructuralControlObservation RStructuralTypedFEControlDriver::observation() const
{ return impl->session.makeObservation(impl->accepted->state,impl->sensors()); }
void RStructuralTypedFEControlDriver::advance(const Callback &callback) { impl->advance(callback); }

namespace
{
constexpr qsizetype maxBundleBytes=128*1024*1024;

QString number(double value)
{
    if (value==std::numeric_limits<double>::infinity()) return "+inf";
    if (value==-std::numeric_limits<double>::infinity()) return "-inf";
    require(std::isfinite(value),"Typed FE configuration contains a non-finite number.");
    return QString::number(value,'g',17);
}

double readNumber(const QJsonValue &value,bool allowInfinity=false)
{
    require(value.isString(),"Typed FE configuration number must be a decimal string.");
    if (allowInfinity && value.toString()=="+inf") return std::numeric_limits<double>::infinity();
    if (allowInfinity && value.toString()=="-inf") return -std::numeric_limits<double>::infinity();
    bool ok=false;
    const double result=value.toString().toDouble(&ok);
    require(ok && std::isfinite(result),"Typed FE configuration has an invalid finite number.");
    return result;
}

uint readUInt(const QJsonValue &value,uint maximum=std::numeric_limits<uint>::max())
{
    const double result=value.toDouble(-1.0);
    require(value.isDouble() && result>=0.0 && result<=maximum && result==std::floor(result),
            "Typed FE configuration has an invalid unsigned integer.");
    return uint(result);
}

void appendBlock(QByteArray &bytes,const QByteArray &block)
{
    require(block.size()<=maxBundleBytes,"Typed FE checkpoint block is too large.");
    const quint32 count=quint32(block.size());
    for (uint i=0;i<4;i++) bytes.append(char((count>>(8*i))&255));
    bytes.append(block);
}

QByteArray takeBlock(const QByteArray &bytes,qsizetype &position)
{
    require(position>=0 && position<=bytes.size()-4,"Truncated typed FE checkpoint block.");
    quint32 count=0;
    for (uint i=0;i<4;i++) count|=quint32(quint8(bytes[position++]))<<(8*i);
    require(count<=maxBundleBytes && count<=bytes.size()-position,"Invalid typed FE checkpoint block size.");
    const auto block=bytes.mid(position,count);
    position+=count;
    return block;
}
}

void RStructuralTypedFEControlDriver::writeCheckpoint(const QString &fileName) const
{
    require(!fileName.isEmpty(),"Typed FE checkpoint requires a file name.");
    const auto &snapshot=*impl->accepted;
    const QByteArray stateBytes=impl->session.serialize(snapshot.state);
    QJsonObject config;
    config["mode"]=impl->mode==Impl::Mode::Cable ? 0 : 1;
    config["seed"]=QString::number(impl->session.deterministicSeed());
    config["dt"]=number(impl->session.fixedDt());
    config["controller"]=impl->session.controllerIdentity();
    QJsonArray channels;
    for (uint i=0;i<impl->bindings.size();i++)
    {
        const auto &channel=impl->session.catalog()[i];
        const auto &binding=impl->bindings[i];
        QJsonObject item;
        item["id"]=channel.channelId;
        item["quantity"]=int(channel.quantity);
        item["minimum"]=number(channel.minimum);
        item["maximum"]=number(channel.maximum);
        item["rate"]=number(channel.rateLimit);
        item["initial"]=number(channel.initialCommand);
        item["delay"]=int(channel.delaySteps);
        item["hold"]=double(channel.holdSteps);
        item["noise"]=number(channel.noiseAmplitude);
        item["element"]=double(binding.elementId);
        item["node"]=double(binding.nodeId);
        item["axis"]=int(binding.axis);
        channels.append(item);
    }
    config["channels"]=channels;
    QJsonArray patches;
    for (const auto &patch : snapshot.piezoInput.patches)
    {
        QJsonObject item;
        item["element"]=double(patch.elementId);
        item["area"]=number(patch.area);
        item["thickness"]=number(patch.thickness);
        item["offset"]=number(patch.centroidOffsetY);
        item["modulus"]=number(patch.youngModulus);
        item["coupling"]=number(patch.coupling);
        item["permittivity"]=number(patch.permittivity);
        item["value"]=number(patch.prescribedValue);
        patches.append(item);
    }
    config["patches"]=patches;
    config["electricMode"]=int(snapshot.piezoInput.electricalMode);
    const auto &input=snapshot.piezoInput;
    QJsonArray tolerances;
    for (double value : {input.translationalAbsoluteTolerance,input.translationalRelativeTolerance,
                         input.rotationalAbsoluteTolerance,input.rotationalRelativeTolerance,
                         input.electricAbsoluteTolerance,input.electricRelativeTolerance})
        tolerances.append(number(value));
    config["tolerances"]=tolerances;
    QTemporaryDir scratch;
    require(scratch.isValid(),"Could not create typed FE checkpoint staging directory.");
    const QString modelPath=scratch.filePath("model.rbm");
    require(snapshot.model.write(modelPath,false)==modelPath,"Typed FE model staging returned an unexpected path.");
    QFile modelFile(modelPath);
    require(modelFile.open(QIODevice::ReadOnly) && modelFile.size()>0 && modelFile.size()<=maxBundleBytes,
            "Could not read a bounded staged typed FE model.");
    const QByteArray modelBytes=modelFile.readAll();
    require(modelBytes.size()==modelFile.size(),"Typed FE staged model read was incomplete.");
    QByteArray payload("RSFEC001",8);
    appendBlock(payload,QJsonDocument(config).toJson(QJsonDocument::Compact));
    appendBlock(payload,modelBytes);
    appendBlock(payload,stateBytes);
    require(payload.size()<=maxBundleBytes-32,"Typed FE checkpoint exceeds 128 MiB.");
    payload.append(QCryptographicHash::hash(payload,QCryptographicHash::Sha256));
    QSaveFile output(fileName);
    require(output.open(QIODevice::WriteOnly) && output.write(payload)==payload.size() && output.commit(),
            "Could not atomically publish typed FE checkpoint.");
}

std::unique_ptr<RStructuralTypedFEControlDriver> RStructuralTypedFEControlDriver::readCheckpoint(
    const QString &fileName,const QString &expectedControllerIdentity)
{
    QFile file(fileName);
    require(file.open(QIODevice::ReadOnly) && file.size()>=52 && file.size()<=maxBundleBytes,
            "Typed FE checkpoint is missing or has an invalid size.");
    const QByteArray bytes=file.readAll();
    require(bytes.size()==file.size(),"Typed FE checkpoint read was incomplete.");
    const QByteArray payload=bytes.left(bytes.size()-32);
    require(payload.left(8)=="RSFEC001" &&
            QCryptographicHash::hash(payload,QCryptographicHash::Sha256)==bytes.right(32),
            "Typed FE checkpoint format or checksum differs.");
    qsizetype position=8;
    const QByteArray configBytes=takeBlock(payload,position);
    const QByteArray modelBytes=takeBlock(payload,position);
    const QByteArray stateBytes=takeBlock(payload,position);
    require(position==payload.size() && configBytes.size()<=4*1024*1024 && !modelBytes.isEmpty(),
            "Typed FE checkpoint has invalid or trailing content.");
    QJsonParseError parseError;
    const auto document=QJsonDocument::fromJson(configBytes,&parseError);
    require(parseError.error==QJsonParseError::NoError && document.isObject(),
            "Typed FE checkpoint configuration is not a JSON object.");
    const QJsonObject config=document.object();
    require(config["controller"].isString() && config["controller"].toString()==expectedControllerIdentity,
            "Typed FE checkpoint controller identity differs.");
    const uint mode=readUInt(config["mode"],1);
    bool seedOk=false;
    require(config["seed"].isString(),"Typed FE seed must be a decimal string.");
    const quint64 seed=config["seed"].toString().toULongLong(&seedOk);
    require(seedOk,"Invalid typed FE deterministic seed.");
    require(config["channels"].isArray(),"Typed FE checkpoint lacks a channel catalog.");
    const auto channelArray=config["channels"].toArray();
    require(!channelArray.empty() && channelArray.size()<=4096,"Typed FE channel count is invalid.");
    std::vector<RStructuralTypedControlChannel> catalog;
    std::vector<RStructuralTypedFEControlBinding> bindings;
    for (const QJsonValue &value : channelArray)
    {
        require(value.isObject(),"Typed FE channel is not an object.");
        const auto item=value.toObject();
        RStructuralTypedControlChannel channel;
        require(item["id"].isString(),"Typed FE channel ID is not a string.");
        channel.channelId=item["id"].toString();
        channel.quantity=Quantity(readUInt(item["quantity"],5));
        channel.minimum=readNumber(item["minimum"],true);
        channel.maximum=readNumber(item["maximum"],true);
        channel.rateLimit=readNumber(item["rate"],true);
        channel.initialCommand=readNumber(item["initial"]);
        channel.delaySteps=readUInt(item["delay"],16384);
        channel.holdSteps=readUInt(item["hold"]);
        channel.noiseAmplitude=readNumber(item["noise"]);
        catalog.push_back(channel);
        bindings.push_back({channel.channelId,readUInt(item["element"]),readUInt(item["node"]),readUInt(item["axis"],2)});
    }
    RStructuralTypedControlSession session(seed,readNumber(config["dt"]),catalog,expectedControllerIdentity);
    const auto state=session.deserialize(stateBytes);
    RStructuralPiezoBeamFEInput input;
    input.electricalMode=RStructuralPiezoBeamElectricalMode(readUInt(config["electricMode"],1));
    require(config["patches"].isArray() && config["patches"].toArray().size()<=4096,
            "Typed FE patch array is invalid.");
    for (const auto &value : config["patches"].toArray())
    {
        require(value.isObject(),"Typed FE patch is not an object.");
        const auto item=value.toObject();
        RStructuralPiezoBeamPatch patch;
        patch.elementId=readUInt(item["element"]);
        patch.area=readNumber(item["area"]);
        patch.thickness=readNumber(item["thickness"]);
        patch.centroidOffsetY=readNumber(item["offset"]);
        patch.youngModulus=readNumber(item["modulus"]);
        patch.coupling=readNumber(item["coupling"]);
        patch.permittivity=readNumber(item["permittivity"]);
        patch.prescribedValue=readNumber(item["value"]);
        input.patches.push_back(patch);
    }
    require(config["tolerances"].isArray() && config["tolerances"].toArray().size()==6,
            "Typed FE tolerance array is invalid.");
    const auto tolerances=config["tolerances"].toArray();
    input.translationalAbsoluteTolerance=readNumber(tolerances[0]);
    input.translationalRelativeTolerance=readNumber(tolerances[1]);
    input.rotationalAbsoluteTolerance=readNumber(tolerances[2]);
    input.rotationalRelativeTolerance=readNumber(tolerances[3]);
    input.electricAbsoluteTolerance=readNumber(tolerances[4]);
    input.electricRelativeTolerance=readNumber(tolerances[5]);
    QTemporaryDir scratch;
    require(scratch.isValid(),"Could not create typed FE restore directory.");
    const QString modelPath=scratch.filePath("model.rbm");
    QFile modelFile(modelPath);
    require(modelFile.open(QIODevice::WriteOnly) && modelFile.write(modelBytes)==modelBytes.size(),
            "Could not stage typed FE restore model.");
    modelFile.close();
    RModel model;
    model.read(modelPath);
    RModel initialModel=model;
    auto initialInput=input;
    for (uint i=0;i<bindings.size();i++)
    {
        require(currentCommand(model,input,bindings[i],catalog[i].quantity)==state.lastAppliedValues[i],
                "Typed FE checkpoint physical input differs from its accepted command.");
        applyCommand(initialModel,initialInput,bindings[i],catalog[i].quantity,catalog[i].initialCommand);
    }
    auto restored=mode==0 ? cable(initialModel,session,bindings) : piezoBeam(initialModel,initialInput,session,bindings);
    auto candidate=std::make_unique<Impl::Snapshot>();
    candidate->model=model;
    candidate->piezoInput=input;
    candidate->state=state;
    if (mode==0)
    {
        candidate->cableResult=RStructuralCableModelAdapter::solve(model);
        RStructuralCableModelAdapter::publish(candidate->model,candidate->cableResult);
    }
    else
    {
        candidate->piezoResult=RStructuralPiezoBeamFEAdapter::solve(model,input);
        RStructuralPiezoBeamFEAdapter::publish(candidate->model,input,candidate->piezoResult);
    }
    restored->impl->accepted.swap(candidate);
    return restored;
}
