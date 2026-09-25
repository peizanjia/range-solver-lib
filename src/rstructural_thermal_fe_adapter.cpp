#include "rstructural_thermal_fe_adapter.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include <QSet>
#include <QCryptographicHash>
#include <QDataStream>
#include <QFile>
#include <QSaveFile>
#include <QTemporaryDir>

#include <rml_model.h>
#include <rml_structural_interaction.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

#include "rsolverheat.h"
#include "rsolverstress.h"

namespace
{
constexpr double referenceKelvin = 293.15;
constexpr qint64 maxCheckpointBytes = 1024LL*1024LL*1024LL;

QString prefix(const RModel &model)
{
    // The generic heat solver temporarily moves reference nodes by the last
    // displacement, then subtracts it; floating roundoff must not invalidate
    // an otherwise unchanged binding. Quantize coordinates to 1e-12 m while
    // hashing every ordered node and connectivity record.
    QByteArray bytes;
    QDataStream stream(&bytes,QIODevice::WriteOnly);
    stream.setByteOrder(QDataStream::BigEndian);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    stream << quint32(model.getNNodes());
    for (unsigned int node=0;node<model.getNNodes();node++)
    {
        const auto &position = model.getNode(node);
        if (!std::isfinite(position.getX()) || !std::isfinite(position.getY()) ||
            !std::isfinite(position.getZ()))
            throw std::invalid_argument("Thermal FE mesh contains a non-finite coordinate.");
        const auto quantize = [](double coordinate) {
            if (std::abs(coordinate) > 1.0e12)
                throw std::invalid_argument("Thermal FE mesh coordinate exceeds ID binding range.");
            const double rounded = std::nearbyint(coordinate*1.0e12)/1.0e12;
            return rounded == 0.0 ? 0.0 : rounded;
        };
        stream << quantize(position.getX()) << quantize(position.getY())
               << quantize(position.getZ());
    }
    stream << quint32(model.getNElements());
    for (unsigned int element=0;element<model.getNElements();element++)
    {
        const auto &member = model.getElement(element);
        stream << quint32(member.getType()) << quint32(member.size());
        for (unsigned int local=0;local<member.size();local++)
            stream << quint32(member.getNodeId(local));
    }
    if (stream.status() != QDataStream::Ok)
        throw std::runtime_error("Could not fingerprint the thermal FE mesh.");
    return QString::fromLatin1(QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex())+
           QStringLiteral(":node:");
}

void require(bool condition,const char *message)
{
    if (!condition) throw std::invalid_argument(message);
}

std::vector<unsigned int> heaterPointElements(const RModel &model,
    const RStructuralThermalControlConfig &config)
{
    require(model.getProblemTaskTree().getProblemTypeMask() ==
                (R_PROBLEM_HEAT | R_PROBLEM_STRESS),
            "Thermal FE control requires a Heat->Stress task mask.");
    require(!model.getProblemSetup().getStructuralDynamicsSetup().getEnabled(),
            "Thermal FE control does not support structural dynamics.");
    require(model.getNNodes() >= 2,"Thermal FE control requires an axial rod mesh.");
    const double referenceY = model.getNode(0).getY();
    const double referenceZ = model.getNode(0).getZ();
    for (unsigned int node=0;node<model.getNNodes();node++)
        require(std::abs(model.getNode(node).getY()-referenceY) <= 1.0e-12 &&
                    std::abs(model.getNode(node).getZ()-referenceZ) <= 1.0e-12,
                "Thermal FE control supports x-aligned 1-D rods only.");
    for (unsigned int element=0;element<model.getNElements();element++)
        require(model.getElement(element).getType() == R_ELEMENT_TRUSS1 ||
                    model.getElement(element).getType() == R_ELEMENT_POINT,
                "Thermal FE control supports TRUSS1 and POINT elements only.");
    const auto &catalog = model.getProblemSetup().getStructuralInteractionSetup();
    const QString meshPrefix = prefix(model);
    std::vector<unsigned int> elements;
    elements.reserve(config.channels.size());
    QSet<unsigned int> usedElements;
    for (const auto &channel : config.channels)
    {
        const RStructuralActuator *actuator = nullptr;
        for (const auto &candidate : catalog.getActuators())
            if (candidate.id == channel.actuatorId)
            {
                require(!actuator,"Duplicate thermal actuator ID in model catalog.");
                actuator = &candidate;
            }
        require(actuator && actuator->type == RStructuralActuatorType::HeatSource &&
                    actuator->variable == R_VARIABLE_HEAT && actuator->unit == QStringLiteral("W") &&
                    actuator->target.scope == RStructuralTargetScope::Node &&
                    actuator->target.nodeId < model.getNNodes(),
                "Thermal control channel requires one node-targeted HeatSource actuator in W.");
        require(channel.heaterNodeId == meshPrefix + QString::number(actuator->target.nodeId),
                "Thermal actuator ID is not bound to this exact mesh and node index.");
        require(channel.minimumPowerWatts >= actuator->minimum &&
                    channel.maximumPowerWatts <= actuator->maximum &&
                    channel.maximumPowerSlewWattsPerSecond <= actuator->rateLimit,
                "Thermal scheduler power limits exceed the model actuator limits.");

        const RStructuralSensor *sensor = nullptr;
        for (const auto &candidate : catalog.getSensors())
            if (candidate.id == channel.temperatureSensorId)
            {
                require(!sensor,"Duplicate thermal sensor ID in model catalog.");
                sensor = &candidate;
            }
        require(sensor && sensor->quantity == RStructuralSensorQuantity::Temperature &&
                    sensor->aggregation == RStructuralSensorAggregation::Point &&
                    sensor->unit == QStringLiteral("K") &&
                    sensor->target.scope == RStructuralTargetScope::Node &&
                    sensor->target.nodeId < model.getNNodes(),
                "Thermal control channel requires one node-targeted point temperature sensor in K.");
        require(channel.temperatureSensorNodeId == meshPrefix + QString::number(sensor->target.nodeId),
                "Thermal sensor ID is not bound to this exact mesh and node index.");

        unsigned int pointElement = RConstants::eod;
        for (unsigned int element=0;element<model.getNElements();element++)
        {
            const auto &candidate = model.getElement(element);
            if (candidate.getType() == R_ELEMENT_POINT &&
                candidate.getNodeId(0) == actuator->target.nodeId)
            {
                require(pointElement == RConstants::eod,
                        "Thermal heater node has more than one POINT element.");
                pointElement = element;
            }
        }
        require(pointElement != RConstants::eod && !usedElements.contains(pointElement),
                "Thermal heater requires an unused POINT element at its target node.");
        usedElements.insert(pointElement);
        elements.push_back(pointElement);
    }
    return elements;
}

void verifyAcceptedPlant(
    const RStructuralThermalControlScheduler::Plant &actual,
    const RStructuralThermalControlScheduler::Plant &accepted)
{
    require(actual.size() == accepted.size(),"Thermal accepted state node count differs from RModel.");
    for (std::size_t node=0;node<actual.size();node++)
    {
        require(actual[node].nodeId == accepted[node].nodeId,
                "Thermal accepted state has a stale or reordered mesh ID.");
        const double temperatureError = std::abs(actual[node].temperatureKelvin-accepted[node].temperatureKelvin);
        const double displacementError = std::abs(actual[node].axialDisplacement-accepted[node].axialDisplacement);
        require(temperatureError <= 1.0e-10*std::max(1.0,std::abs(actual[node].temperatureKelvin)) &&
                    displacementError <= 1.0e-10*std::max(1.0,std::abs(actual[node].axialDisplacement)),
                "Thermal accepted state and RModel fields disagree.");
    }
}

void verifyAcceptedClock(const RModel &model,
                         const RStructuralThermalControlState &accepted)
{
    const double modelTime = model.getTimeSolver().getComputedTime();
    require(std::isfinite(modelTime) && std::isfinite(accepted.acceptedTimeSeconds) &&
                std::abs(modelTime-accepted.acceptedTimeSeconds) <=
                    1.0e-10*std::max({1.0,std::abs(modelTime),
                                     std::abs(accepted.acceptedTimeSeconds)}),
            "Thermal FE model time differs from the accepted controller state.");
}

void verifyModelScheduler(const RModel &model,
                          const RStructuralThermalControlScheduler &scheduler)
{
    if (!model.getProblemSetup().getStructuralThermalControlSetup().enabled)
        return;
    const auto persisted = RStructuralThermalFEAdapter::schedulerFromModel(model);
    require(persisted.configurationSha256() == scheduler.configurationSha256(),
            "Thermal checkpoint scheduler differs from persisted model control setup.");
}
}

QString RStructuralThermalFEAdapter::nodeId(const RModel &model,unsigned int nodeIndex)
{
    require(nodeIndex < model.getNNodes(),"Thermal FE node ID addresses an invalid node.");
    return prefix(model)+QString::number(nodeIndex);
}

RStructuralThermalControlScheduler RStructuralThermalFEAdapter::schedulerFromModel(
    const RModel &model)
{
    const auto &setup = model.getProblemSetup().getStructuralThermalControlSetup();
    require(setup.enabled,"Thermal FE control setup is disabled in RModel.");
    setup.validate();
    RStructuralThermalControlConfig config;
    config.thermalSubstepSeconds = setup.thermalSubstepSeconds;
    config.thermalSubstepsPerMechanicalMacrostep =
        setup.thermalSubstepsPerMechanicalMacrostep;
    config.mechanicalMacrostepsPerControllerSample =
        setup.mechanicalMacrostepsPerControllerSample;
    config.sparseObservationNodeIds = setup.sparseObservationNodeIds;
    config.channels.reserve(setup.channels.size());
    for (const auto &source : setup.channels)
    {
        RStructuralThermalControlChannelConfig channel;
        channel.actuatorId = source.actuatorId;
        channel.heaterNodeId = source.heaterNodeId;
        channel.temperatureSensorId = source.temperatureSensorId;
        channel.temperatureSensorNodeId = source.temperatureSensorNodeId;
        channel.turnOnBelowKelvin = source.turnOnBelowKelvin;
        channel.turnOffAtOrAboveKelvin = source.turnOffAtOrAboveKelvin;
        channel.onPowerWatts = source.onPowerWatts;
        channel.minimumPowerWatts = source.minimumPowerWatts;
        channel.maximumPowerWatts = source.maximumPowerWatts;
        channel.maximumPowerSlewWattsPerSecond = source.maximumPowerSlewUnlimited
            ? std::numeric_limits<double>::infinity()
            : source.maximumPowerSlewWattsPerSecond;
        config.channels.push_back(std::move(channel));
    }
    RStructuralThermalControlScheduler scheduler(std::move(config));
    (void)heaterPointElements(model,scheduler.getConfig());
    return scheduler;
}

RStructuralThermalControlScheduler::Plant RStructuralThermalFEAdapter::capture(const RModel &model)
{
    const unsigned int temperaturePosition = model.findVariable(R_VARIABLE_TEMPERATURE);
    const unsigned int displacementPosition = model.findVariable(R_VARIABLE_DISPLACEMENT);
    const RVariable *temperature = temperaturePosition == RConstants::eod
        ? nullptr : &model.getVariable(temperaturePosition);
    const RVariable *displacement = displacementPosition == RConstants::eod
        ? nullptr : &model.getVariable(displacementPosition);
    require(!temperature || (temperature->getApplyType() == R_VARIABLE_APPLY_NODE &&
                             temperature->getNVectors() == 1 &&
                             temperature->getNValues() == model.getNNodes()),
            "Thermal FE adapter requires one nodal temperature component.");
    require(!displacement || (displacement->getApplyType() == R_VARIABLE_APPLY_NODE &&
                              displacement->getNVectors() == 3 &&
                              displacement->getNValues() == model.getNNodes()),
            "Thermal FE adapter requires three nodal displacement components.");
    const QString meshPrefix = prefix(model);
    RStructuralThermalControlScheduler::Plant plant;
    plant.reserve(model.getNNodes());
    for (unsigned int node=0;node<model.getNNodes();node++)
    {
        const double kelvin = temperature ? temperature->getValue(0,node) : referenceKelvin;
        const double axial = displacement ? displacement->getValue(0,node) : 0.0;
        if (displacement)
            require(std::abs(displacement->getValue(1,node)) <= 1.0e-12 &&
                        std::abs(displacement->getValue(2,node)) <= 1.0e-12,
                    "Thermal FE control cannot omit transverse displacement state.");
        require(std::isfinite(kelvin) && kelvin >= 0.0 && std::isfinite(axial),
                "Thermal FE state has an invalid temperature or displacement.");
        plant.push_back({meshPrefix+QString::number(node),kelvin,axial});
    }
    std::sort(plant.begin(),plant.end(),[](const auto &a,const auto &b) {
        return a.nodeId < b.nodeId;
    });
    return plant;
}

RStructuralThermalControlScheduler::Result RStructuralThermalFEAdapter::advance(
    RModel &model,const RStructuralThermalControlScheduler &scheduler,
    const RStructuralThermalControlState &accepted)
{
    const auto &config = scheduler.getConfig();
    const auto pointElements = heaterPointElements(model,config);
    verifyAcceptedPlant(capture(model),accepted.plant);
    verifyAcceptedClock(model,accepted);
    RModel candidateModel(model);
    RSolverSharedData candidateSharedData;
    RSolverHeat heat(&candidateModel,QString(),QString(),candidateSharedData);
    heat.setReferenceLineThermalMass(true);
    RSolverStress stress(&candidateModel,QString(),QString(),candidateSharedData,false);
    unsigned int substep = 0;
    const double acceptedTime = model.getTimeSolver().getComputedTime();
    const auto result = scheduler.advance(accepted,[&](const auto &input,
        const auto &heaters,double dt,bool solveMechanical) {
        verifyAcceptedPlant(capture(candidateModel),input);
        require(heaters.size() == pointElements.size(),
                "Thermal scheduler returned an incomplete heater command vector.");
        std::vector<std::pair<uint,double>> powers;
        powers.reserve(heaters.size());
        for (std::size_t i=0;i<heaters.size();i++)
        {
            require(heaters[i].actuatorId == config.channels[i].actuatorId,
                    "Thermal scheduler heater order differs from FE binding.");
            powers.emplace_back(pointElements[i],heaters[i].watts);
        }
        RTimeSolver &time = candidateModel.getTimeSolver();
        time.setEnabled(true);
        time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);
        time.setInputTimeStepSize(dt);
        time.setTimes({acceptedTime + double(substep+1)*dt});
        time.setCurrentTimeStep(0);
        time.setComputedTime(time.getCurrentTime());
        heat.setControlledPointHeat(powers);
        // Each invocation uses a fresh solver. Recompute its geometric scales
        // on every FE substep, while retaining the last accepted temperature
        // except for a genuinely new run. RSolverGeneric treats its first
        // execution and first physical time step as separate concepts.
        const bool newRun = accepted.completedMechanicalMacrosteps == 0 && substep == 0;
        const bool originalRestart = candidateModel.getProblemSetup().getRestart();
        if (!newRun)
            candidateModel.getProblemSetup().setRestart(true);
        heat.run(true,0);
        candidateModel.getProblemSetup().setRestart(originalRestart);
        require(heat.hasConverged(),"Controlled thermal FE substep did not converge.");
        if (solveMechanical)
        {
            // This bounded adapter uses quasistatic mechanical equilibrium at
            // each macro boundary. Keep the heat substep clock, but do not
            // route the legacy time solver into a mechanical inertia term.
            time.setEnabled(false);
            stress.run(true,0);
            require(stress.hasConverged(),"Controlled mechanical FE macrostep did not converge.");
            time.setEnabled(true);
        }
        ++substep;
        return RStructuralThermalControlSubstepResult{capture(candidateModel),heaters};
    });
    require(prefix(candidateModel) == prefix(model),
            "Thermal FE macrostep changed the bound reference mesh.");
    model = candidateModel;
    return result;
}

void RStructuralThermalFEAdapter::writeCheckpoint(
    const QString &fileName,const RModel &model,
    const RStructuralThermalControlScheduler &scheduler,
    const RStructuralThermalControlState &accepted)
{
    require(!fileName.isEmpty(),"Thermal FE checkpoint requires a file name.");
    (void)heaterPointElements(model,scheduler.getConfig());
    verifyModelScheduler(model,scheduler);
    verifyAcceptedPlant(capture(model),accepted.plant);
    verifyAcceptedClock(model,accepted);
    const QByteArray stateBytes = scheduler.serialize(accepted);
    QTemporaryDir scratch;
    if (!scratch.isValid())
        throw std::runtime_error("Could not create a thermal checkpoint staging directory.");
    const QString modelPath = scratch.filePath(
        QStringLiteral("model.") + RModel::getDefaultFileExtension(true));
    model.write(modelPath,false);
    QFile modelFile(modelPath);
    if (!modelFile.open(QIODevice::ReadOnly) ||
        modelFile.size() <= 0 || modelFile.size() > maxCheckpointBytes)
        throw std::runtime_error("Could not read the staged thermal FE model.");
    const QByteArray modelBytes = modelFile.readAll();
    if (modelBytes.isEmpty())
        throw std::runtime_error("The staged thermal FE model is empty.");

    QByteArray payload;
    QDataStream stream(&payload,QIODevice::WriteOnly);
    stream.setVersion(QDataStream::Qt_6_0);
    stream.setByteOrder(QDataStream::BigEndian);
    stream << quint32(0x54464543) << quint32(1) << modelBytes << stateBytes;
    if (stream.status() != QDataStream::Ok)
        throw std::runtime_error("Could not encode the thermal FE checkpoint.");
    const QByteArray digest = QCryptographicHash::hash(payload,QCryptographicHash::Sha256);
    if (payload.size() > maxCheckpointBytes-digest.size())
        throw std::invalid_argument("Thermal FE checkpoint exceeds its readback size limit.");
    QSaveFile output(fileName);
    if (!output.open(QIODevice::WriteOnly) ||
        output.write(payload) != payload.size() ||
        output.write(digest) != digest.size() || !output.commit())
        throw std::runtime_error("Could not atomically publish the thermal FE checkpoint.");
}

namespace
{
void readCheckpointImpl(
    const QString &fileName,RModel &model,
    const RStructuralThermalControlScheduler *requiredScheduler,
    RStructuralThermalControlState &accepted)
{
    QFile input(fileName);
    if (!input.open(QIODevice::ReadOnly) || input.size() < 40 ||
        input.size() > maxCheckpointBytes)
        throw std::invalid_argument("Thermal FE checkpoint is missing or has an invalid size.");
    const QByteArray bytes = input.readAll();
    if (bytes.size() != input.size())
        throw std::runtime_error("Could not read the complete thermal FE checkpoint.");
    const QByteArray payload = bytes.left(bytes.size()-32);
    const QByteArray digest = bytes.right(32);
    if (QCryptographicHash::hash(payload,QCryptographicHash::Sha256) != digest)
        throw std::invalid_argument("Thermal FE checkpoint checksum differs from its content.");
    QDataStream stream(payload);
    stream.setVersion(QDataStream::Qt_6_0);
    stream.setByteOrder(QDataStream::BigEndian);
    quint32 magic = 0,version = 0;
    QByteArray modelBytes,stateBytes;
    stream >> magic >> version >> modelBytes >> stateBytes;
    if (stream.status() != QDataStream::Ok || !stream.atEnd() ||
        magic != 0x54464543 || version != 1 || modelBytes.isEmpty() || stateBytes.isEmpty())
        throw std::invalid_argument("Thermal FE checkpoint structure or version is invalid.");
    QTemporaryDir scratch;
    if (!scratch.isValid())
        throw std::runtime_error("Could not create a thermal checkpoint restore directory.");
    const QString modelPath = scratch.filePath(
        QStringLiteral("model.") + RModel::getDefaultFileExtension(true));
    QFile modelFile(modelPath);
    if (!modelFile.open(QIODevice::WriteOnly) ||
        modelFile.write(modelBytes) != modelBytes.size())
        throw std::runtime_error("Could not stage the thermal FE model for restoration.");
    modelFile.close();
    RModel restoredModel;
    restoredModel.read(modelPath);
    std::optional<RStructuralThermalControlScheduler> modelScheduler;
    if (!requiredScheduler)
        modelScheduler.emplace(RStructuralThermalFEAdapter::schedulerFromModel(restoredModel));
    const auto &scheduler = requiredScheduler ? *requiredScheduler : *modelScheduler;
    const auto restoredState = scheduler.deserialize(stateBytes);
    (void)heaterPointElements(restoredModel,scheduler.getConfig());
    verifyModelScheduler(restoredModel,scheduler);
    verifyAcceptedPlant(RStructuralThermalFEAdapter::capture(restoredModel),
                        restoredState.plant);
    verifyAcceptedClock(restoredModel,restoredState);
    model = restoredModel;
    accepted = restoredState;
}
}

void RStructuralThermalFEAdapter::readCheckpoint(
    const QString &fileName,RModel &model,
    const RStructuralThermalControlScheduler &scheduler,
    RStructuralThermalControlState &accepted)
{
    readCheckpointImpl(fileName,model,&scheduler,accepted);
}

void RStructuralThermalFEAdapter::readConfiguredCheckpoint(
    const QString &fileName,RModel &model,
    RStructuralThermalControlState &accepted)
{
    readCheckpointImpl(fileName,model,nullptr,accepted);
}
