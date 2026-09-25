#include "rstructural_coupled_field_snapshot.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include <QCryptographicHash>
#include <QDataStream>
#include <QHash>
#include <QSet>

#include <rml_model.h>
#include <rml_variable.h>

#include "rstructural_temperature_field.h"

namespace
{
void validateIds(const std::vector<QString> &ids,std::size_t expected,const char *name)
{
    if (ids.size() != expected)
        throw std::invalid_argument(std::string(name)+" count differs from model storage.");
    QSet<QString> seen;
    for (const QString &id : ids)
    {
        if (id.isEmpty() || id != id.trimmed() || seen.contains(id))
            throw std::invalid_argument(std::string(name)+" IDs must be unique, nonempty and unpadded.");
        seen.insert(id);
    }
}

QByteArray fingerprint(const RModel &model,const std::vector<QString> &nodeIds,
                       const std::vector<QString> &elementIds)
{
    validateIds(nodeIds,model.getNNodes(),"Physical node");
    validateIds(elementIds,model.getNElements(),"Physical element");
    const std::vector<double> placeholder(model.getNNodes(),293.15);
    const QByteArray nodeAndMesh = RStructuralTemperatureField::capture(
        model,nodeIds,placeholder).meshDigest();
    std::vector<QByteArray> bindings;
    bindings.reserve(model.getNElements());
    for (uint i=0;i<model.getNElements();i++)
    {
        const RElement &element = model.getElement(i);
        QByteArray entry;
        QDataStream stream(&entry,QIODevice::WriteOnly);
        stream.setVersion(QDataStream::Qt_6_0);
        stream.setByteOrder(QDataStream::LittleEndian);
        stream << elementIds[i] << quint32(element.getType()) << quint32(element.size());
        for (uint local=0;local<element.size();local++)
        {
            const uint node = element.getNodeId(local);
            if (node >= nodeIds.size())
                throw std::invalid_argument("Coupled field element has an invalid node index.");
            stream << nodeIds[node];
        }
        if (stream.status() != QDataStream::Ok)
            throw std::runtime_error("Could not encode a coupled field element binding.");
        bindings.push_back(std::move(entry));
    }
    std::sort(bindings.begin(),bindings.end());
    QByteArray bytes;
    QDataStream stream(&bytes,QIODevice::WriteOnly);
    stream.setVersion(QDataStream::Qt_6_0);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream << quint32(1) << nodeAndMesh << quint64(bindings.size());
    for (const auto &entry : bindings) stream << entry;
    if (stream.status() != QDataStream::Ok)
        throw std::runtime_error("Could not fingerprint the coupled physical mesh.");
    return QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex();
}

const RVariable *field(const RModel &model,RVariableType type,
                       RVariableApplyType apply,uint count,bool mandatory = false)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        if (mandatory)
            throw std::invalid_argument("Coupled field requires a persisted nodal temperature variable.");
        return nullptr;
    }
    const RVariable &variable = model.getVariable(position);
    if (variable.getApplyType() != apply || variable.getNValues() != count ||
        variable.getNVectors() == 0)
        throw std::invalid_argument("Coupled field variable has an unexpected scope or shape.");
    return &variable;
}

std::array<double,3> vector3(const RVariable &variable,uint index)
{
    if (variable.getNVectors() != 3)
        throw std::invalid_argument("Coupled field vector must have three components.");
    std::array<double,3> result{};
    for (uint component=0;component<3;component++)
    {
        result[component] = variable.getValue(component,index);
        if (!std::isfinite(result[component]))
            throw std::invalid_argument("Coupled field contains a non-finite vector component.");
    }
    return result;
}

std::vector<double> values(const RVariable *variable,uint index)
{
    if (!variable) return {};
    std::vector<double> result;
    result.reserve(variable->getNVectors());
    for (uint component=0;component<variable->getNVectors();component++)
    {
        const double value = variable->getValue(component,index);
        if (!std::isfinite(value))
            throw std::invalid_argument("Coupled field contains a non-finite element result.");
        result.push_back(value);
    }
    return result;
}

template <typename Record>
std::vector<Record> select(const std::vector<Record> &all,const std::vector<QString> &ids)
{
    QHash<QString,const Record*> byId;
    for (const auto &record : all) byId.insert(record.id,&record);
    QSet<QString> seen;
    std::vector<Record> selected;
    selected.reserve(ids.size());
    for (const QString &id : ids)
    {
        if (seen.contains(id) || !byId.contains(id))
            throw std::invalid_argument("Coupled sparse selection has a duplicate or unknown physical ID.");
        seen.insert(id);
        selected.push_back(*byId.value(id));
    }
    std::sort(selected.begin(),selected.end(),[](const auto &left,const auto &right) {
        return left.id < right.id;
    });
    return selected;
}

std::vector<QString> persistedNodeIds(const RModel &model)
{
    std::vector<QString> ids;
    ids.reserve(model.getNNodes());
    for (uint i=0;i<model.getNNodes();i++)
        ids.push_back(QString::number(model.getNodePhysicalId(i)));
    return ids;
}

std::vector<QString> persistedElementIds(const RModel &model)
{
    std::vector<QString> ids;
    ids.reserve(model.getNElements());
    for (uint i=0;i<model.getNElements();i++)
        ids.push_back(QString::number(model.getElementPhysicalId(i)));
    return ids;
}
}

RStructuralCoupledFieldSnapshot RStructuralCoupledFieldSnapshot::capture(
    const RModel &model,const RStructuralCoupledFieldSnapshot *previous)
{
    return capture(model,persistedNodeIds(model),persistedElementIds(model),previous);
}

RStructuralCoupledFieldSnapshot RStructuralCoupledFieldSnapshot::capture(
    const RModel &model,const std::vector<QString> &nodeIds,
    const std::vector<QString> &elementIds,
    const RStructuralCoupledFieldSnapshot *previous)
{
    RStructuralCoupledFieldSnapshot snapshot;
    snapshot.digest = fingerprint(model,nodeIds,elementIds);
    snapshot.time = model.getTimeSolver().getComputedTime();
    if (!std::isfinite(snapshot.time))
        throw std::invalid_argument("Coupled field time must be finite seconds.");
    if (previous && (previous->digest != snapshot.digest ||
                     previous->time >= snapshot.time))
        throw std::invalid_argument("Coupled field rate requires an earlier snapshot on the same physical mesh.");
    const RVariable *temperature = field(model,R_VARIABLE_TEMPERATURE,
        R_VARIABLE_APPLY_NODE,model.getNNodes(),true);
    if (temperature->getNVectors() != 1)
        throw std::invalid_argument("Coupled field temperature must have one nodal component.");
    const RVariable *displacement = field(model,R_VARIABLE_DISPLACEMENT,
        R_VARIABLE_APPLY_NODE,model.getNNodes());
    const RVariable *velocity = field(model,R_VARIABLE_VELOCITY,
        R_VARIABLE_APPLY_NODE,model.getNNodes());
    const RVariable *acceleration = field(model,R_VARIABLE_ACCELERATION,
        R_VARIABLE_APPLY_NODE,model.getNNodes());
    QHash<QString,double> priorTemperature;
    if (previous)
        for (const auto &record : previous->nodes)
            priorTemperature.insert(record.id,record.temperatureKelvin);
    snapshot.nodes.reserve(nodeIds.size());
    for (uint i=0;i<model.getNNodes();i++)
    {
        Node node;
        node.id = nodeIds[i];
        node.temperatureKelvin = temperature->getValue(0,i);
        if (!std::isfinite(node.temperatureKelvin) || node.temperatureKelvin < 0.0)
            throw std::invalid_argument("Coupled field requires finite absolute kelvin.");
        if (previous)
        {
            if (!priorTemperature.contains(node.id))
                throw std::invalid_argument("Previous coupled field lacks a physical node ID.");
            node.temperatureRateKelvinPerSecond =
                (node.temperatureKelvin-priorTemperature.value(node.id))/
                (snapshot.time-previous->time);
        }
        if (displacement) node.displacementMetres = vector3(*displacement,i);
        if (velocity) node.velocityMetresPerSecond = vector3(*velocity,i);
        if (acceleration) node.accelerationMetresPerSecondSquared = vector3(*acceleration,i);
        snapshot.nodes.push_back(std::move(node));
    }
    std::sort(snapshot.nodes.begin(),snapshot.nodes.end(),[](const auto &a,const auto &b) {
        return a.id < b.id;
    });

    const RVariable *heatFlux = field(model,R_VARIABLE_HEAT_FLUX,
        R_VARIABLE_APPLY_ELEMENT,model.getNElements());
    const RVariable *heatSource = field(model,R_VARIABLE_HEAT,
        R_VARIABLE_APPLY_ELEMENT,model.getNElements());
    const RVariableType stressTypes[6] = {
        R_VARIABLE_STRESS_X,R_VARIABLE_STRESS_Y,R_VARIABLE_STRESS_Z,
        R_VARIABLE_STRESS_YZ,R_VARIABLE_STRESS_XZ,R_VARIABLE_STRESS_XY
    };
    std::array<const RVariable*,6> stressComponents{};
    uint presentComponents = 0;
    for (uint component=0;component<stressComponents.size();component++)
    {
        stressComponents[component] = field(model,stressTypes[component],
            R_VARIABLE_APPLY_ELEMENT,model.getNElements());
        if (stressComponents[component])
        {
            if (stressComponents[component]->getNVectors() != 1)
                throw std::invalid_argument("Coupled field stress component must be scalar.");
            presentComponents++;
        }
    }
    if (presentComponents != 0 && presentComponents != stressComponents.size())
        throw std::invalid_argument("Coupled field requires a complete six-component stress tensor.");
    snapshot.elements.reserve(elementIds.size());
    for (uint i=0;i<model.getNElements();i++)
    {
        std::vector<double> stress;
        if (presentComponents)
        {
            stress.reserve(stressComponents.size());
            for (const RVariable *component : stressComponents)
                stress.push_back(values(component,i).front());
        }
        snapshot.elements.push_back({elementIds[i],values(heatFlux,i),
                                     values(heatSource,i),std::move(stress)});
    }
    std::sort(snapshot.elements.begin(),snapshot.elements.end(),[](const auto &a,const auto &b) {
        return a.id < b.id;
    });
    return snapshot;
}

void RStructuralCoupledFieldSnapshot::verifyMesh(
    const RModel &model,const std::vector<QString> &nodeIds,
    const std::vector<QString> &elementIds) const
{
    if (fingerprint(model,nodeIds,elementIds) != digest)
        throw std::invalid_argument("Coupled field snapshot is bound to a different physical mesh.");
}

void RStructuralCoupledFieldSnapshot::verifyMesh(const RModel &model) const
{
    verifyMesh(model,persistedNodeIds(model),persistedElementIds(model));
}

std::vector<RStructuralCoupledFieldSnapshot::Node>
RStructuralCoupledFieldSnapshot::selectNodes(const std::vector<QString> &ids) const
{
    return select(nodes,ids);
}

std::vector<RStructuralCoupledFieldSnapshot::Element>
RStructuralCoupledFieldSnapshot::selectElements(const std::vector<QString> &ids) const
{
    return select(elements,ids);
}
