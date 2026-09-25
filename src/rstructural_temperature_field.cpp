#include "rstructural_temperature_field.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <QCryptographicHash>
#include <QDataStream>
#include <QHash>
#include <QSet>

#include <rml_model.h>

namespace
{
void validateIds(const RModel &model,const std::vector<QString> &ids)
{
    if (ids.size() != model.getNNodes())
        throw std::invalid_argument("Temperature physical node ID count does not match the mesh.");
    QSet<QString> seen;
    for (const QString &id : ids)
    {
        if (id.isEmpty() || id != id.trimmed() || seen.contains(id))
            throw std::invalid_argument("Temperature physical node IDs must be unique, nonempty and unpadded.");
        seen.insert(id);
    }
}

void writeNode(QDataStream &stream,const QString &id,const RNode &node)
{
    if (!std::isfinite(node.getX()) || !std::isfinite(node.getY()) ||
        !std::isfinite(node.getZ()))
        throw std::invalid_argument("Temperature field mesh contains a non-finite coordinate.");
    stream << id.toUtf8() << node.getX() << node.getY() << node.getZ();
}
}

QByteArray RStructuralTemperatureField::fingerprint(
    const RModel &model,const std::vector<QString> &nodeIds)
{
    validateIds(model,nodeIds);
    std::vector<uint> nodeOrder(model.getNNodes());
    for (uint i=0;i<model.getNNodes();i++) nodeOrder[i] = i;
    std::sort(nodeOrder.begin(),nodeOrder.end(),[&](uint left,uint right) {
        return nodeIds[left] < nodeIds[right];
    });

    QByteArray encoded;
    QDataStream stream(&encoded,QIODevice::WriteOnly);
    stream.setVersion(QDataStream::Qt_6_0);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
    stream << quint32(1) << quint64(model.getNNodes()) << quint64(model.getNElements());
    for (uint node : nodeOrder) writeNode(stream,nodeIds[node],model.getNode(node));

    // An element is encoded with its native endpoint order. Sort the encoded
    // records so moving elements to other vector slots cannot change binding.
    std::vector<QByteArray> elements;
    elements.reserve(model.getNElements());
    for (uint elementId=0;elementId<model.getNElements();elementId++)
    {
        const RElement &element = model.getElement(elementId);
        QByteArray entry;
        QDataStream elementStream(&entry,QIODevice::WriteOnly);
        elementStream.setVersion(QDataStream::Qt_6_0);
        elementStream.setByteOrder(QDataStream::LittleEndian);
        elementStream << quint32(element.getType()) << quint64(element.size());
        for (uint local=0;local<element.size();local++)
        {
            const uint nodeId = element.getNodeId(local);
            if (nodeId >= nodeIds.size())
                throw std::invalid_argument("Temperature field mesh has an invalid element node ID.");
            elementStream << nodeIds[nodeId].toUtf8();
        }
        if (elementStream.status() != QDataStream::Ok)
            throw std::invalid_argument("Temperature element identity encoding failed.");
        elements.push_back(std::move(entry));
    }
    std::sort(elements.begin(),elements.end());
    for (const QByteArray &element : elements) stream << element;
    if (stream.status() != QDataStream::Ok)
        throw std::invalid_argument("Temperature mesh identity encoding failed.");
    return QCryptographicHash::hash(encoded,QCryptographicHash::Sha256).toHex();
}

RStructuralTemperatureField RStructuralTemperatureField::capture(
    const RModel &model,const std::vector<QString> &nodeIds,
    const std::vector<double> &nodalKelvin)
{
    if (nodalKelvin.size() != model.getNNodes())
        throw std::invalid_argument("Temperature sample count does not match mesh nodes.");
    RStructuralTemperatureField field;
    field.digest = fingerprint(model,nodeIds);
    field.records.reserve(nodeIds.size());
    for (uint node=0;node<model.getNNodes();node++)
    {
        const double temperature = nodalKelvin[node];
        if (!std::isfinite(temperature) || temperature < 0.0)
            throw std::invalid_argument("Temperature samples must be finite absolute kelvin.");
        field.records.push_back({nodeIds[node],temperature});
    }
    return field;
}

std::vector<double> RStructuralTemperatureField::mapTo(
    const RModel &model,const std::vector<QString> &nodeIds) const
{
    if (fingerprint(model,nodeIds) != this->digest)
        throw std::invalid_argument("Temperature field is bound to a different physical mesh.");
    if (records.size() != nodeIds.size())
        throw std::invalid_argument("Temperature field record count is inconsistent.");
    QHash<QString,double> byId;
    for (const Record &record : records)
    {
        if (byId.contains(record.nodeId) || !std::isfinite(record.kelvin) || record.kelvin < 0.0)
            throw std::invalid_argument("Temperature field contains duplicate or invalid records.");
        byId.insert(record.nodeId,record.kelvin);
    }
    std::vector<double> ordered;
    ordered.reserve(nodeIds.size());
    for (const QString &id : nodeIds)
    {
        if (!byId.contains(id))
            throw std::invalid_argument("Temperature field lacks a destination physical node ID.");
        ordered.push_back(byId.value(id));
    }
    return ordered;
}
