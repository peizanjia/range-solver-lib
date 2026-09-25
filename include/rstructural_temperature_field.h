#ifndef RSTRUCTURAL_TEMPERATURE_FIELD_H
#define RSTRUCTURAL_TEMPERATURE_FIELD_H

#include <QByteArray>
#include <QString>

#include <vector>

class RModel;

//! Temperature snapshot keyed by caller-assigned, stable physical node IDs.
//!
//! The binding hashes exact coordinates and oriented connectivity after
//! replacing model indices with physical IDs. Node/element reorder is safe;
//! geometry, connectivity or an ID substitution is rejected. The caller must
//! persist/reconstruct the same physical IDs when transferring between models.
class RStructuralTemperatureField
{
    public:
        struct Record
        {
            QString nodeId;
            double kelvin = 0.0;
        };

        static RStructuralTemperatureField capture(
            const RModel &model,
            const std::vector<QString> &nodeIds,
            const std::vector<double> &nodalKelvin);

        //! Return temperatures in the destination model's node-index order.
        std::vector<double> mapTo(
            const RModel &model,
            const std::vector<QString> &nodeIds) const;

        const QByteArray &meshDigest() const { return digest; }
        const std::vector<Record> &values() const { return records; }

    private:
        QByteArray digest;
        std::vector<Record> records;
        static QByteArray fingerprint(const RModel &model,
                                      const std::vector<QString> &nodeIds);
};

#endif // RSTRUCTURAL_TEMPERATURE_FIELD_H
