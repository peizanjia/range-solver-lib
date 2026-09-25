#ifndef RSTRUCTURAL_COUPLED_FIELD_SNAPSHOT_H
#define RSTRUCTURAL_COUPLED_FIELD_SNAPSHOT_H

#include <QString>
#include <QByteArray>

#include <array>
#include <optional>
#include <vector>

class RModel;

//! Read-only, physical-ID-keyed thermal/mechanical state on one exact mesh.
//!
//! Missing optional fields stay absent; they are never reported as zero. The
//! caller supplies persistent physical node and element IDs. Reordering mesh
//! storage is permitted when those IDs and their topology remain unchanged.
class RStructuralCoupledFieldSnapshot
{
    public:
        struct Node
        {
            QString id;
            double temperatureKelvin = 0.0;
            std::optional<double> temperatureRateKelvinPerSecond;
            std::optional<std::array<double,3>> displacementMetres;
            std::optional<std::array<double,3>> velocityMetresPerSecond;
            std::optional<std::array<double,3>> accelerationMetresPerSecondSquared;
        };

        struct Element
        {
            QString id;
            std::vector<double> heatFluxWattsPerSquareMetre;
            std::vector<double> heatSourceWatts;
            //! Solver components [xx,yy,zz,yz,xz,xy]; global for volumes,
            //! local element frame for surfaces and lines. Empty if absent.
            std::vector<double> stressPascals;
        };

        static RStructuralCoupledFieldSnapshot capture(
            const RModel &model,const std::vector<QString> &nodeIds,
            const std::vector<QString> &elementIds,
            const RStructuralCoupledFieldSnapshot *previous = nullptr);

        void verifyMesh(const RModel &model,const std::vector<QString> &nodeIds,
                        const std::vector<QString> &elementIds) const;
        std::vector<Node> selectNodes(const std::vector<QString> &ids) const;
        std::vector<Element> selectElements(const std::vector<QString> &ids) const;
        const std::vector<Node> &allNodes() const { return nodes; }
        const std::vector<Element> &allElements() const { return elements; }
        double timeSeconds() const { return time; }
        const QByteArray &meshDigest() const { return digest; }

    private:
        QByteArray digest;
        double time = 0.0;
        std::vector<Node> nodes;
        std::vector<Element> elements;
};

#endif // RSTRUCTURAL_COUPLED_FIELD_SNAPSHOT_H
