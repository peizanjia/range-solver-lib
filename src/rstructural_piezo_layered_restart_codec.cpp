#include "rstructural_piezo_layered_restart_codec.h"
#include <rbl_hash_sha256.h>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace {
constexpr uint bound=16*1024*1024, vectorBound=1024*1024;
void require(bool b) {if(!b)throw std::invalid_argument("Invalid layered canonical restart bytes");}
struct Archive {
    QByteArray bytes; bool reading=false, comparing=false, equal=true; uint offset=0;
    const QByteArray *comparison=nullptr; uint comparisonSize=0;
    void compareBytes(const char *data,uint size) {
        require(offset<=bound && size<=bound-offset);
        const uint available=offset<comparisonSize?comparisonSize-offset:0;
        const uint common=std::min(size,available);
        const char *other=common?comparison->constData()+offset:nullptr;
        if(common && std::memcmp(data,other,size_t(common))!=0)equal=false;
        if(common!=size)equal=false;
        offset+=size;
    }
    void integer(quint64 &v) {
        if(reading){
            require(offset+8u<=uint(bytes.size()));
            v=qFromLittleEndian<quint64>(bytes.constData()+offset);offset+=8u;
        }else if(comparing){
            char word[8];qToLittleEndian<quint64>(v,word);compareBytes(word,8);
        }else{
            char word[8];qToLittleEndian<quint64>(v,word);bytes.append(word,8);
            require(uint(bytes.size())<=bound);
        }
    }
    void count(uint &v,uint maximum=vectorBound) {quint64 x=v;integer(x);require(x<=maximum);v=uint(x);}
    void real(double &v) {static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE binary64 required");quint64 x=0;if(!reading)std::memcpy(&x,&v,8);integer(x);if(reading)std::memcpy(&v,&x,8);require(std::isfinite(v));}
    void boolean(bool &v) {quint64 x=v?1:0;integer(x);require(x<=1);v=x!=0;}
    void blob(QByteArray &v) {require(v.size()<=qsizetype(bound));uint n=uint(v.size());count(n,bound);if(reading){require(n<=uint(bytes.size())-offset);v=bytes.mid(int(offset),int(n));offset+=n;}else if(comparing){compareBytes(v.constData(),n);}else{bytes.append(v);require(uint(bytes.size())<=bound);}}
    void text(QString &v) {QByteArray b=v.toUtf8();blob(b);if(reading){v=QString::fromUtf8(b);require(v.toUtf8()==b);}else require(QString::fromUtf8(b)==v);}
    template<class E> void enumeration(E &v,uint max) {uint n=uint(v);count(n,max);v=E(n);}
    template<class V> void vector(V &v) {
        require(v.size()<=vectorBound);uint n=uint(v.size());count(n);
        if(reading)v.resize(n);
#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        if(!reading){
            static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE binary64 required");
            // Validate every element even after an earlier comparison mismatch.
            for(double x:v)require(std::isfinite(x));
            const uint size=n*8u; // count() limits n to vectorBound.
            if(comparing)compareBytes(reinterpret_cast<const char *>(v.data()),size);
            else{
                require(uint(bytes.size())<=bound && size<=bound-uint(bytes.size()));
                if(size)bytes.append(reinterpret_cast<const char *>(v.data()),int(size));
            }
            return;
        }
#endif
        // Decoder and big-endian hosts retain the canonical per-word path.
        for(double &x:v)real(x);
    }
    template<class V,class F> void entries(V &v,F f) {uint n=uint(v.size());count(n,65536);if(reading)v.resize(n);for(auto &x:v)f(x);}
};
void physical(Archive &a,RStructuralPiezoLayeredDynamicsState &s) {
    a.text(s.sourceFingerprint);a.text(s.forcingIdentity);a.blob(s.driverIdentity);a.blob(s.integrityDigest);
    auto &c=s.integrator;auto &r=c.state;
    a.real(r.time);a.vector(r.q);a.vector(r.v);a.vector(r.a);a.vector(r.externalForce);
    for(double *x:{&r.kineticEnergy,&r.strainEnergy,&r.externalWork,&r.dampingDissipation,&r.energyBalanceResidual,&c.fixedStep,&c.initialMechanicalEnergy})a.real(*x);
    for(auto *v:{&s.heldActuatorVoltages,&s.initialActuatorVoltages,&s.perturbationForces,&s.physicalDisplacements,&s.physicalVelocities,&s.physicalAccelerations,&s.physicalReactions,&s.electrical.internalForce,&s.electrical.deltaInternalForce})a.vector(*v);
    a.entries(s.electrical.ports,[&](RStructuralPiezoElectricalPortResponse &p){a.text(p.portId);a.enumeration(p.role,3);for(double *x:{&p.voltage,&p.charge,&p.deltaVoltage,&p.deltaCharge,&p.amplifierVoltage,&p.deltaAmplifierVoltage})a.real(*x);});
    a.entries(s.readouts,[&](RStructuralPiezoLayeredDynamicsReadout &p){a.text(p.portId);a.text(p.unit);a.enumeration(p.mode,4);a.boolean(p.actionRecord);a.real(p.value);a.real(p.deltaValue);});
    for(double *x:{&s.kineticEnergy,&s.storedEnergy,&s.initialEnergy,&s.backgroundExternalWork,&s.electricalWork,&s.viscousDissipation,&s.energyBalanceResidual,&s.freeForceResidual})a.real(*x);
}
void record(Archive &a,RStructuralPiezoLayeredRunState &s) {
    QString magic="range-fea/layered-native-restart/binary64/v1";a.text(magic);require(magic=="range-fea/layered-native-restart/binary64/v1");
    a.text(s.sourceFingerprint);a.text(s.forcingIdentity);a.text(s.controllerIdentity);a.text(s.backendIdentity);
    a.blob(s.configurationIdentity);a.blob(s.opaqueControllerState);a.blob(s.lastFrameBytes);
    for(uint *x:{&s.targetSteps,&s.completedFE,&s.controlEveryFESteps,&s.controlTick,&s.controlPhase})a.count(*x);
    a.real(s.initialTime);a.real(s.dt);a.real(s.nominalControlTs);
    a.entries(s.envelopes,[&](RStructuralPiezoLayeredVoltageEnvelope &e){a.real(e.minimum);a.real(e.maximum);a.real(e.maximumRate);require(e.minimum<=e.maximum && e.maximumRate>=0);});
    a.vector(s.heldRequested);a.vector(s.lastApplied);
    a.entries(s.actions,[&](RStructuralPiezoLayeredVoltageAction &v){a.text(v.portId);for(double *x:{&v.requested,&v.amplitudeLimited,&v.previousApplied,&v.applied})a.real(*x);a.boolean(v.amplitudeClipped);a.boolean(v.rateClipped);});
    physical(a,s.physical);
    require(s.targetSteps>0 && s.completedFE<=s.targetSteps && s.dt>0 && s.controlEveryFESteps>0 && s.controlPhase==s.completedFE%s.controlEveryFESteps);
}
}
QByteArray RStructuralPiezoLayeredRestartCodec::encode(const RStructuralPiezoLayeredRunState &state) {
    auto copy=state;Archive a;
    // Bounded capacity hint; record() retains every authoritative format guard.
    qsizetype capacity=4096;
    for(const auto *blob:{&state.configurationIdentity,&state.opaqueControllerState,&state.lastFrameBytes})
        capacity+=std::min(blob->size(),qsizetype(bound)-capacity);
    a.bytes.reserve(capacity);
    record(a,copy);auto digest=RHashSha256::hash(a.bytes);
    require(uint(a.bytes.size()+digest.size())<=bound);a.bytes.append(digest);return a.bytes;
}
bool RStructuralPiezoLayeredRestartCodec::canonicalEquals(const RStructuralPiezoLayeredRunState &expected,
                                                            const RStructuralPiezoLayeredRunState &actual) {
    const QByteArray actualBytes=encode(actual);
    require(actualBytes.size()>=32 && uint(actualBytes.size())<=bound);
    const uint actualBodySize=uint(actualBytes.size())-32u;
    auto copy=expected;Archive a;a.comparing=true;a.comparison=&actualBytes;a.comparisonSize=actualBodySize;
    record(a,copy);
    require(a.offset<=bound-32u);
    return a.equal && a.offset==actualBodySize;
}
RStructuralPiezoLayeredRunState RStructuralPiezoLayeredRestartCodec::decode(const QByteArray &bytes) {
    require(bytes.size()>=32 && bytes.size()<=bound);const auto payload=bytes.left(bytes.size()-32);
    require(RHashSha256::hash(payload)==bytes.right(32));
    Archive a;a.reading=true;a.bytes=payload;RStructuralPiezoLayeredRunState s;record(a,s);require(a.offset==uint(payload.size()));return s;
}
