#include "rstructural_piezo_axial_restart_codec.h"
#include <rbl_hash_sha256.h>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>
#include <stdexcept>

namespace {
constexpr uint bound=32*1024*1024, vectorBound=1024*1024;
void require(bool b) {if(!b)throw std::invalid_argument("Invalid axial canonical restart bytes");}
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
        if(reading){require(n<=(uint(bytes.size())-offset)/8);v.resize(n);}
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
    template<class V,class F> void entries(V &v,F f) {require(v.size()<=65536);uint n=uint(v.size());count(n,65536);if(reading){require(n<=(uint(bytes.size())-offset)/8);v.resize(n);}for(auto &x:v)f(x);}
};
void physical(Archive &a,RStructuralTrussDynamicsState &s) {
    a.real(s.time);a.blob(s.systemIdentity);
    for(auto *v:{&s.displacements,&s.velocities,&s.accelerations,&s.reactions})a.vector(*v);
    a.entries(s.appliedInputs.memberInputs,[&](RStructuralTrussInput &input){
        uint kind=uint(input.index());a.count(kind,3);
        double v=0;
        if(!a.reading)std::visit([&](const auto &x){using T=std::decay_t<decltype(x)>;if constexpr(!std::is_same_v<T,RStructuralTrussNoInput>)v=x.value;},input);
        if(kind)a.real(v);
        if(a.reading){switch(kind){case 0:input=RStructuralTrussNoInput{};break;case 1:input=RStructuralTrussVoltage{v};break;case 2:input=RStructuralTrussCharge{v};break;case 3:input=RStructuralTrussCableLength{v};break;}}
    });
    a.vector(s.appliedInputs.externalForces);a.vector(s.appliedInputs.rightExternalForces);
    for(double *x:{&s.storedEnergy,&s.kineticEnergy,&s.cumulativeExternalWork,&s.cumulativeSupportWork,&s.cumulativeElectricalWork,&s.cumulativeLengthWork,&s.initialTotalEnergy,&s.cumulativeDampingDissipation})a.real(*x);
    require(s.appliedInputs.rightExternalForces.empty());
}
void record(Archive &a,RStructuralPiezoAxialRunState &s) {
    QString magic="range-fea/axial-native-restart/binary64/v1";a.text(magic);require(magic=="range-fea/axial-native-restart/binary64/v1");
    a.text(s.sourceFingerprint);a.text(s.forcingIdentity);a.text(s.controllerIdentity);a.text(s.backendIdentity);
    a.blob(s.configurationIdentity);a.blob(s.opaqueControllerState);a.blob(s.lastFrameBytes);
    for(uint *x:{&s.targetSteps,&s.completedFE,&s.controlEveryFESteps,&s.controlTick,&s.controlPhase})a.count(*x);
    a.real(s.initialTime);a.real(s.dt);a.real(s.nominalControlTs);
    a.entries(s.envelopes,[&](RStructuralPiezoAxialVoltageEnvelope &e){a.real(e.minimum);a.real(e.maximum);a.real(e.maximumRate);require(e.minimum<=e.maximum && e.maximumRate>=0);});
    a.vector(s.heldRequested);a.vector(s.lastApplied);
    a.entries(s.actions,[&](RStructuralPiezoAxialVoltageAction &v){a.text(v.portId);for(double *x:{&v.requested,&v.amplitudeLimited,&v.previousApplied,&v.applied})a.real(*x);a.boolean(v.amplitudeClipped);a.boolean(v.rateClipped);});
    physical(a,s.physical);physical(a,s.initialPhysicalReference);
    a.entries(s.readouts,[&](RStructuralPiezoNativeReadoutSample &p){a.text(p.portId);a.enumeration(p.mode,4);a.text(p.unit);a.boolean(p.actionRecord);for(double *x:{&p.voltage,&p.charge,&p.deltaVoltage,&p.deltaCharge,&p.value,&p.deltaValue})a.real(*x);});
    require(s.targetSteps>0 && s.completedFE<=s.targetSteps && s.dt>0 && s.controlEveryFESteps>0 && s.controlPhase==s.completedFE%s.controlEveryFESteps);
}
}
QByteArray RStructuralPiezoAxialRestartCodec::encode(const RStructuralPiezoAxialRunState &state) {
    auto copy=state;Archive a;
    // Bounded capacity hint; record() retains every authoritative format guard.
    qsizetype capacity=4096;
    for(const auto *blob:{&state.configurationIdentity,&state.opaqueControllerState,&state.lastFrameBytes})
        capacity+=std::min(blob->size(),qsizetype(bound)-capacity);
    a.bytes.reserve(capacity);
    record(a,copy);auto digest=RHashSha256::hash(a.bytes);
    require(uint(a.bytes.size()+digest.size())<=bound);a.bytes.append(digest);return a.bytes;
}
bool RStructuralPiezoAxialRestartCodec::canonicalEquals(const RStructuralPiezoAxialRunState &expected,
                                                        const RStructuralPiezoAxialRunState &actual) {
    const QByteArray actualBytes=encode(actual);
    require(actualBytes.size()>=32 && uint(actualBytes.size())<=bound);
    const uint actualBodySize=uint(actualBytes.size())-32u;
    auto copy=expected;Archive a;a.comparing=true;a.comparison=&actualBytes;a.comparisonSize=actualBodySize;
    record(a,copy);
    require(a.offset<=bound-32u);
    return a.equal && a.offset==actualBodySize;
}
RStructuralPiezoAxialRunState RStructuralPiezoAxialRestartCodec::decode(const QByteArray &bytes) {
    require(bytes.size()>=32 && bytes.size()<=bound);const auto payload=bytes.left(bytes.size()-32);
    require(RHashSha256::hash(payload)==bytes.right(32));
    Archive a;a.reading=true;a.bytes=payload;RStructuralPiezoAxialRunState s;record(a,s);require(a.offset==uint(payload.size()));return s;
}
