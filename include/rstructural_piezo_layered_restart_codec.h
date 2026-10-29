#ifndef RSTRUCTURAL_PIEZO_LAYERED_RESTART_CODEC_H
#define RSTRUCTURAL_PIEZO_LAYERED_RESTART_CODEC_H
#include "rstructural_piezo_layered_run_state.h"
//! Canonical bounded little-endian IEEE binary64 codec. Integrity is not authentication.
//! Decode returns a new object; driver and controller semantic recertification belong to the facade.
class RStructuralPiezoLayeredRestartCodec
{
public:
    static QByteArray encode(const RStructuralPiezoLayeredRunState &state);
    static bool canonicalEquals(const RStructuralPiezoLayeredRunState &expected,
                                const RStructuralPiezoLayeredRunState &actual);
    static RStructuralPiezoLayeredRunState decode(const QByteArray &bytes);
};
#endif
