#ifndef RSTRUCTURAL_PIEZO_AXIAL_RESTART_CODEC_H
#define RSTRUCTURAL_PIEZO_AXIAL_RESTART_CODEC_H
#include "rstructural_piezo_axial_run_state.h"
//! Bounded canonical little-endian binary64 archive. Semantic restore belongs to the facade.
class RStructuralPiezoAxialRestartCodec
{
public:
    static QByteArray encode(const RStructuralPiezoAxialRunState &state);
    static bool canonicalEquals(const RStructuralPiezoAxialRunState &expected,
                                const RStructuralPiezoAxialRunState &actual);
    static RStructuralPiezoAxialRunState decode(const QByteArray &bytes);
};
#endif
