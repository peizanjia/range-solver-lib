#ifndef RSTRUCTURAL_PIEZO_AXIAL_RUN_STATE_H
#define RSTRUCTURAL_PIEZO_AXIAL_RUN_STATE_H
#include "rstructural_piezo_native_readout.h"

struct RStructuralPiezoAxialVoltageEnvelope
{
    double minimum=0,maximum=0,maximumRate=0;
};
struct RStructuralPiezoAxialVoltageAction
{
    QString portId;
    double requested=0,amplitudeLimited=0,previousApplied=0,applied=0;
    bool amplitudeClipped=false,rateClipped=false;
};
//! Complete accepted state. Compatibility hashes are not authentication.
struct RStructuralPiezoAxialRunState
{
    RStructuralTrussDynamicsState physical,initialPhysicalReference;
    std::vector<RStructuralPiezoNativeReadoutSample> readouts;
    QString sourceFingerprint,forcingIdentity,controllerIdentity,backendIdentity;
    QByteArray configurationIdentity,opaqueControllerState,lastFrameBytes;
    uint targetSteps=0,completedFE=0,controlEveryFESteps=1,controlTick=0,controlPhase=0;
    double initialTime=0,dt=0,nominalControlTs=0;
    std::vector<RStructuralPiezoAxialVoltageEnvelope> envelopes;
    std::vector<double> heldRequested,lastApplied;
    std::vector<RStructuralPiezoAxialVoltageAction> actions;
};
#endif
