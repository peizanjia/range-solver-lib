#ifndef RSTRUCTURAL_PIEZO_LAYERED_RUN_STATE_H
#define RSTRUCTURAL_PIEZO_LAYERED_RUN_STATE_H
#include "rstructural_piezo_layered_dynamics.h"

struct RStructuralPiezoLayeredVoltageEnvelope
{
    double minimum=0,maximum=0,maximumRate=0;
};
struct RStructuralPiezoLayeredVoltageAction
{
    QString portId;
    double requested=0,amplitudeLimited=0,previousApplied=0,applied=0;
    bool amplitudeClipped=false,rateClipped=false;
};
//! Complete accepted runtime. Electrical channels remain a named port table in physical.
struct RStructuralPiezoLayeredRunState
{
    RStructuralPiezoLayeredDynamicsState physical;
    QString sourceFingerprint,forcingIdentity,controllerIdentity,backendIdentity;
    QByteArray configurationIdentity,opaqueControllerState,lastFrameBytes;
    uint targetSteps=0,completedFE=0,controlEveryFESteps=1,controlTick=0,controlPhase=0;
    double initialTime=0,dt=0,nominalControlTs=0;
    std::vector<RStructuralPiezoLayeredVoltageEnvelope> envelopes;
    std::vector<double> heldRequested,lastApplied;
    std::vector<RStructuralPiezoLayeredVoltageAction> actions;
};
#endif
