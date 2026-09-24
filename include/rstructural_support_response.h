#ifndef RSTRUCTURAL_SUPPORT_RESPONSE_H
#define RSTRUCTURAL_SUPPORT_RESPONSE_H

#include <QString>

#include <rml_structural_interaction.h>

//! Accepted-step support-on-structure response in physical units.
//! A nonempty motionId is exactly the stable ID of a prescribed motion.
//! Empty motionId denotes a stationary global point support.
struct RStructuralSupportSample
{
    QString motionId;
    RStructuralDofRef target;
    double time = 0.0;
    double reaction = 0.0;
    double velocity = 0.0;
    double power = 0.0;
    double cumulativeWork = 0.0;
};

#endif // RSTRUCTURAL_SUPPORT_RESPONSE_H
