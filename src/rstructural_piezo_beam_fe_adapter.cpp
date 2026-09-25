#include "rstructural_piezo_beam_fe_adapter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include <QByteArray>
#include <QCryptographicHash>

#include <rbl_error.h>
#include <rbl_rvector.h>

#include <rml_boundary_condition.h>
#include <rml_element.h>
#include <rml_environment_condition.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_material_property.h>
#include <rml_matrix_solver_conf.h>
#include <rml_model.h>
#include <rml_point.h>
#include <rml_sparse_matrix.h>
#include <rml_variable.h>

#include "rmatrixsolver.h"

namespace
{
constexpr double referenceTemperatureKelvin = 293.15;
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();

[[noreturn]] void fail(const QString &message)
{
    const QByteArray utf8 = message.toUtf8();
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",utf8.constData());
}

double stableNorm(const std::vector<double> &values,const QString &description)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            fail(description+QStringLiteral(" contains a non-finite value."));
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        fail(description+QStringLiteral(" norm is non-finite."));
    return norm;
}

double gammaBound(std::size_t operationCount)
{
    const double scaledCount = double(operationCount)*std::numeric_limits<double>::epsilon();
    if (!std::isfinite(scaledCount) || scaledCount >= 1.0)
        fail(QStringLiteral("Beam piezo roundoff operation count is invalid."));
    return scaledCount/(1.0-scaledCount);
}

double roundoffBound(std::size_t operationCount,double magnitude)
{
    if (!std::isfinite(magnitude) || magnitude < 0.0)
        fail(QStringLiteral("Beam piezo roundoff scale is invalid."));
    return gammaBound(operationCount)*magnitude+
        double(operationCount)*std::numeric_limits<double>::denorm_min();
}

void appendUInt(QByteArray &bytes,uint value)
{
    for (unsigned int byte=0;byte<4;byte++)
        bytes.append(char((value>>(8*byte))&0xffu));
}

void appendSize(QByteArray &bytes,std::size_t value)
{
    static_assert(sizeof(std::size_t)<=sizeof(quint64),"Fingerprint size exceeds 64 bits.");
    const quint64 encoded = quint64(value);
    for (unsigned int byte=0;byte<8;byte++)
        bytes.append(char((encoded>>(8*byte))&0xffu));
}

void appendDouble(QByteArray &bytes,double value)
{
    static_assert(sizeof(double)==sizeof(quint64),"Fingerprint requires binary64 doubles.");
    quint64 bits = 0;
    std::memcpy(&bits,&value,sizeof(bits));
    for (unsigned int byte=0;byte<8;byte++)
        bytes.append(char((bits>>(8*byte))&0xffu));
}

struct BeamMember
{
    uint elementId = 0;
    uint node0 = 0;
    uint node1 = 0;
    double length = 0.0;
    double area = 0.0;
    double secondMomentZ = 0.0;
    double youngModulus = 0.0;
    double density = 0.0;
    double gravityX = 0.0;
    double gravityY = 0.0;
};

struct PatchData
{
    RStructuralPiezoBeamPatch patch;
    std::size_t memberIndex = 0;
    std::array<std::array<double,6>,6> stiffness{};
    std::array<double,6> coupling{};
    double capacitance = 0.0;
};

struct ValidatedModel
{
    std::vector<BeamMember> members;
    std::vector<PatchData> patches;
    std::vector<double> prescribed;
    std::vector<bool> prescribedMask;
    std::vector<double> mechanicalLoads;
    std::vector<double> mechanicalLoadMagnitude;
    std::vector<std::size_t> mechanicalLoadTerms;
    QString fingerprint;
};

void accumulateMagnitude(std::vector<double> &magnitudes,
                         std::vector<std::size_t> &termCounts,
                         std::size_t position,
                         double term,
                         const QString &description)
{
    const double magnitude = std::fabs(term);
    const double updated = magnitudes[position]+magnitude;
    if (!std::isfinite(magnitude) || !std::isfinite(updated) ||
        termCounts[position]==std::numeric_limits<std::size_t>::max())
        fail(description+QStringLiteral(" scale or operation count overflowed."));
    magnitudes[position] = updated;
    ++termCounts[position];
}

void addMechanicalLoad(ValidatedModel &data,std::size_t dof,double value)
{
    const double updated = data.mechanicalLoads[dof]+value;
    if (!std::isfinite(value) || !std::isfinite(updated))
        fail(QStringLiteral("BEAM1 piezo load assembly overflowed."));
    data.mechanicalLoads[dof] = updated;
    accumulateMagnitude(data.mechanicalLoadMagnitude,data.mechanicalLoadTerms,dof,value,
                        QStringLiteral("BEAM1 piezo load"));
}

void addSymmetric(std::array<std::array<double,6>,6> &matrix,
                  std::size_t row,std::size_t column,double value)
{
    matrix[row][column] += value;
    if (row != column)
        matrix[column][row] += value;
}

std::array<std::size_t,6> localDofs(const BeamMember &member)
{
    return {{3*std::size_t(member.node0),
             3*std::size_t(member.node0)+1,
             3*std::size_t(member.node0)+2,
             3*std::size_t(member.node1),
             3*std::size_t(member.node1)+1,
             3*std::size_t(member.node1)+2}};
}

std::array<std::array<double,6>,6> hostStiffness(const BeamMember &member)
{
    std::array<std::array<double,6>,6> stiffness{};
    const double axial = member.youngModulus*member.area/member.length;
    addSymmetric(stiffness,0,0,axial);
    addSymmetric(stiffness,3,3,axial);
    addSymmetric(stiffness,0,3,-axial);

    const double length = member.length;
    const double factor = member.youngModulus*member.secondMomentZ/(length*length*length);
    const std::array<std::size_t,4> dofs = {{1,2,4,5}};
    const double terms[4][4] = {
        {12.0,6.0*length,-12.0,6.0*length},
        {6.0*length,4.0*length*length,-6.0*length,2.0*length*length},
        {-12.0,-6.0*length,12.0,-6.0*length},
        {6.0*length,2.0*length*length,-6.0*length,4.0*length*length}
    };
    for (std::size_t row=0;row<4;row++)
        for (std::size_t column=row;column<4;column++)
            addSymmetric(stiffness,dofs[row],dofs[column],factor*terms[row][column]);
    return stiffness;
}

void buildPatchConstitutive(PatchData &data,const BeamMember &member)
{
    constexpr double inverseSqrtThree = 0.57735026918962576451;
    const double points[2] = {0.5*(1.0-inverseSqrtThree),
                              0.5*(1.0+inverseSqrtThree)};
    const double integrationWeight = 0.5*member.length;
    for (double xi : points)
    {
        std::array<double,6> strainRow{};
        strainRow[0] = -1.0/member.length;
        strainRow[3] = 1.0/member.length;
        const double curvature[4] = {
            (-6.0+12.0*xi)/(member.length*member.length),
            (-4.0+6.0*xi)/member.length,
            (6.0-12.0*xi)/(member.length*member.length),
            (-2.0+6.0*xi)/member.length
        };
        const std::size_t bendingDofs[4] = {1,2,4,5};
        for (std::size_t local=0;local<4;local++)
            strainRow[bendingDofs[local]] -= data.patch.centroidOffsetY*curvature[local];

        const double mechanicalScale = data.patch.youngModulus*data.patch.area*integrationWeight;
        const double couplingScale = data.patch.coupling*data.patch.area/
            data.patch.thickness*integrationWeight;
        for (std::size_t row=0;row<6;row++)
        {
            data.coupling[row] += couplingScale*strainRow[row];
            for (std::size_t column=0;column<6;column++)
                data.stiffness[row][column] += mechanicalScale*
                    strainRow[row]*strainRow[column];
        }
    }

    data.capacitance = data.patch.permittivity*data.patch.area*member.length/
        (data.patch.thickness*data.patch.thickness);
    if (!std::isfinite(data.capacitance) || data.capacitance <= 0.0)
        fail(QStringLiteral("BEAM1 piezo patch %1 has invalid capacitance.").arg(data.patch.elementId));
}

void appendDoubleVector(QByteArray &bytes,const std::vector<double> &values)
{
    appendSize(bytes,values.size());
    for (double value : values) appendDouble(bytes,value);
}

QString resultCertificate(const RStructuralPiezoBeamFEResult &result)
{
    QByteArray bytes("range-fea/piezo/beam-result/v1");
    appendDoubleVector(bytes,result.displacements);
    appendDoubleVector(bytes,result.rotations);
    appendDoubleVector(bytes,result.internalForces);
    appendDoubleVector(bytes,result.patchVoltages);
    appendDoubleVector(bytes,result.patchCharges);
    appendDouble(bytes,result.enthalpy);
    appendDouble(bytes,result.translationalResidualNorm);
    appendDouble(bytes,result.rotationalResidualNorm);
    appendDouble(bytes,result.translationalRoundoffBound);
    appendDouble(bytes,result.rotationalRoundoffBound);
    appendDouble(bytes,result.electricResidualNorm);
    appendDouble(bytes,result.electricRoundoffBound);
    const QByteArray fingerprint = result.modelFingerprint.toUtf8();
    appendSize(bytes,std::size_t(fingerprint.size()));
    bytes.append(fingerprint);
    return QString::fromLatin1(
        QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex());
}

void setPrescribed(ValidatedModel &data,std::size_t dof,double value)
{
    if (!std::isfinite(value))
        fail(QStringLiteral("BEAM1 piezo mechanical boundary value is non-finite."));
    if (data.prescribedMask[dof])
    {
        const double tolerance = 1.0e-12*std::max({1.0,std::fabs(value),
                                                   std::fabs(data.prescribed[dof])});
        if (std::fabs(data.prescribed[dof]-value)>tolerance)
            fail(QStringLiteral("BEAM1 piezo model contains contradictory point constraints."));
        return;
    }
    data.prescribedMask[dof] = true;
    data.prescribed[dof] = value;
}

template<typename Condition>
double componentValue(const Condition &condition,RVariableType type,double time)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod || !condition.getComponent(position).getEnabled())
        return 0.0;
    const RConditionComponent &component = condition.getComponent(position);
    if (component.empty())
        fail(QStringLiteral("BEAM1 piezo boundary component '%1' has no values.")
                 .arg(RVariable::getName(type)));
    return component.get(time);
}

ValidatedModel inspectModel(const RModel &model,
                            const RStructuralPiezoBeamFEInput &input)
{
    if (model.getNNodes()<2 || model.getNNodes()>std::numeric_limits<uint>::max()/3)
        fail(QStringLiteral("BEAM1 piezo adapter requires a representable model with at least two nodes."));
    if (input.electricalMode != RStructuralPiezoBeamElectricalMode::PrescribedVoltage &&
        input.electricalMode != RStructuralPiezoBeamElectricalMode::PrescribedCharge)
        fail(QStringLiteral("BEAM1 piezo electrical mode is unsupported."));
    if (input.patches.empty())
        fail(QStringLiteral("BEAM1 piezo adapter requires at least one bonded patch."));
    if (input.patches.size()>std::numeric_limits<uint>::max())
        fail(QStringLiteral("BEAM1 piezo patch count exceeds the sparse solver index limit."));
    const double tolerances[6] = {input.translationalAbsoluteTolerance,
                                  input.translationalRelativeTolerance,
                                  input.rotationalAbsoluteTolerance,
                                  input.rotationalRelativeTolerance,
                                  input.electricAbsoluteTolerance,
                                  input.electricRelativeTolerance};
    for (double tolerance : tolerances)
        if (!std::isfinite(tolerance) || tolerance<0.0)
            fail(QStringLiteral("BEAM1 piezo residual tolerances must be finite and nonnegative."));
    if (input.translationalAbsoluteTolerance==0.0 && input.translationalRelativeTolerance==0.0)
        fail(QStringLiteral("BEAM1 piezo translational residual gate cannot be disabled."));
    if (input.rotationalAbsoluteTolerance==0.0 && input.rotationalRelativeTolerance==0.0)
        fail(QStringLiteral("BEAM1 piezo rotational residual gate cannot be disabled."));
    if (input.electricAbsoluteTolerance==0.0 && input.electricRelativeTolerance==0.0)
        fail(QStringLiteral("BEAM1 piezo electrical residual gate cannot be disabled."));

    std::map<uint,std::size_t> memberPosition;
    std::vector<bool> lineMemberSeen(model.getNElements(),false);
    ValidatedModel data;
    for (uint lineId=0;lineId<model.getNLines();lineId++)
    {
        const RLine &line = model.getLine(lineId);
        if (line.empty()) continue;
        const double area = line.getCrossArea();
        const double iz = line.getSecondMomentZ();
        const RR3Vector &localY = line.getLocalYDirection();
        const double localYNorm = localY.length();
        if (!std::isfinite(area) || area<=0.0 || !std::isfinite(iz) || iz<=0.0 ||
            !std::isfinite(localYNorm) || localYNorm<=0.0 || localY[1]<=0.0 ||
            std::fabs(localY[0])>1.0e-12*localYNorm ||
            std::fabs(localY[2])>1.0e-12*localYNorm)
            fail(QStringLiteral("BEAM1 piezo currently requires positive SI A/Iz and local Y aligned with global +Y on line '%1'.")
                     .arg(line.getName()));

        const RMaterial &material = line.getMaterial();
        const uint modulusPosition = material.findPosition(RMaterialProperty::ModulusOfElasiticity);
        if (modulusPosition>=material.size())
            fail(QStringLiteral("BEAM1 piezo line '%1' has no Young modulus material property.")
                     .arg(line.getName()));
        const RMaterialProperty &modulusProperty = material.get(modulusPosition);
        if (modulusProperty.empty())
            fail(QStringLiteral("BEAM1 piezo line '%1' has an empty Young modulus table.")
                     .arg(line.getName()));
        const double modulus = modulusProperty.get(referenceTemperatureKelvin);
        if (!std::isfinite(modulus) || modulus<=0.0)
            fail(QStringLiteral("BEAM1 piezo line '%1' has an invalid Young modulus at 293.15 K.")
                     .arg(line.getName()));
        double density = 0.0;
        const uint densityPosition = material.findPosition(RMaterialProperty::Density);
        if (densityPosition<material.size())
        {
            const RMaterialProperty &densityProperty = material.get(densityPosition);
            if (densityProperty.empty())
                fail(QStringLiteral("BEAM1 piezo line '%1' has an empty density table.")
                         .arg(line.getName()));
            density = densityProperty.get(referenceTemperatureKelvin);
            if (!std::isfinite(density) || density<=0.0)
                fail(QStringLiteral("BEAM1 piezo line '%1' has an invalid density at 293.15 K.")
                         .arg(line.getName()));
        }
        const REnvironmentCondition &gravity =
            line.getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        const RVariableType gravityTypes[3] = {R_VARIABLE_G_ACCELERATION_X,
                                                R_VARIABLE_G_ACCELERATION_Y,
                                                R_VARIABLE_G_ACCELERATION_Z};
        const double gravityX = componentValue(gravity,gravityTypes[0],
                                                model.getTimeSolver().getCurrentTime());
        const double gravityY = componentValue(gravity,gravityTypes[1],
                                                model.getTimeSolver().getCurrentTime());
        const double gravityZ = componentValue(gravity,gravityTypes[2],
                                                model.getTimeSolver().getCurrentTime());
        if (!std::isfinite(gravityX) || !std::isfinite(gravityY) || !std::isfinite(gravityZ) ||
            gravityZ!=0.0)
            fail(QStringLiteral("BEAM1 piezo supports only finite in-plane line gravity; set out-of-plane gravity to zero."));
        if ((gravityX!=0.0 || gravityY!=0.0) && density<=0.0)
            fail(QStringLiteral("BEAM1 piezo in-plane gravity requires a positive line-material density."));

        for (uint position=0;position<line.size();position++)
        {
            const uint elementId = line.get(position);
            if (elementId>=model.getNElements() || lineMemberSeen[elementId])
                fail(QStringLiteral("BEAM1 piezo line groups contain an invalid or repeated element ID."));
            const RElement &element = model.getElement(elementId);
            if (element.getType()!=R_ELEMENT_BEAM1 || element.size()!=2)
                fail(QStringLiteral("BEAM1 piezo adapter accepts only native two-node BEAM1 line elements."));
            lineMemberSeen[elementId] = true;
            const uint node0 = element.getNodeId(0);
            const uint node1 = element.getNodeId(1);
            if (node0>=model.getNNodes() || node1>=model.getNNodes() || node0==node1)
                fail(QStringLiteral("BEAM1 piezo element %1 has invalid endpoint IDs.").arg(elementId));
            const RNode &first = model.getNode(node0);
            const RNode &second = model.getNode(node1);
            const double coordinates[6] = {first.getX(),first.getY(),first.getZ(),
                                           second.getX(),second.getY(),second.getZ()};
            for (double coordinate : coordinates)
                if (!std::isfinite(coordinate))
                    fail(QStringLiteral("BEAM1 piezo element %1 has non-finite coordinates.").arg(elementId));
            const double dx = second.getX()-first.getX();
            const double dy = second.getY()-first.getY();
            const double dz = second.getZ()-first.getZ();
            const double length = std::hypot(std::hypot(dx,dy),dz);
            const double geometryTolerance = 1.0e-12*std::max(1.0,length);
            if (!std::isfinite(length) || length<=0.0 || dx<=0.0 ||
                std::fabs(dy)>geometryTolerance || std::fabs(dz)>geometryTolerance ||
                std::fabs(first.getY())>geometryTolerance || std::fabs(second.getY())>geometryTolerance ||
                std::fabs(first.getZ())>geometryTolerance || std::fabs(second.getZ())>geometryTolerance)
                fail(QStringLiteral("BEAM1 piezo element %1 is not a straight +X member on the global XY reference line.")
                         .arg(elementId));
            BeamMember member;
            member.elementId = elementId;
            member.node0 = node0;
            member.node1 = node1;
            member.length = length;
            member.area = area;
            member.secondMomentZ = iz;
            member.youngModulus = modulus;
            member.density = density;
            member.gravityX = gravityX;
            member.gravityY = gravityY;
            memberPosition.emplace(elementId,data.members.size());
            data.members.push_back(member);
        }
    }
    if (data.members.empty() || model.getNElements(R_ELEMENT_BEAM1)!=data.members.size())
        fail(QStringLiteral("BEAM1 piezo requires every native BEAM1 element to belong to a line group."));
    for (uint elementId=0;elementId<model.getNElements();elementId++)
    {
        const RElementType type = model.getElement(elementId).getType();
        if (type!=R_ELEMENT_BEAM1 && type!=R_ELEMENT_POINT)
            fail(QStringLiteral("BEAM1 piezo does not support additional %1 elements.")
                     .arg(RElement::getName(type)));
        if (type==R_ELEMENT_BEAM1 && !lineMemberSeen[elementId])
            fail(QStringLiteral("BEAM1 piezo element %1 has no line-group section data.").arg(elementId));
    }

    std::sort(data.members.begin(),data.members.end(),[&model](const BeamMember &first,const BeamMember &second)
    {
        return model.getNode(first.node0).getX()<model.getNode(second.node0).getX();
    });
    memberPosition.clear();
    std::vector<uint> incidence(model.getNNodes(),0);
    for (std::size_t position=0;position<data.members.size();position++)
    {
        const BeamMember &member = data.members[position];
        memberPosition.emplace(member.elementId,position);
        ++incidence[member.node0];
        ++incidence[member.node1];
        if (position>0 && data.members[position-1].node1!=member.node0)
            fail(QStringLiteral("BEAM1 piezo requires one connected, non-branching beam chain ordered from +X."));
    }
    if (incidence[data.members.front().node0]!=1 || incidence[data.members.back().node1]!=1)
        fail(QStringLiteral("BEAM1 piezo model must have exactly two chain endpoints."));
    for (std::size_t position=1;position+1<data.members.size();position++)
        if (incidence[data.members[position].node0]!=2 || incidence[data.members[position].node1]!=2)
            fail(QStringLiteral("BEAM1 piezo model contains a disconnected or branched beam chain."));
    for (uint node=0;node<model.getNNodes();node++)
        if (incidence[node]==0)
            fail(QStringLiteral("BEAM1 piezo model contains unused node %1.").arg(node));

    data.prescribed.assign(3*std::size_t(model.getNNodes()),0.0);
    data.prescribedMask.assign(data.prescribed.size(),false);
    data.mechanicalLoads.assign(data.prescribed.size(),0.0);
    data.mechanicalLoadMagnitude.assign(data.prescribed.size(),0.0);
    data.mechanicalLoadTerms.assign(data.prescribed.size(),0);
    const double time = model.getTimeSolver().getCurrentTime();
    if (!std::isfinite(time))
        fail(QStringLiteral("BEAM1 piezo model time is non-finite."));
    const RVariableType displacementTypes[3] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    const RVariableType rotationTypes[3] = {
        R_VARIABLE_ROTATION_X,R_VARIABLE_ROTATION_Y,R_VARIABLE_ROTATION_Z};
    const RVariableType forceTypes[3] = {R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    const RVariableType momentTypes[3] = {R_VARIABLE_MOMENT_X,R_VARIABLE_MOMENT_Y,R_VARIABLE_MOMENT_Z};
    for (uint pointId=0;pointId<model.getNPoints();pointId++)
    {
        const RPoint &point = model.getPoint(pointId);
        for (uint pointElementPosition=0;pointElementPosition<point.size();pointElementPosition++)
        {
            const RElement &pointElement = model.getElement(point.get(pointElementPosition));
            if (pointElement.getType()!=R_ELEMENT_POINT || pointElement.size()!=1)
                fail(QStringLiteral("BEAM1 piezo point groups must contain one-node POINT elements."));
            const uint node = pointElement.getNodeId(0);
            if (node>=model.getNNodes() || incidence[node]==0)
                fail(QStringLiteral("BEAM1 piezo point condition references a node outside the beam chain."));
            for (uint bcIndex=0;bcIndex<point.getNBoundaryConditions();bcIndex++)
            {
                const RBoundaryCondition &condition = point.getBoundaryCondition(bcIndex);
                if (!(RBoundaryCondition::getProblemTypeMask(condition.getType())&R_PROBLEM_STRESS)) continue;
                if (condition.getType()==R_BOUNDARY_CONDITION_DISPLACEMENT)
                {
                    for (uint component=0;component<3;component++)
                    {
                        const uint position = condition.findComponentPosition(displacementTypes[component]);
                        if (position==RConstants::eod || !condition.getComponent(position).getEnabled()) continue;
                        const double value = componentValue(condition,displacementTypes[component],time);
                        if (!std::isfinite(value))
                            fail(QStringLiteral("BEAM1 piezo displacement condition is non-finite."));
                        if (component<2)
                            setPrescribed(data,3*std::size_t(node)+component,value);
                        else if (value!=0.0)
                            fail(QStringLiteral("BEAM1 piezo rejects nonzero out-of-plane displacement constraints."));
                    }
                }
                else if (condition.getType()==R_BOUNDARY_CONDITION_ROTATION)
                {
                    for (uint component=0;component<3;component++)
                    {
                        const uint position = condition.findComponentPosition(rotationTypes[component]);
                        if (position==RConstants::eod || !condition.getComponent(position).getEnabled()) continue;
                        const double value = componentValue(condition,rotationTypes[component],time);
                        if (!std::isfinite(value))
                            fail(QStringLiteral("BEAM1 piezo rotation condition is non-finite."));
                        if (component==2)
                            setPrescribed(data,3*std::size_t(node)+2,value);
                        else if (value!=0.0)
                            fail(QStringLiteral("BEAM1 piezo rejects nonzero out-of-plane rotation constraints."));
                    }
                }
                else if (condition.getType()==R_BOUNDARY_CONDITION_FORCE ||
                         condition.getType()==R_BOUNDARY_CONDITION_MOMENT)
                {
                    const bool force = condition.getType()==R_BOUNDARY_CONDITION_FORCE;
                    const RVariableType *types = force ? forceTypes : momentTypes;
                    for (uint component=0;component<3;component++)
                    {
                        const double value = componentValue(condition,types[component],time);
                        if (!std::isfinite(value))
                            fail(QStringLiteral("BEAM1 piezo force or moment condition is non-finite."));
                        const bool active = force ? component<2 : component==2;
                        if (!active && value!=0.0)
                            fail(QStringLiteral("BEAM1 piezo rejects nonzero out-of-plane force or moment."));
                        if (active)
                            addMechanicalLoad(data,
                                3*std::size_t(node)+(force ? component : 2),value);
                    }
                }
                else
                {
                    fail(QStringLiteral("BEAM1 piezo direct adapter supports only point displacement, rotation, force and moment conditions."));
                }
            }
        }
    }
    for (uint lineId=0;lineId<model.getNLines();lineId++)
    {
        const RLine &line = model.getLine(lineId);
        for (uint bcIndex=0;bcIndex<line.getNBoundaryConditions();bcIndex++)
        {
            if (RBoundaryCondition::getProblemTypeMask(
                    line.getBoundaryCondition(bcIndex).getType())&R_PROBLEM_STRESS)
                fail(QStringLiteral("BEAM1 piezo direct adapter accepts mechanical conditions only on point groups."));
        }
    }

    // Consistent Euler-Bernoulli nodal loads for constant in-plane gravity.
    for (const BeamMember &member : data.members)
    {
        const double axialLoadPerLength = member.density*member.area*member.gravityX;
        const double transverseLoadPerLength = member.density*member.area*member.gravityY;
        const double halfLength = 0.5*member.length;
        addMechanicalLoad(data,3*std::size_t(member.node0),axialLoadPerLength*halfLength);
        addMechanicalLoad(data,3*std::size_t(member.node1),axialLoadPerLength*halfLength);
        addMechanicalLoad(data,3*std::size_t(member.node0)+1,transverseLoadPerLength*halfLength);
        addMechanicalLoad(data,3*std::size_t(member.node1)+1,transverseLoadPerLength*halfLength);
        const double consistentMoment = transverseLoadPerLength*member.length*member.length/12.0;
        addMechanicalLoad(data,3*std::size_t(member.node0)+2,consistentMoment);
        addMechanicalLoad(data,3*std::size_t(member.node1)+2,-consistentMoment);
    }

    std::set<uint> patchIds;
    data.patches.reserve(input.patches.size());
    for (const RStructuralPiezoBeamPatch &patch : input.patches)
    {
        const auto found = memberPosition.find(patch.elementId);
        if (found==memberPosition.end() || !patchIds.insert(patch.elementId).second)
            fail(QStringLiteral("BEAM1 piezo patch IDs must identify unique native beam elements."));
        if (!std::isfinite(patch.area) || patch.area<=0.0 ||
            !std::isfinite(patch.thickness) || patch.thickness<=0.0 ||
            !std::isfinite(patch.centroidOffsetY) ||
            !std::isfinite(patch.youngModulus) || patch.youngModulus<=0.0 ||
            !std::isfinite(patch.coupling) || patch.coupling==0.0 ||
            !std::isfinite(patch.permittivity) || patch.permittivity<=0.0 ||
            !std::isfinite(patch.prescribedValue))
            fail(QStringLiteral("BEAM1 piezo patch %1 has invalid SI area, thickness, offset, material or boundary data.")
                     .arg(patch.elementId));
        PatchData patchData;
        patchData.patch = patch;
        patchData.memberIndex = found->second;
        const BeamMember &member = data.members[found->second];
        if (member.gravityX!=0.0 || member.gravityY!=0.0)
            fail(QStringLiteral("BEAM1 piezo rejects gravity on patched elements because patch mass loading is not configured."));
        buildPatchConstitutive(patchData,member);
        data.patches.push_back(patchData);
    }

    QByteArray fingerprint("range-fea/piezo/beam-fe/direct/v2");
    appendUInt(fingerprint,uint(model.getNNodes()));
    appendUInt(fingerprint,uint(data.members.size()));
    for (uint node=0;node<model.getNNodes();node++)
    {
        appendDouble(fingerprint,model.getNode(node).getX());
        appendDouble(fingerprint,model.getNode(node).getY());
        appendDouble(fingerprint,model.getNode(node).getZ());
    }
    for (const BeamMember &member : data.members)
    {
        appendUInt(fingerprint,member.elementId);
        appendUInt(fingerprint,member.node0);
        appendUInt(fingerprint,member.node1);
        appendDouble(fingerprint,member.length);
        appendDouble(fingerprint,member.area);
        appendDouble(fingerprint,member.secondMomentZ);
        appendDouble(fingerprint,member.youngModulus);
        appendDouble(fingerprint,member.density);
        appendDouble(fingerprint,member.gravityX);
        appendDouble(fingerprint,member.gravityY);
    }
    for (std::size_t dof=0;dof<data.prescribed.size();dof++)
    {
        appendUInt(fingerprint,data.prescribedMask[dof] ? 1u : 0u);
        appendDouble(fingerprint,data.prescribed[dof]);
        appendDouble(fingerprint,data.mechanicalLoads[dof]);
        appendDouble(fingerprint,data.mechanicalLoadMagnitude[dof]);
        appendSize(fingerprint,data.mechanicalLoadTerms[dof]);
    }
    appendUInt(fingerprint,static_cast<uint>(input.electricalMode));
    appendDouble(fingerprint,input.translationalAbsoluteTolerance);
    appendDouble(fingerprint,input.translationalRelativeTolerance);
    appendDouble(fingerprint,input.rotationalAbsoluteTolerance);
    appendDouble(fingerprint,input.rotationalRelativeTolerance);
    appendDouble(fingerprint,input.electricAbsoluteTolerance);
    appendDouble(fingerprint,input.electricRelativeTolerance);
    appendUInt(fingerprint,uint(data.patches.size()));
    for (const PatchData &patch : data.patches)
    {
        appendUInt(fingerprint,patch.patch.elementId);
        appendDouble(fingerprint,patch.patch.area);
        appendDouble(fingerprint,patch.patch.thickness);
        appendDouble(fingerprint,patch.patch.centroidOffsetY);
        appendDouble(fingerprint,patch.patch.youngModulus);
        appendDouble(fingerprint,patch.patch.coupling);
        appendDouble(fingerprint,patch.patch.permittivity);
        appendDouble(fingerprint,patch.patch.prescribedValue);
    }
    data.fingerprint = QString::fromLatin1(
        QCryptographicHash::hash(fingerprint,QCryptographicHash::Sha256).toHex());
    return data;
}

void addLocalMatrix(RSparseMatrix &matrix,
                    const std::array<std::size_t,6> &dofs,
                    const std::array<std::array<double,6>,6> &local)
{
    for (std::size_t row=0;row<6;row++)
        for (std::size_t column=0;column<6;column++)
        {
            const double value = local[row][column];
            if (!std::isfinite(value))
                fail(QStringLiteral("BEAM1 piezo tangent contains a non-finite value."));
            if (value!=0.0)
                matrix.addValue(uint(dofs[row]),uint(dofs[column]),value);
        }
}

void accumulateLocalProducts(std::vector<double> &magnitudes,
                             std::vector<std::size_t> &termCounts,
                             const std::array<std::size_t,6> &dofs,
                             const std::array<std::array<double,6>,6> &local,
                             const std::vector<double> &solution)
{
    for (std::size_t row=0;row<6;row++)
        for (std::size_t column=0;column<6;column++)
        {
            const double coefficient = local[row][column];
            if (coefficient==0.0) continue;
            const double term = coefficient*solution[dofs[column]];
            if (!std::isfinite(term))
                fail(QStringLiteral("BEAM1 piezo local matrix product is non-finite."));
            accumulateMagnitude(magnitudes,termCounts,dofs[row],term,
                                QStringLiteral("BEAM1 piezo matrix-product"));
        }
}

std::size_t roundoffOperationCount(std::size_t matrixTerms,
                                   std::size_t loadTerms,
                                   std::size_t sparseProductTerms)
{
    const std::size_t max = std::numeric_limits<std::size_t>::max();
    if (matrixTerms>(max-8)/4 || loadTerms>(max-8)/2 || sparseProductTerms>(max-8)/2)
        fail(QStringLiteral("BEAM1 piezo roundoff operation count overflowed."));
    const std::size_t matrixOperations = 4*matrixTerms;
    const std::size_t loadOperations = 2*loadTerms;
    const std::size_t productOperations = 2*sparseProductTerms;
    if (matrixOperations>max-loadOperations-8 ||
        matrixOperations+loadOperations+8>max-productOperations)
        fail(QStringLiteral("BEAM1 piezo roundoff operation count overflowed."));
    return matrixOperations+loadOperations+productOperations+8;
}

std::vector<double> multiply(const RSparseMatrix &matrix,const std::vector<double> &values)
{
    if (matrix.getNRows()!=values.size())
        fail(QStringLiteral("BEAM1 piezo sparse matrix/vector dimensions do not match."));
    RRVector vector(uint(values.size()),0.0);
    for (std::size_t i=0;i<values.size();i++) vector[uint(i)] = values[i];
    RRVector product;
    RSparseMatrix::mlt(matrix,vector,product);
    if (product.getNRows()!=values.size())
        fail(QStringLiteral("BEAM1 piezo matrix product has an invalid dimension."));
    return std::vector<double>(product.begin(),product.end());
}

RVariable makeOutputVariable(const RModel &model,RVariableType type,
                             RVariableApplyType applyType,uint vectors,uint values,
                             const std::vector<double> &data,const QString &units)
{
    if (data.size()!=std::size_t(vectors)*values)
        fail(QStringLiteral("BEAM1 piezo output variable has inconsistent dimensions."));
    RVariable variable;
    const uint existing = model.findVariable(type);
    if (existing==RConstants::eod)
        variable = RVariable(type,applyType);
    else
        variable = model.getVariable(existing);
    variable.setType(type);
    variable.setApplyType(applyType);
    variable.setUnits(units);
    variable.resize(vectors,values);
    for (uint vector=0;vector<vectors;vector++)
        for (uint value=0;value<values;value++)
            variable.setValue(vector,value,data[std::size_t(vector)*values+value]);
    return variable;
}

}

void RStructuralPiezoBeamFEAdapter::validate(
    const RModel &model,const RStructuralPiezoBeamFEInput &input)
{
    (void)inspectModel(model,input);
}

RStructuralPiezoBeamFEResult RStructuralPiezoBeamFEAdapter::solve(
    const RModel &model,const RStructuralPiezoBeamFEInput &input)
{
    const ValidatedModel data = inspectModel(model,input);
    const std::size_t nNodes = model.getNNodes();
    const std::size_t mechanicalDofs = 3*nNodes;
    const bool chargeMode = input.electricalMode==
        RStructuralPiezoBeamElectricalMode::PrescribedCharge;
    const std::size_t totalDofs = mechanicalDofs+(chargeMode ? data.patches.size() : 0);
    if (totalDofs>std::numeric_limits<uint>::max())
        fail(QStringLiteral("BEAM1 piezo system exceeds the sparse solver index limit."));

    RSparseMatrix tangent;
    tangent.setNRows(uint(totalDofs));
    std::vector<double> rhs(totalDofs,0.0);
    for (std::size_t dof=0;dof<mechanicalDofs;dof++) rhs[dof] = data.mechanicalLoads[dof];

    for (const BeamMember &member : data.members)
    {
        const auto dofs = localDofs(member);
        addLocalMatrix(tangent,dofs,hostStiffness(member));
    }
    for (std::size_t patchIndex=0;patchIndex<data.patches.size();patchIndex++)
    {
        const PatchData &patch = data.patches[patchIndex];
        const BeamMember &member = data.members[patch.memberIndex];
        const auto dofs = localDofs(member);
        addLocalMatrix(tangent,dofs,patch.stiffness);
        if (chargeMode)
        {
            const uint voltageDof = uint(mechanicalDofs+patchIndex);
            for (std::size_t local=0;local<6;local++)
            {
                const uint mechanicalDof = uint(dofs[local]);
                const double value = patch.coupling[local];
                if (value!=0.0)
                {
                    tangent.addValue(mechanicalDof,voltageDof,value);
                    tangent.addValue(voltageDof,mechanicalDof,value);
                }
            }
            tangent.addValue(voltageDof,voltageDof,-patch.capacitance);
            rhs[mechanicalDofs+patchIndex] = -patch.patch.prescribedValue;
        }
        else
        {
            for (std::size_t local=0;local<6;local++)
                rhs[dofs[local]] -= patch.coupling[local]*patch.patch.prescribedValue;
        }
    }

    std::vector<double> prescribed(totalDofs,0.0);
    std::vector<bool> prescribedMask(totalDofs,false);
    std::copy(data.prescribed.begin(),data.prescribed.end(),prescribed.begin());
    std::copy(data.prescribedMask.begin(),data.prescribedMask.end(),prescribedMask.begin());
    std::vector<std::size_t> globalToFree(totalDofs,noFreeDof);
    std::vector<std::size_t> freeToGlobal;
    for (std::size_t dof=0;dof<totalDofs;dof++)
    {
        if (!prescribedMask[dof])
        {
            globalToFree[dof] = freeToGlobal.size();
            freeToGlobal.push_back(dof);
        }
    }
    if (freeToGlobal.size()>std::numeric_limits<uint>::max())
        fail(QStringLiteral("BEAM1 piezo free system exceeds the sparse solver index limit."));

    RSparseMatrix freeTangent;
    freeTangent.setNRows(uint(freeToGlobal.size()));
    std::vector<double> freeRhs(freeToGlobal.size(),0.0);
    for (std::size_t freeRow=0;freeRow<freeToGlobal.size();freeRow++)
    {
        const std::size_t globalRow = freeToGlobal[freeRow];
        freeRhs[freeRow] = rhs[globalRow];
        const RSparseVector<double> &entries = tangent.getVector(uint(globalRow));
        for (uint position=0;position<entries.size();position++)
        {
            const std::size_t globalColumn = entries.getIndex(position);
            const double value = entries.getValue(position);
            if (!std::isfinite(value))
                fail(QStringLiteral("BEAM1 piezo assembled tangent is non-finite."));
            if (prescribedMask[globalColumn])
                freeRhs[freeRow] -= value*prescribed[globalColumn];
            else
            {
                freeTangent.addValue(uint(freeRow),uint(globalToFree[globalColumn]),value);
            }
        }
        if (!std::isfinite(freeRhs[freeRow]))
            fail(QStringLiteral("BEAM1 piezo prescribed-value elimination overflowed."));
    }

    std::vector<double> fullSolution(prescribed);
    if (!freeToGlobal.empty())
    {
        std::vector<double> diagonalScale(freeToGlobal.size(),0.0);
        RSparseMatrix scaledTangent;
        scaledTangent.setNRows(uint(freeToGlobal.size()));
        RRVector scaledRhs(uint(freeToGlobal.size()),0.0);
        for (std::size_t row=0;row<freeToGlobal.size();row++)
        {
            const double diagonal = freeTangent.findValue(uint(row),uint(row));
            if (!std::isfinite(diagonal) || diagonal==0.0)
                fail(QStringLiteral("BEAM1 piezo free tangent has a zero or non-finite diagonal."));
            diagonalScale[row] = 1.0/std::sqrt(std::fabs(diagonal));
            if (!std::isfinite(diagonalScale[row]) || diagonalScale[row]<=0.0)
                fail(QStringLiteral("BEAM1 piezo symmetric scaling is invalid."));
            scaledRhs[uint(row)] = diagonalScale[row]*freeRhs[row];
            if (!std::isfinite(scaledRhs[uint(row)]))
                fail(QStringLiteral("BEAM1 piezo scaled right-hand side is non-finite."));
        }
        for (std::size_t row=0;row<freeToGlobal.size();row++)
        {
            const RSparseVector<double> &entries = freeTangent.getVector(uint(row));
            for (uint position=0;position<entries.size();position++)
            {
                const uint column = entries.getIndex(position);
                const double value = entries.getValue(position)*diagonalScale[row]*diagonalScale[column];
                if (!std::isfinite(value))
                    fail(QStringLiteral("BEAM1 piezo scaled tangent is non-finite."));
                if (value!=0.0)
                    scaledTangent.addValue(uint(row),column,value);
            }
        }

        RMatrixSolverConf solverConf = model.getMatrixSolverConf(RMatrixSolverConf::GMRES);
        if (solverConf.getNInnerIterations()==0 || solverConf.getNOuterIterations()==0 ||
            !std::isfinite(solverConf.getSolverCvgValue()) || solverConf.getSolverCvgValue()<=0.0)
            fail(QStringLiteral("BEAM1 piezo GMRES configuration is invalid."));
        solverConf.setNInnerIterations(unsigned(std::min<std::size_t>(
            solverConf.getNInnerIterations(),freeToGlobal.size())));
        RMatrixSolver solver(solverConf);
        solver.disableConvergenceLogFile();
        RRVector scaledSolution;
        solver.solve(scaledTangent,scaledRhs,scaledSolution,R_MATRIX_PRECONDITIONER_NONE,1);
        if (scaledSolution.getNRows()!=freeToGlobal.size())
            fail(QStringLiteral("BEAM1 piezo GMRES returned an invalid solution dimension."));
        for (std::size_t free=0;free<freeToGlobal.size();free++)
        {
            const double value = diagonalScale[free]*scaledSolution[uint(free)];
            if (!std::isfinite(value))
                fail(QStringLiteral("BEAM1 piezo GMRES returned a non-finite solution."));
            fullSolution[freeToGlobal[free]] = value;
        }
    }

    const std::vector<double> freshProduct = multiply(tangent,fullSolution);
    std::vector<double> translationalResidual;
    std::vector<double> rotationalResidual;
    std::vector<double> freeTranslationalLoads;
    std::vector<double> freeRotationalLoads;
    std::vector<double> freeTranslationalInternal;
    std::vector<double> freeRotationalInternal;
    std::vector<double> translationalRoundoffRows;
    std::vector<double> rotationalRoundoffRows;
    std::vector<double> electricResidual;
    std::vector<double> prescribedCharges;
    std::vector<double> recoveredCharges;
    std::vector<double> electricRoundoffRows;
    std::vector<double> voltageForce(mechanicalDofs,0.0);
    std::vector<double> voltageForceMagnitude(mechanicalDofs,0.0);
    std::vector<std::size_t> voltageForceTerms(mechanicalDofs,0);
    std::vector<double> matrixProductMagnitude(mechanicalDofs,0.0);
    std::vector<std::size_t> matrixProductTerms(mechanicalDofs,0);
    for (const BeamMember &member : data.members)
        accumulateLocalProducts(matrixProductMagnitude,matrixProductTerms,
                                localDofs(member),hostStiffness(member),fullSolution);
    for (const PatchData &patch : data.patches)
    {
        const BeamMember &member = data.members[patch.memberIndex];
        const auto dofs = localDofs(member);
        accumulateLocalProducts(matrixProductMagnitude,matrixProductTerms,
                                dofs,patch.stiffness,fullSolution);
    }
    if (!chargeMode)
    {
        for (const PatchData &patch : data.patches)
        {
            const BeamMember &member = data.members[patch.memberIndex];
            const auto dofs = localDofs(member);
            for (std::size_t local=0;local<6;local++)
            {
                const double term = patch.coupling[local]*patch.patch.prescribedValue;
                const double updated = voltageForce[dofs[local]]+term;
                if (!std::isfinite(term) || !std::isfinite(updated))
                    fail(QStringLiteral("BEAM1 piezo voltage-load assembly overflowed."));
                voltageForce[dofs[local]] = updated;
                accumulateMagnitude(voltageForceMagnitude,voltageForceTerms,
                                   dofs[local],term,QStringLiteral("BEAM1 piezo voltage load"));
            }
        }
    }
    else
    {
        for (std::size_t patchIndex=0;patchIndex<data.patches.size();patchIndex++)
        {
            const PatchData &patch = data.patches[patchIndex];
            const BeamMember &member = data.members[patch.memberIndex];
            const auto dofs = localDofs(member);
            const double voltage = fullSolution[mechanicalDofs+patchIndex];
            for (std::size_t local=0;local<6;local++)
            {
                const double term = patch.coupling[local]*voltage;
                if (!std::isfinite(term))
                    fail(QStringLiteral("BEAM1 piezo electromechanical product is non-finite."));
                accumulateMagnitude(matrixProductMagnitude,matrixProductTerms,
                                    dofs[local],term,
                                    QStringLiteral("BEAM1 piezo electromechanical matrix"));
            }
        }
    }
    for (std::size_t dof=0;dof<mechanicalDofs;dof++)
    {
        if (!prescribedMask[dof])
        {
            const double residual = freshProduct[dof]-rhs[dof];
            const double internal = freshProduct[dof]+voltageForce[dof];
            if (!std::isfinite(residual) || !std::isfinite(internal))
                fail(QStringLiteral("BEAM1 piezo fresh mechanical residual is non-finite."));
            const std::size_t sparseTerms = tangent.getVector(uint(dof)).size();
            if (data.mechanicalLoadTerms[dof]>
                std::numeric_limits<std::size_t>::max()-voltageForceTerms[dof])
                fail(QStringLiteral("BEAM1 piezo load roundoff operation count overflowed."));
            const std::size_t loadTerms = data.mechanicalLoadTerms[dof]+voltageForceTerms[dof];
            const double rowMagnitude = matrixProductMagnitude[dof]+
                data.mechanicalLoadMagnitude[dof]+voltageForceMagnitude[dof];
            const double rowBound = roundoffBound(
                roundoffOperationCount(matrixProductTerms[dof],loadTerms,sparseTerms),
                rowMagnitude);
            if (dof%3==2)
            {
                rotationalResidual.push_back(residual);
                freeRotationalLoads.push_back(rhs[dof]);
                freeRotationalInternal.push_back(internal);
                rotationalRoundoffRows.push_back(rowBound);
            }
            else
            {
                translationalResidual.push_back(residual);
                freeTranslationalLoads.push_back(rhs[dof]);
                freeTranslationalInternal.push_back(internal);
                translationalRoundoffRows.push_back(rowBound);
            }
        }
    }
    for (std::size_t patchIndex=0;patchIndex<data.patches.size();patchIndex++)
    {
        const PatchData &patch = data.patches[patchIndex];
        double voltage = patch.patch.prescribedValue;
        if (chargeMode) voltage = fullSolution[mechanicalDofs+patchIndex];
        const BeamMember &member = data.members[patch.memberIndex];
        const auto dofs = localDofs(member);
        double couplingWork = 0.0;
        double absoluteCouplingWork = 0.0;
        std::size_t couplingTerms = 0;
        for (std::size_t local=0;local<6;local++)
        {
            const double term = patch.coupling[local]*fullSolution[dofs[local]];
            const double updated = couplingWork+term;
            const double updatedMagnitude = absoluteCouplingWork+std::fabs(term);
            if (!std::isfinite(term) || !std::isfinite(updated) ||
                !std::isfinite(updatedMagnitude) ||
                (patch.coupling[local]!=0.0 && couplingTerms==
                    std::numeric_limits<std::size_t>::max()))
                fail(QStringLiteral("BEAM1 piezo patch coupling work overflowed."));
            couplingWork = updated;
            absoluteCouplingWork = updatedMagnitude;
            if (patch.coupling[local]!=0.0) ++couplingTerms;
        }
        const double capacitanceTerm = patch.capacitance*voltage;
        const double charge = capacitanceTerm-couplingWork;
        if (!std::isfinite(charge) || !std::isfinite(voltage))
            fail(QStringLiteral("BEAM1 piezo recovered voltage or charge is non-finite."));
        if (chargeMode)
        {
            const double residual = charge-patch.patch.prescribedValue;
            const double magnitude = std::fabs(capacitanceTerm)+absoluteCouplingWork+
                std::fabs(patch.patch.prescribedValue);
            if (!std::isfinite(residual) || !std::isfinite(magnitude))
                fail(QStringLiteral("BEAM1 piezo fresh electric residual is non-finite."));
            electricResidual.push_back(residual);
            electricRoundoffRows.push_back(roundoffBound(
                roundoffOperationCount(couplingTerms,1,0),magnitude));
            prescribedCharges.push_back(patch.patch.prescribedValue);
            recoveredCharges.push_back(charge);
        }
    }

    RStructuralPiezoBeamFEResult result;
    result.displacements.assign(3*nNodes,0.0);
    result.rotations.assign(3*nNodes,0.0);
    result.internalForces.assign(6*nNodes,0.0);
    result.patchVoltages.resize(data.patches.size(),0.0);
    result.patchCharges.resize(data.patches.size(),0.0);
    for (std::size_t node=0;node<nNodes;node++)
    {
        result.displacements[3*node] = fullSolution[3*node];
        result.displacements[3*node+1] = fullSolution[3*node+1];
        result.rotations[3*node+2] = fullSolution[3*node+2];
    }
    for (std::size_t dof=0;dof<mechanicalDofs;dof++)
    {
        const std::size_t node = dof/3;
        const std::size_t component = dof%3;
        const std::size_t resultComponent = component==2 ? 5 : component;
        result.internalForces[6*node+resultComponent] =
            freshProduct[dof]+voltageForce[dof];
    }
    result.translationalResidualNorm = stableNorm(
        translationalResidual,QStringLiteral("BEAM1 piezo fresh translational residual"));
    result.rotationalResidualNorm = stableNorm(
        rotationalResidual,QStringLiteral("BEAM1 piezo fresh rotational residual"));
    result.translationalRoundoffBound = stableNorm(
        translationalRoundoffRows,QStringLiteral("BEAM1 piezo translational roundoff bound"));
    result.rotationalRoundoffBound = stableNorm(
        rotationalRoundoffRows,QStringLiteral("BEAM1 piezo rotational roundoff bound"));
    result.electricResidualNorm = stableNorm(electricResidual,
                                             QStringLiteral("BEAM1 piezo fresh electric residual"));
    result.electricRoundoffBound = stableNorm(
        electricRoundoffRows,QStringLiteral("BEAM1 piezo electric roundoff bound"));
    const double translationalScale = std::max(
        stableNorm(freeTranslationalLoads,QStringLiteral("BEAM1 piezo free translational loads")),
        stableNorm(freeTranslationalInternal,QStringLiteral("BEAM1 piezo free translational internal forces")));
    const double rotationalScale = std::max(
        stableNorm(freeRotationalLoads,QStringLiteral("BEAM1 piezo free rotational loads")),
        stableNorm(freeRotationalInternal,QStringLiteral("BEAM1 piezo free rotational internal moments")));
    const double translationalTolerance = input.translationalAbsoluteTolerance+
        input.translationalRelativeTolerance*translationalScale+
        result.translationalRoundoffBound;
    const double rotationalTolerance = input.rotationalAbsoluteTolerance+
        input.rotationalRelativeTolerance*rotationalScale+
        result.rotationalRoundoffBound;
    const double electricScale = std::max(stableNorm(prescribedCharges,
                                                      QStringLiteral("BEAM1 piezo prescribed charges")),
                                           stableNorm(recoveredCharges,
                                                      QStringLiteral("BEAM1 piezo recovered charges")));
    const double electricTolerance = input.electricAbsoluteTolerance+
        input.electricRelativeTolerance*electricScale+
        result.electricRoundoffBound;
    if (!std::isfinite(translationalTolerance) || !std::isfinite(rotationalTolerance) ||
        !std::isfinite(electricTolerance) ||
        result.translationalResidualNorm>translationalTolerance ||
        result.rotationalResidualNorm>rotationalTolerance ||
        (chargeMode && result.electricResidualNorm>electricTolerance))
        fail(QStringLiteral("BEAM1 piezo solution failed its fresh N/N*m/C residual certification."));

    double enthalpy = 0.0;
    for (const BeamMember &member : data.members)
    {
        const auto dofs = localDofs(member);
        const auto base = hostStiffness(member);
        for (std::size_t row=0;row<6;row++)
            for (std::size_t column=0;column<6;column++)
                enthalpy += 0.5*fullSolution[dofs[row]]*base[row][column]*fullSolution[dofs[column]];
    }
    for (std::size_t patchIndex=0;patchIndex<data.patches.size();patchIndex++)
    {
        const PatchData &patch = data.patches[patchIndex];
        const BeamMember &member = data.members[patch.memberIndex];
        const auto dofs = localDofs(member);
        const double voltage = chargeMode
            ? fullSolution[mechanicalDofs+patchIndex] : patch.patch.prescribedValue;
        double quadratic = 0.0;
        double couplingWork = 0.0;
        for (std::size_t row=0;row<6;row++)
        {
            couplingWork += patch.coupling[row]*fullSolution[dofs[row]];
            for (std::size_t column=0;column<6;column++)
                quadratic += 0.5*fullSolution[dofs[row]]*patch.stiffness[row][column]*
                    fullSolution[dofs[column]];
        }
        enthalpy += quadratic+voltage*couplingWork-0.5*patch.capacitance*voltage*voltage;
        result.patchVoltages[patchIndex] = voltage;
        result.patchCharges[patchIndex] = patch.capacitance*voltage-couplingWork;
    }
    if (!std::isfinite(enthalpy))
        fail(QStringLiteral("BEAM1 piezo enthalpy is non-finite."));
    result.enthalpy = enthalpy;
    result.modelFingerprint = data.fingerprint;
    result.resultSeal = resultCertificate(result);
    return result;
}

void RStructuralPiezoBeamFEAdapter::publish(
    RModel &model,const RStructuralPiezoBeamFEInput &input,
    const RStructuralPiezoBeamFEResult &result)
{
    const ValidatedModel data = inspectModel(model,input);
    const std::size_t nNodes = model.getNNodes();
    if (result.displacements.size()!=3*nNodes || result.rotations.size()!=3*nNodes ||
        result.internalForces.size()!=6*nNodes ||
        result.patchVoltages.size()!=data.patches.size() ||
        result.patchCharges.size()!=data.patches.size() ||
        result.modelFingerprint!=data.fingerprint ||
        !std::isfinite(result.enthalpy) ||
        !std::isfinite(result.translationalResidualNorm) || result.translationalResidualNorm<0.0 ||
        !std::isfinite(result.rotationalResidualNorm) || result.rotationalResidualNorm<0.0 ||
        !std::isfinite(result.translationalRoundoffBound) || result.translationalRoundoffBound<0.0 ||
        !std::isfinite(result.rotationalRoundoffBound) || result.rotationalRoundoffBound<0.0 ||
        !std::isfinite(result.electricResidualNorm) || result.electricResidualNorm<0.0 ||
        !std::isfinite(result.electricRoundoffBound) || result.electricRoundoffBound<0.0)
        fail(QStringLiteral("BEAM1 piezo candidate dimensions or model binding are invalid."));
    for (double value : result.displacements)
        if (!std::isfinite(value)) fail(QStringLiteral("BEAM1 piezo displacement is non-finite."));
    for (double value : result.rotations)
        if (!std::isfinite(value)) fail(QStringLiteral("BEAM1 piezo rotation is non-finite."));
    for (double value : result.internalForces)
        if (!std::isfinite(value)) fail(QStringLiteral("BEAM1 piezo internal force is non-finite."));
    for (double value : result.patchVoltages)
        if (!std::isfinite(value)) fail(QStringLiteral("BEAM1 piezo voltage is non-finite."));
    for (double value : result.patchCharges)
        if (!std::isfinite(value)) fail(QStringLiteral("BEAM1 piezo charge is non-finite."));

    if (result.resultSeal.isEmpty() || result.resultSeal!=resultCertificate(result))
        fail(QStringLiteral("BEAM1 piezo candidate differs from its trusted solve certificate."));

    std::vector<double> displacementOutput(3*nNodes,0.0);
    std::vector<double> rotationOutput(3*nNodes,0.0);
    for (std::size_t node=0;node<nNodes;node++)
        for (std::size_t component=0;component<3;component++)
        {
            displacementOutput[component*nNodes+node] = result.displacements[3*node+component];
            rotationOutput[component*nNodes+node] = result.rotations[3*node+component];
        }
    std::vector<double> potentialOutput(model.getNElements(),0.0);
    std::vector<double> chargeOutput(model.getNElements(),0.0);
    for (std::size_t patch=0;patch<data.patches.size();patch++)
    {
        const uint element = data.patches[patch].patch.elementId;
        potentialOutput[element] = result.patchVoltages[patch];
        chargeOutput[element] = result.patchCharges[patch];
    }

    // Construct every output before atomically replacing the model result set.
    const std::vector<RVariable> outputs = {
        makeOutputVariable(model,R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3,uint(nNodes),
                           displacementOutput,QStringLiteral("m")),
        makeOutputVariable(model,R_VARIABLE_ROTATION,R_VARIABLE_APPLY_NODE,3,uint(nNodes),
                           rotationOutput,QStringLiteral("rad")),
        makeOutputVariable(model,R_VARIABLE_ELECTRIC_POTENTIAL,R_VARIABLE_APPLY_ELEMENT,1,
                           model.getNElements(),potentialOutput,
                           QStringLiteral("V (patch voltage difference)")),
        makeOutputVariable(model,R_VARIABLE_ELECTRIC_CHARGE,R_VARIABLE_APPLY_ELEMENT,1,
                           model.getNElements(),chargeOutput,QStringLiteral("C"))
    };
    model.upsertVariablesAtomic(outputs);
}
