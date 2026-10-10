#include "rstructural_piezo_model_adapter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <QByteArray>
#include <QCryptographicHash>

#include <rbl_error.h>
#include <rbl_rvector.h>

#include <rml_matrix_solver_conf.h>
#include <rml_line.h>
#include <rml_material_property.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_piezo_setup.h>
#include <rml_variable.h>

#include "rmatrixpreconditioner.h"
#include "rmatrixsolver.h"
#include "rstructural_piezo_coupled_1d.h"
#include "rstructural_piezo_kernel.h"
#include "rstructural_static_truss_model.h"

namespace
{
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();

[[noreturn]] void fail(const QString &message)
{
    const QByteArray utf8 = message.toUtf8();
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",utf8.constData());
}

double stableNorm(const std::vector<double> &values)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            fail("Piezo static solve encountered a non-finite vector value.");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        fail("Piezo static solve vector norm is non-finite.");
    return norm;
}

// Count-delimited little-endian integers and exact binary64 values preserve
// endpoint order and prevent concatenation ambiguities. This matches the beam
// adapter's local provenance convention; it is not a cryptographic authority.
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
    static_assert(sizeof(double)==sizeof(quint64) && std::numeric_limits<double>::is_iec559,
                  "Fingerprint requires IEEE binary64 doubles.");
    quint64 bits = 0;
    std::memcpy(&bits,&value,sizeof(bits));
    for (unsigned int byte=0;byte<8;byte++)
        bytes.append(char((bits>>(8*byte))&0xffu));
}

void appendString(QByteArray &bytes,const QString &value)
{
    const QByteArray encoded = value.toUtf8();
    appendSize(bytes,std::size_t(encoded.size()));
    bytes.append(encoded);
}

void appendDoubleVector(QByteArray &bytes,const std::vector<double> &values)
{
    appendSize(bytes,values.size());
    for (double value : values) appendDouble(bytes,value);
}

void appendBoolVector(QByteArray &bytes,const std::vector<bool> &values)
{
    appendSize(bytes,values.size());
    for (bool value : values) appendUInt(bytes,value ? 1u : 0u);
}

QString digest(const QByteArray &bytes)
{
    return QString::fromLatin1(QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex());
}

QString resultCertificate(const RStructuralPiezoModelResult &result,
                          const QString &inputFingerprint)
{
    QByteArray bytes("range-fea/piezo/truss-result/v1");
    appendString(bytes,inputFingerprint);
    appendDoubleVector(bytes,result.displacements);
    appendDoubleVector(bytes,result.internalForces);
    appendDoubleVector(bytes,result.memberAxialForces);
    appendDoubleVector(bytes,result.memberRightCharges);
    appendBoolVector(bytes,result.memberElectricalPortAvailable);
    appendDouble(bytes,result.residualNorm);
    appendDouble(bytes,result.enthalpy);
    appendUInt(bytes,static_cast<uint>(result.electricalMode));
    appendDoubleVector(bytes,result.nodalPotentials);
    appendDoubleVector(bytes,result.nodalConductorCharges);
    appendDouble(bytes,result.electricResidualNorm);
    return digest(bytes);
}

struct PiezoModelData
{
    RStructuralStaticTrussModelData truss;
    //! Active port pointers in sorted member order; nullptr means purely mechanical.
    std::vector<const RStructuralPiezoPort*> portsByMember;
};

struct CoupledPiezoModelData
{
    RStructuralPiezoCoupled1DInput input;
    std::array<double,3> axis = {{0.0,0.0,0.0}};
};

PiezoModelData buildPiezoModelData(const RModel &model)
{
    const RProblemSetup &problemSetup = model.getProblemSetup();
    const RStructuralPiezoSetup &setup = problemSetup.getStructuralPiezoSetup();
    if (!setup.enabled)
        fail("Piezo static adapter requires an enabled structural piezo setup.");
    if (problemSetup.getStructuralCableSetup().enabled)
        fail("Structural cable and piezo configurations cannot be enabled together.");

    setup.validate(model);
    PiezoModelData data;
    const bool thermalPiezo = model.getProblemTaskTree().getProblemTypeMask() ==
        (R_PROBLEM_HEAT | R_PROBLEM_STRESS);
    if (thermalPiezo && setup.electricalMode != RStructuralPiezoElectricalMode::NodalPotential)
        fail("Heat->Piezo supports only the nodal-potential coupled rod mode.");
    data.truss = RStructuralStaticTrussModel::build(model,false,thermalPiezo);
    if (!setup.passiveElementIds.empty() &&
        setup.electricalMode != RStructuralPiezoElectricalMode::LegacyTerminalVoltage)
        fail("Mixed active/passive trusses require LegacyTerminalVoltage mode.");
    if (setup.ports.size()+setup.passiveElementIds.size() != data.truss.members.size())
        fail("Piezo setup must classify every TRUSS1 member with one terminal port or as passive.");

    std::map<uint,const RStructuralPiezoPort*> portByElement;
    for (const RStructuralPiezoPort &port : setup.ports)
    {
        if (!std::isfinite(port.coupling) || !std::isfinite(port.permittivity) ||
            port.permittivity <= 0.0 || !std::isfinite(port.voltage))
        {
            fail(QString("Piezo port '%1' has invalid coupling, permittivity or voltage.").arg(port.id));
        }
        if (!portByElement.emplace(port.elementId,&port).second)
            fail("Piezo setup has more than one port for a TRUSS1 member.");
    }

    std::set<uint> passiveElements;
    for (uint elementId : setup.passiveElementIds)
    {
        if (!passiveElements.insert(elementId).second || portByElement.count(elementId))
            fail("Piezo active/passive member classification is duplicated or overlaps.");
        if (elementId >= model.getNElements() ||
            model.getElement(elementId).getType() != R_ELEMENT_TRUSS1)
            fail("Passive piezo classification requires a valid TRUSS1 element.");
    }
    data.portsByMember.reserve(data.truss.members.size());
    for (const RStructuralStaticTrussMember &member : data.truss.members)
    {
        const auto found = portByElement.find(member.elementId);
        if (found == portByElement.end() && !passiveElements.count(member.elementId))
            fail(QString("TRUSS1 element %1 has no explicit active/passive classification.").arg(member.elementId));
        data.portsByMember.push_back(found == portByElement.end() ? nullptr : found->second);
    }
    return data;
}

double dot(const std::array<double,3> &first,const std::array<double,3> &second)
{
    return first[0]*second[0]+first[1]*second[1]+first[2]*second[2];
}

std::array<double,3> subtract(const std::array<double,3> &first,
                              const std::array<double,3> &second)
{
    return {{first[0]-second[0],first[1]-second[1],first[2]-second[2]}};
}

double vectorNorm(const std::array<double,3> &value)
{
    return std::hypot(std::hypot(value[0],value[1]),value[2]);
}

double integratedThermalStrain(const RMaterialProperty &alpha,double currentKelvin)
{
    constexpr double referenceKelvin = 293.15;
    if (alpha.empty() || !std::isfinite(currentKelvin) || currentKelvin < 0.0)
        fail("Heat->Piezo requires a finite nodal absolute temperature and a nonempty alpha table.");
    const double low = std::min(referenceKelvin,currentKelvin);
    const double high = std::max(referenceKelvin,currentKelvin);
    double previous = low;
    double integral = 0.0;
    for (uint point=0;point<alpha.size();point++)
    {
        const double knot = alpha.getKey(point);
        if (!std::isfinite(knot)) fail("Heat->Piezo alpha table has a non-finite temperature key.");
        if (knot <= previous || knot >= high) continue;
        const double left = alpha.get(previous);
        const double right = alpha.get(knot);
        if (!std::isfinite(left) || !std::isfinite(right))
            fail("Heat->Piezo alpha table has a non-finite coefficient.");
        integral += 0.5*(left+right)*(knot-previous);
        previous = knot;
    }
    const double left = alpha.get(previous);
    const double right = alpha.get(high);
    if (!std::isfinite(left) || !std::isfinite(right))
        fail("Heat->Piezo alpha table has a non-finite coefficient.");
    integral += 0.5*(left+right)*(high-previous);
    const double strain = (currentKelvin >= referenceKelvin ? 1.0 : -1.0)*integral;
    if (!std::isfinite(strain)) fail("Heat->Piezo integrated thermal strain is non-finite.");
    return strain;
}

CoupledPiezoModelData buildCoupledPiezoModelData(const RModel &model,
                                                const PiezoModelData &data,
                                                bool useCurrentThermalField = false)
{
    if (data.truss.nodes.empty() || data.truss.members.empty())
        fail("Nodal-potential piezo mode requires at least one node and TRUSS1 member.");

    CoupledPiezoModelData coupled;
    const RStructuralStaticTrussMember &axisMember = data.truss.members.front();
    const auto &axisStart = data.truss.nodes[axisMember.node1];
    const auto &axisEnd = data.truss.nodes[axisMember.node2];
    for (std::size_t component=0;component<3;component++)
        coupled.axis[component] = (axisEnd[component]-axisStart[component])/axisMember.referenceLength;

    const std::size_t nNodes = data.truss.nodes.size();
    RStructuralPiezoCoupled1DInput &input = coupled.input;
    input.coordinates.resize(nNodes,0.0);
    input.mechanicalLoads.assign(nNodes,0.0);
    input.prescribedConductorCharges.assign(nNodes,0.0);
    input.displacementPrescribed.assign(nNodes,false);
    input.prescribedDisplacements.assign(nNodes,0.0);
    input.potentialPrescribed.assign(nNodes,false);
    input.prescribedPotentials.assign(nNodes,0.0);

    const RVariable *temperature = nullptr;
    if (useCurrentThermalField &&
        model.getProblemTaskTree().getProblemTypeMask() == (R_PROBLEM_HEAT | R_PROBLEM_STRESS))
    {
        const uint position = model.findVariable(R_VARIABLE_TEMPERATURE);
        if (position == RConstants::eod)
            fail("Heat->Piezo requires the converged nodal heat result before the piezo solve.");
        temperature = &model.getVariable(position);
        if (temperature->getApplyType() != R_VARIABLE_APPLY_NODE ||
            temperature->getNVectors() != 1 || temperature->getNValues() != model.getNNodes())
            fail("Heat->Piezo requires one absolute-kelvin temperature per current mesh node.");
    }

    double geometryScale = 1.0;
    for (const auto &node : data.truss.nodes)
        geometryScale = std::max(geometryScale,vectorNorm(subtract(node,axisStart)));
    for (std::size_t node=0;node<nNodes;node++)
    {
        const auto relative = subtract(data.truss.nodes[node],axisStart);
        const double coordinate = dot(relative,coupled.axis);
        const std::array<double,3> offAxis = {{
            relative[0]-coordinate*coupled.axis[0],
            relative[1]-coordinate*coupled.axis[1],
            relative[2]-coordinate*coupled.axis[2]
        }};
        if (vectorNorm(offAxis) > 1.0e-10*geometryScale)
        {
            fail(QString("Nodal-potential piezo mode requires collinear TRUSS1 geometry; node %1 lies off the member axis.")
                     .arg(node));
        }
        input.coordinates[node] = coordinate;
    }

    std::vector<bool> incident(nNodes,false);
    input.members.reserve(data.truss.members.size());
    for (std::size_t memberPosition=0;memberPosition<data.truss.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.truss.members[memberPosition];
        const RStructuralPiezoPort &port = *data.portsByMember[memberPosition];
        const double projectedLength = std::fabs(input.coordinates[member.node2]-input.coordinates[member.node1]);
        const double geometryTolerance = 1.0e-10*std::max(geometryScale,member.referenceLength);
        if (std::fabs(projectedLength-member.referenceLength) > geometryTolerance ||
            projectedLength <= geometryTolerance)
        {
            fail(QString("Nodal-potential piezo element %1 is not a nonzero member on the common axis.")
                     .arg(member.elementId));
        }
        RStructuralPiezoCoupled1DMember coupledMember;
        coupledMember.node0 = member.node1;
        coupledMember.node1 = member.node2;
        coupledMember.area = member.area;
        coupledMember.youngModulus = member.youngModulus;
        coupledMember.coupling = port.coupling;
        coupledMember.permittivity = port.permittivity;
        if (temperature)
        {
            const double first = temperature->getValue(0,member.node1);
            const double second = temperature->getValue(0,member.node2);
            if (!std::isfinite(first) || !std::isfinite(second) || first < 0.0 || second < 0.0)
                fail("Heat->Piezo encountered invalid nodal absolute temperature.");
            const double meanTemperature = 0.5*(first+second);
            const RMaterial *material = nullptr;
            for (uint lineId=0;lineId<model.getNLines();lineId++)
            {
                const RLine &line = model.getLine(lineId);
                for (uint local=0;local<line.size();local++)
                    if (line.get(local) == member.elementId)
                        material = &line.getMaterial();
            }
            if (!material) fail("Heat->Piezo member has no line material.");
            const uint alphaPosition = material->findPosition(
                RMaterialProperty::ThermalExpansionCoefficient);
            const uint modulusPosition = material->findPosition(
                RMaterialProperty::ModulusOfElasiticity);
            if (alphaPosition >= material->size() || modulusPosition >= material->size())
                fail("Heat->Piezo member requires Young modulus and thermal expansion tables.");
            coupledMember.thermalStrain = integratedThermalStrain(
                material->get(alphaPosition),meanTemperature);
            coupledMember.youngModulus = material->get(modulusPosition).get(meanTemperature);
            if (!std::isfinite(coupledMember.youngModulus) || coupledMember.youngModulus <= 0.0)
                fail("Heat->Piezo member has invalid Young modulus at current temperature.");
        }
        input.members.push_back(coupledMember);
        incident[member.node1] = true;
        incident[member.node2] = true;
    }
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!incident[node])
        {
            fail(QString("Nodal-potential piezo mode does not support isolated model node %1.").arg(node));
        }
    }

    const double axisComponentTolerance = 1.0e-14;
    for (std::size_t node=0;node<nNodes;node++)
    {
        const std::array<double,3> load = {{
            data.truss.deadNodalLoads[3*node],
            data.truss.deadNodalLoads[3*node+1],
            data.truss.deadNodalLoads[3*node+2]
        }};
        const double axialLoad = dot(load,coupled.axis);
        const std::array<double,3> transverseLoad = {{
            load[0]-axialLoad*coupled.axis[0],
            load[1]-axialLoad*coupled.axis[1],
            load[2]-axialLoad*coupled.axis[2]
        }};
        const double loadNorm = vectorNorm(load);
        if (vectorNorm(transverseLoad) > 1.0e-12*std::max(1.0,loadNorm))
        {
            fail(QString("Nodal-potential piezo mode rejects transverse active force at node %1.").arg(node));
        }
        input.mechanicalLoads[node] = axialLoad;

        bool anyPrescribedComponent = false;
        bool hasAxialPrescription = false;
        double axialDisplacement = 0.0;
        for (std::size_t component=0;component<3;component++)
        {
            const std::size_t dof = 3*node+component;
            if (!data.truss.prescribedMask[dof])
                continue;
            anyPrescribedComponent = true;
            const double value = data.truss.prescribedDisplacements[dof];
            if (std::fabs(coupled.axis[component]) <= axisComponentTolerance)
            {
                if (std::fabs(value) > 1.0e-12*std::max(1.0,std::fabs(value)))
                {
                    fail(QString("Nodal-potential piezo mode rejects transverse prescribed displacement at node %1.")
                             .arg(node));
                }
                continue;
            }
            const double candidate = value/coupled.axis[component];
            if (!std::isfinite(candidate))
                fail(QString("Nodal-potential piezo displacement at node %1 is non-finite.").arg(node));
            if (hasAxialPrescription &&
                std::fabs(candidate-axialDisplacement) >
                    1.0e-12*std::max(1.0,std::max(std::fabs(candidate),std::fabs(axialDisplacement))))
            {
                fail(QString("Nodal-potential piezo displacement conditions at node %1 do not define one axial translation.")
                         .arg(node));
            }
            axialDisplacement = candidate;
            hasAxialPrescription = true;
        }
        if (anyPrescribedComponent && hasAxialPrescription)
        {
            for (std::size_t component=0;component<3;component++)
            {
                const std::size_t dof = 3*node+component;
                if (data.truss.prescribedMask[dof] &&
                    std::fabs(data.truss.prescribedDisplacements[dof]-
                              axialDisplacement*coupled.axis[component]) >
                        1.0e-12*std::max(1.0,std::fabs(data.truss.prescribedDisplacements[dof])))
                {
                    fail(QString("Nodal-potential piezo prescribed displacement at node %1 is not collinear with the member axis.")
                             .arg(node));
                }
            }
            input.displacementPrescribed[node] = true;
            input.prescribedDisplacements[node] = axialDisplacement;
        }
    }

    const RStructuralPiezoSetup &setup = model.getProblemSetup().getStructuralPiezoSetup();
    for (const RStructuralPiezoNodeBoundary &boundary : setup.nodalBoundaryConditions)
    {
        if (boundary.nodeId >= nNodes)
            fail("Piezo nodal electrical condition references an invalid model node.");
        if (boundary.type == RStructuralPiezoBoundaryType::Potential)
        {
            input.potentialPrescribed[boundary.nodeId] = true;
            input.prescribedPotentials[boundary.nodeId] = boundary.value;
        }
        else if (boundary.type == RStructuralPiezoBoundaryType::ConductorCharge)
        {
            input.prescribedConductorCharges[boundary.nodeId] = boundary.value;
        }
        else
        {
            fail(QString("Piezo node %1 has an unsupported electrical boundary type.").arg(boundary.nodeId));
        }
    }

    const RMatrixSolverConf &gmres = model.getMatrixSolverConf(RMatrixSolverConf::GMRES);
    input.gmresRestart = gmres.getNInnerIterations();
    input.gmresCycles = gmres.getNOuterIterations();
    input.gmresRelativeTolerance = gmres.getSolverCvgValue();
    return coupled;
}

QString modelInputFingerprint(const RModel &model,const PiezoModelData &data)
{
    QByteArray bytes("range-fea/piezo/truss-input/v1");
    appendUInt(bytes,model.getNElements());
    appendUInt(bytes,model.getProblemTaskTree().getProblemTypeMask());
    appendSize(bytes,data.truss.nodes.size());
    for (const auto &node : data.truss.nodes)
        for (double coordinate : node) appendDouble(bytes,coordinate);
    appendSize(bytes,data.truss.members.size());
    for (const auto &member : data.truss.members)
    {
        appendUInt(bytes,member.elementId);
        appendUInt(bytes,member.node1);
        appendUInt(bytes,member.node2);
        appendDouble(bytes,member.youngModulus);
        appendDouble(bytes,member.area);
        appendDouble(bytes,member.referenceLength);
    }
    appendBoolVector(bytes,data.truss.prescribedMask);
    appendDoubleVector(bytes,data.truss.prescribedDisplacements);
    appendDoubleVector(bytes,data.truss.deadNodalLoads);
    const auto &setup = model.getProblemSetup().getStructuralPiezoSetup();
    appendString(bytes,setup.meshFingerprint);
    appendUInt(bytes,static_cast<uint>(setup.electricalMode));
    appendSize(bytes,setup.ports.size());
    for (const auto &port : setup.ports)
    {
        appendString(bytes,port.id);
        appendUInt(bytes,port.elementId);
        appendDouble(bytes,port.coupling);
        appendDouble(bytes,port.permittivity);
        appendDouble(bytes,port.voltage);
    }
    appendSize(bytes,setup.passiveElementIds.size());
    for (uint elementId : setup.passiveElementIds) appendUInt(bytes,elementId);
    appendSize(bytes,setup.nodalBoundaryConditions.size());
    for (const auto &boundary : setup.nodalBoundaryConditions)
    {
        appendUInt(bytes,boundary.nodeId);
        appendUInt(bytes,static_cast<uint>(boundary.type));
        appendDouble(bytes,boundary.value);
    }
    if (setup.electricalMode == RStructuralPiezoElectricalMode::NodalPotential)
    {
        // Rebuild the same actual coupled inputs, including temperature-evaluated
        // E and integrated thermal strain, rather than hashing only a profile.
        const auto coupled = buildCoupledPiezoModelData(model,data,true);
        const auto &input = coupled.input;
        for (double component : coupled.axis) appendDouble(bytes,component);
        appendDoubleVector(bytes,input.coordinates);
        appendSize(bytes,input.members.size());
        for (const auto &member : input.members)
        {
            appendUInt(bytes,member.node0);
            appendUInt(bytes,member.node1);
            appendDouble(bytes,member.area);
            appendDouble(bytes,member.youngModulus);
            appendDouble(bytes,member.coupling);
            appendDouble(bytes,member.permittivity);
            appendDouble(bytes,member.thermalStrain);
        }
        appendDoubleVector(bytes,input.mechanicalLoads);
        appendDoubleVector(bytes,input.prescribedConductorCharges);
        appendBoolVector(bytes,input.displacementPrescribed);
        appendDoubleVector(bytes,input.prescribedDisplacements);
        appendBoolVector(bytes,input.potentialPrescribed);
        appendDoubleVector(bytes,input.prescribedPotentials);
        appendUInt(bytes,input.gmresRestart);
        appendUInt(bytes,input.gmresCycles);
        appendDouble(bytes,input.gmresRelativeTolerance);
        appendDouble(bytes,input.mechanicalAbsoluteTolerance);
        appendDouble(bytes,input.mechanicalRelativeTolerance);
        appendDouble(bytes,input.electricAbsoluteTolerance);
        appendDouble(bytes,input.electricRelativeTolerance);
        if (model.getProblemTaskTree().getProblemTypeMask() == (R_PROBLEM_HEAT | R_PROBLEM_STRESS))
        {
            // Temperature is a consumed input here, even though it was produced
            // by a preceding heat solve. Hash each raw node value as well as the
            // effective member strain so distinct fields with equal averages
            // cannot reuse one another's response certificate.
            const auto &temperature = model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE));
            appendUInt(bytes,temperature.getNValues());
            for (uint node=0;node<temperature.getNValues();node++)
                appendDouble(bytes,temperature.getValue(0,node));
        }
    }
    else
    {
        const auto &cg = model.getMatrixSolverConf(RMatrixSolverConf::CG);
        appendUInt(bytes,cg.getNInnerIterations());
        appendUInt(bytes,cg.getNOuterIterations());
        appendDouble(bytes,cg.getSolverCvgValue());
    }
    return digest(bytes);
}

struct FreeSystem
{
    RSparseMatrix tangent;
    std::vector<double> rhs;
    std::vector<std::size_t> globalToFree;
    std::vector<std::size_t> freeToGlobal;
    double freeExternalLoadNorm = 0.0;
};

FreeSystem assembleFreeSystem(const PiezoModelData &data)
{
    const std::size_t dofCount = data.truss.prescribedMask.size();
    if (dofCount > std::numeric_limits<uint>::max())
        fail("Piezo model has too many DOFs for the sparse solver index type.");

    FreeSystem system;
    system.globalToFree.assign(dofCount,noFreeDof);
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!data.truss.prescribedMask[dof])
        {
            system.globalToFree[dof] = system.freeToGlobal.size();
            system.freeToGlobal.push_back(dof);
        }
    }
    if (system.freeToGlobal.size() > std::numeric_limits<uint>::max())
        fail("Piezo model has too many free DOFs for the sparse solver index type.");

    system.tangent.setNRows(uint(system.freeToGlobal.size()));
    system.rhs.assign(system.freeToGlobal.size(),0.0);
    for (std::size_t free=0;free<system.freeToGlobal.size();free++)
    {
        const double load = data.truss.deadNodalLoads[system.freeToGlobal[free]];
        if (!std::isfinite(load))
            fail("Piezo model contains a non-finite free nodal load.");
        system.rhs[free] = load;
        system.freeExternalLoadNorm = std::hypot(system.freeExternalLoadNorm,load);
    }
    if (!std::isfinite(system.freeExternalLoadNorm))
        fail("Piezo model free-load norm is non-finite.");

    for (std::size_t memberPosition=0;memberPosition<data.truss.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.truss.members[memberPosition];
        const RStructuralPiezoPort *port = data.portsByMember[memberPosition];
        const auto &x1 = data.truss.nodes[member.node1];
        const auto &x2 = data.truss.nodes[member.node2];
        const double inverseLength = 1.0/member.referenceLength;
        const std::array<double,3> direction = {{
            (x2[0]-x1[0])*inverseLength,
            (x2[1]-x1[1])*inverseLength,
            (x2[2]-x1[2])*inverseLength
        }};
        const std::array<double,6> b = {{
            -direction[0],-direction[1],-direction[2],
             direction[0], direction[1], direction[2]
        }};
        const std::array<std::size_t,6> globalDofs = {{
            std::size_t(member.node1)*3,std::size_t(member.node1)*3+1,std::size_t(member.node1)*3+2,
            std::size_t(member.node2)*3,std::size_t(member.node2)*3+1,std::size_t(member.node2)*3+2
        }};
        const double axialStiffness = member.youngModulus*member.area*inverseLength;
        const double voltageForce = port ? port->coupling*member.area*inverseLength*port->voltage : 0.0;
        if (!std::isfinite(axialStiffness) || axialStiffness <= 0.0 ||
            !std::isfinite(voltageForce))
        {
            fail(QString("Piezo element %1 has non-finite axial or coupling stiffness.").arg(member.elementId));
        }

        for (std::size_t row=0;row<6;row++)
        {
            const std::size_t freeRow = system.globalToFree[globalDofs[row]];
            if (freeRow == noFreeDof)
                continue;

            // Internal force is K*u + b*(e*A/L)*V, so move the prescribed-
            // voltage term to the external side with a negative sign.
            system.rhs[freeRow] -= b[row]*voltageForce;
            if (!std::isfinite(system.rhs[freeRow]))
                fail("Piezo equivalent voltage load is non-finite.");

            for (std::size_t column=0;column<6;column++)
            {
                const double value = axialStiffness*b[row]*b[column];
                if (!std::isfinite(value))
                    fail("Piezo element tangent contains a non-finite entry.");
                if (value == 0.0)
                    continue;
                const std::size_t freeColumn = system.globalToFree[globalDofs[column]];
                if (freeColumn == noFreeDof)
                {
                    system.rhs[freeRow] -= value*data.truss.prescribedDisplacements[globalDofs[column]];
                }
                else
                {
                    const double updated = system.tangent.findValue(uint(freeRow),uint(freeColumn))+value;
                    if (!std::isfinite(updated))
                        fail("Piezo sparse tangent assembly is non-finite.");
                    system.tangent.addValue(uint(freeRow),uint(freeColumn),value);
                }
                if (!std::isfinite(system.rhs[freeRow]))
                    fail("Piezo prescribed-displacement load shift is non-finite.");
            }
        }
    }
    return system;
}

RStructuralPiezoModelResult recoverResult(const PiezoModelData &data,
                                          const std::vector<double> &displacements)
{
    RStructuralPiezoModelResult result;
    result.displacements = displacements;
    result.internalForces.assign(displacements.size(),0.0);
    result.memberAxialForces.resize(data.truss.members.size(),0.0);
    result.memberRightCharges.resize(data.truss.members.size(),0.0);
    result.memberElectricalPortAvailable.resize(data.truss.members.size(),false);

    for (std::size_t memberPosition=0;memberPosition<data.truss.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.truss.members[memberPosition];
        const RStructuralPiezoPort *port = data.portsByMember[memberPosition];
        const auto &x1 = data.truss.nodes[member.node1];
        const auto &x2 = data.truss.nodes[member.node2];
        const double inverseLength = 1.0/member.referenceLength;
        const std::array<double,3> direction = {{
            (x2[0]-x1[0])*inverseLength,
            (x2[1]-x1[1])*inverseLength,
            (x2[2]-x1[2])*inverseLength
        }};
        const std::array<double,6> b = {{
            -direction[0],-direction[1],-direction[2],
             direction[0], direction[1], direction[2]
        }};
        const std::array<std::size_t,6> globalDofs = {{
            std::size_t(member.node1)*3,std::size_t(member.node1)*3+1,std::size_t(member.node1)*3+2,
            std::size_t(member.node2)*3,std::size_t(member.node2)*3+1,std::size_t(member.node2)*3+2
        }};
        double elongation = 0.0;
        for (std::size_t local=0;local<6;local++)
            elongation += b[local]*displacements[globalDofs[local]];
        if (!std::isfinite(elongation))
            fail(QString("Piezo element %1 has non-finite axial elongation.").arg(member.elementId));

        double mechanicalForce = 0.0;
        double memberEnthalpy = 0.0;
        if (port)
        {
            const RStructuralPiezoKernelResult response = RStructuralPiezoKernel::evaluate(
                elongation,port->voltage,member.youngModulus,port->coupling,
                port->permittivity,member.area,member.referenceLength);
            mechanicalForce = response.mechanicalForce;
            memberEnthalpy = response.enthalpy;
            result.memberRightCharges[memberPosition] = -response.dielectricDisplacementFlux;
            result.memberElectricalPortAvailable[memberPosition] = true;
        }
        else
        {
            // Passive bars have no electrical kernel, material or branch.
            mechanicalForce = member.youngModulus*member.area*inverseLength*elongation;
            memberEnthalpy = 0.5*mechanicalForce*elongation;
        }
        if (!std::isfinite(mechanicalForce) || !std::isfinite(memberEnthalpy))
            fail("Piezo mechanical member recovery is non-finite.");
        result.memberAxialForces[memberPosition] = mechanicalForce;
        result.enthalpy += memberEnthalpy;
        if (!std::isfinite(result.enthalpy))
            fail("Total piezo electric enthalpy is non-finite.");

        for (std::size_t local=0;local<6;local++)
        {
            const std::size_t dof = globalDofs[local];
            result.internalForces[dof] += b[local]*mechanicalForce;
            if (!std::isfinite(result.internalForces[dof]))
                fail("Piezo internal nodal force is non-finite.");
        }
    }

    std::vector<double> freeResidual;
    freeResidual.reserve(displacements.size());
    for (std::size_t dof=0;dof<displacements.size();dof++)
    {
        if (!data.truss.prescribedMask[dof])
        {
            freeResidual.push_back(result.internalForces[dof]-data.truss.deadNodalLoads[dof]);
        }
    }
    result.residualNorm = stableNorm(freeResidual);
    return result;
}

void storeElectricCharge(RModel &model,const std::vector<double> &elementCharges)
{
    const uint position = model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod
        ? model.addVariable(R_VARIABLE_ELECTRIC_CHARGE)
        : model.findVariable(R_VARIABLE_ELECTRIC_CHARGE);
    RVariable &charge = model.getVariable(position);
    charge.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    charge.resize(1,uint(elementCharges.size()));
    double minimum = 0.0;
    double maximum = 0.0;
    for (std::size_t element=0;element<elementCharges.size();element++)
    {
        const double value = elementCharges[element];
        charge.setValue(0,uint(element),value);
        minimum = std::min(minimum,value);
        maximum = std::max(maximum,value);
    }
    charge.getVariableData().setMinMaxDisplayValue(minimum,maximum);
}

void storeNodalScalar(RModel &model,RVariableType type,const std::vector<double> &values)
{
    const uint position = model.findVariable(type) == RConstants::eod
        ? model.addVariable(type)
        : model.findVariable(type);
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(1,uint(values.size()));
    double minimum = 0.0;
    double maximum = 0.0;
    for (std::size_t node=0;node<values.size();node++)
    {
        const double value = values[node];
        if (!std::isfinite(value))
            fail(QString("Piezo nodal output '%1' contains a non-finite value.")
                     .arg(RVariable::getName(type)));
        variable.setValue(0,uint(node),value);
        minimum = std::min(minimum,value);
        maximum = std::max(maximum,value);
    }
    variable.getVariableData().setMinMaxDisplayValue(minimum,maximum);
}

RStructuralPiezoModelResult solveNodalPotential(const RModel &model,
                                                const PiezoModelData &data)
{
    CoupledPiezoModelData coupled = buildCoupledPiezoModelData(model,data,true);
    RStructuralPiezoCoupled1DResult solved;
    try
    {
        solved = RStructuralPiezoCoupled1D::solve(coupled.input);
    }
    catch (const std::exception &exception)
    {
        fail(QString::fromUtf8(exception.what()));
    }

    RStructuralPiezoModelResult result;
    const std::size_t nNodes = data.truss.nodes.size();
    result.displacements.assign(3*nNodes,0.0);
    result.internalForces.assign(3*nNodes,0.0);
    result.memberAxialForces = solved.memberAxialForces;
    result.memberRightCharges = solved.memberRightCharges;
    result.memberElectricalPortAvailable.assign(data.truss.members.size(),true);
    result.nodalPotentials = solved.potentials;
    result.nodalConductorCharges = solved.conductorCharges;
    result.electricalMode = RStructuralPiezoElectricalMode::NodalPotential;
    result.residualNorm = solved.mechanicalResidualNorm;
    result.electricResidualNorm = solved.electricResidualNorm;
    result.enthalpy = solved.enthalpy;

    if (solved.displacements.size() != nNodes || solved.internalForces.size() != nNodes)
        fail("Coupled piezo core returned mechanical arrays with invalid dimensions.");
    for (std::size_t node=0;node<nNodes;node++)
    {
        for (std::size_t component=0;component<3;component++)
        {
            result.displacements[3*node+component] = solved.displacements[node]*coupled.axis[component];
            result.internalForces[3*node+component] = solved.internalForces[node]*coupled.axis[component];
        }
    }
    return result;
}
}

void RStructuralPiezoModelAdapter::validate(const RModel &model)
{
    if (!model.getProblemSetup().getStructuralPiezoSetup().enabled)
        return;
    const PiezoModelData data = buildPiezoModelData(model);
    if (model.getProblemSetup().getStructuralPiezoSetup().electricalMode ==
        RStructuralPiezoElectricalMode::NodalPotential)
    {
        (void)buildCoupledPiezoModelData(model,data);
    }
}

RStructuralPiezoModelResult RStructuralPiezoModelAdapter::solve(const RModel &model)
{
    const PiezoModelData data = buildPiezoModelData(model);
    if (model.getProblemSetup().getStructuralPiezoSetup().electricalMode ==
        RStructuralPiezoElectricalMode::NodalPotential)
    {
        auto result = solveNodalPotential(model,data);
        result.inputFingerprint = modelInputFingerprint(model,data);
        result.resultSeal = resultCertificate(result,result.inputFingerprint);
        return result;
    }
    FreeSystem system = assembleFreeSystem(data);
    const RMatrixSolverConf &solverConf = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    const double configuredTolerance = solverConf.getSolverCvgValue();
    if (!std::isfinite(configuredTolerance) || configuredTolerance < 0.0)
        fail("Piezo CG convergence tolerance must be finite and nonnegative.");

    std::vector<double> displacement(data.truss.prescribedDisplacements);
    if (!system.freeToGlobal.empty())
    {
        for (std::size_t row=0;row<system.freeToGlobal.size();row++)
        {
            const double diagonal = system.tangent.findValue(uint(row),uint(row));
            if (!std::isfinite(diagonal) || diagonal <= 0.0)
            {
                fail("Piezo free tangent has a nonpositive diagonal; it is singular or contains an unconstrained DOF.");
            }
        }

        RRVector rhs(uint(system.rhs.size()),0.0);
        for (std::size_t row=0;row<system.rhs.size();row++) rhs[uint(row)] = system.rhs[row];
        RMatrixSolver linearSolver(solverConf);
        linearSolver.disableConvergenceLogFile();
        RRVector freeSolution;
        linearSolver.solve(system.tangent,rhs,freeSolution,
                           R_MATRIX_PRECONDITIONER_JACOBI,1);
        if (freeSolution.getNRows() != rhs.getNRows())
            fail("Piezo CG returned a solution with invalid dimensions.");

        std::vector<double> freeValues(freeSolution.begin(),freeSolution.end());
        const double effectiveRhsNorm = stableNorm(system.rhs);
        const double linearTolerance =
            std::max(1.0e-8,configuredTolerance)*
            std::max(effectiveRhsNorm,system.freeExternalLoadNorm)+1.0e-12;
        if (!std::isfinite(linearTolerance))
            fail("Piezo linear residual tolerance is non-finite.");

        RRVector product;
        RSparseMatrix::mlt(system.tangent,freeSolution,product);
        if (product.getNRows() != rhs.getNRows())
            fail("Piezo fresh linear residual has invalid dimensions.");
        std::vector<double> freshLinearResidual(system.rhs.size(),0.0);
        for (std::size_t row=0;row<system.rhs.size();row++)
        {
            if (!std::isfinite(freeValues[row]) || !std::isfinite(product[uint(row)]))
                fail("Piezo CG returned a non-finite solution or matrix product.");
            freshLinearResidual[row] = system.rhs[row]-product[uint(row)];
        }
        if (stableNorm(freshLinearResidual) > linearTolerance)
            fail("Piezo CG failed fresh b-Au residual certification.");

        for (std::size_t free=0;free<system.freeToGlobal.size();free++)
            displacement[system.freeToGlobal[free]] = freeValues[free];
    }

    RStructuralPiezoModelResult result = recoverResult(data,displacement);
    const double effectiveRhsNorm = stableNorm(system.rhs);
    const double physicalTolerance = std::max(1.0e-8,configuredTolerance)*
        std::max(effectiveRhsNorm,system.freeExternalLoadNorm)+1.0e-12;
    if (!std::isfinite(physicalTolerance) || result.residualNorm > physicalTolerance)
        fail("Piezo recovered physical free-force residual failed certification.");
    result.inputFingerprint = modelInputFingerprint(model,data);
    result.resultSeal = resultCertificate(result,result.inputFingerprint);
    return result;
}

void RStructuralPiezoModelAdapter::clearResults(RModel &model)
{
    RStructuralStaticTrussModel::clearResults(model);
    const uint potentialPosition = model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL);
    if (potentialPosition != RConstants::eod)
        model.removeVariable(potentialPosition);
}

void RStructuralPiezoModelAdapter::publish(
    RModel &model,const RStructuralPiezoModelResult &result)
{
    const PiezoModelData data = buildPiezoModelData(model);
    const RStructuralPiezoElectricalMode configuredMode =
        model.getProblemSetup().getStructuralPiezoSetup().electricalMode;
    if (result.electricalMode != configuredMode)
        fail("Piezo result electrical mode does not match the current model setup.");
    const bool nodalPotentialMode =
        configuredMode == RStructuralPiezoElectricalMode::NodalPotential;
    const std::size_t dofCount = std::size_t(model.getNNodes())*3;
    if (result.displacements.size() != dofCount || result.internalForces.size() != dofCount ||
        result.memberAxialForces.size() != data.truss.members.size() ||
        result.memberRightCharges.size() != data.truss.members.size() ||
        result.memberElectricalPortAvailable.size() != data.truss.members.size() ||
        (nodalPotentialMode &&
         (result.nodalPotentials.size() != model.getNNodes() ||
          result.nodalConductorCharges.size() != model.getNNodes() ||
          !std::isfinite(result.electricResidualNorm) || result.electricResidualNorm < 0.0)) ||
        !std::isfinite(result.residualNorm) || result.residualNorm < 0.0 ||
        !std::isfinite(result.enthalpy))
    {
        fail("Piezo converged response dimensions or scalar values do not match the current model.");
    }
    for (double value : result.displacements)
        if (!std::isfinite(value)) fail("Piezo displacement contains a non-finite value.");
    for (double value : result.internalForces)
        if (!std::isfinite(value)) fail("Piezo internal force contains a non-finite value.");
    for (double value : result.memberAxialForces)
        if (!std::isfinite(value)) fail("Piezo axial force contains a non-finite value.");
    for (double value : result.memberRightCharges)
        if (!std::isfinite(value)) fail("Piezo right-electrode charge contains a non-finite value.");
    if (nodalPotentialMode)
    {
        for (double value : result.nodalPotentials)
            if (!std::isfinite(value)) fail("Piezo nodal potential contains a non-finite value.");
        for (double value : result.nodalConductorCharges)
            if (!std::isfinite(value)) fail("Piezo nodal conductor charge contains a non-finite value.");
    }

    for (std::size_t member=0;member<data.truss.members.size();member++)
    {
        const bool active = data.portsByMember[member] != nullptr;
        if (result.memberElectricalPortAvailable[member] != active)
            fail("Piezo result active/passive electrical port identity does not match the current setup.");
        if (!active && result.memberRightCharges[member] != 0.0)
            fail("Passive member charge must remain a zero placeholder, not an electrical output.");
    }

    if (result.inputFingerprint.isEmpty() ||
        result.inputFingerprint != modelInputFingerprint(model,data))
        fail("Piezo result input provenance does not match the current model physical inputs.");
    if (result.resultSeal.isEmpty() ||
        result.resultSeal != resultCertificate(result,result.inputFingerprint))
        fail("Piezo result seal is invalid: public solved result fields were modified or not adapter-owned.");

    std::vector<double> elementCharges(model.getNElements(),0.0);
    for (std::size_t member=0;member<data.truss.members.size();member++)
    {
        const RStructuralStaticTrussMember &trussMember = data.truss.members[member];
        if (trussMember.elementId >= elementCharges.size())
            fail("Piezo charge result references an invalid element ID.");
        const double stress = result.memberAxialForces[member]/trussMember.area;
        if (!std::isfinite(stress))
            fail(QString("Piezo stress recovery is non-finite for element %1.").arg(trussMember.elementId));
        if (!nodalPotentialMode && data.portsByMember[member])
            elementCharges[trussMember.elementId] = result.memberRightCharges[member];
    }
    for (double value : elementCharges)
        if (!std::isfinite(value)) fail("Piezo element charge output contains a non-finite value.");

    // All dimensions, recovered quantities and element slots are certified
    // before the shared mechanics helper or charge variable mutates the model.
    RStructuralStaticTrussModel::publishMechanical(
        model,result.displacements,result.internalForces,result.memberAxialForces);
    if (nodalPotentialMode)
    {
        storeNodalScalar(model,R_VARIABLE_ELECTRIC_POTENTIAL,result.nodalPotentials);
        storeNodalScalar(model,R_VARIABLE_ELECTRIC_CHARGE,result.nodalConductorCharges);
    }
    else
    {
        storeElectricCharge(model,elementCharges);
    }
}
