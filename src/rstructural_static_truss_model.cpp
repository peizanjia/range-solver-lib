#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include <rbl_error.h>

#include <rml_element_group.h>
#include <rml_line.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_structural_control_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_structural_interaction.h>
#include <rml_structural_piezo_setup.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

#include "rstructural_static_truss_model.h"

namespace
{

constexpr double referenceTemperature = 293.15;

[[noreturn]] void fail(const QString &message)
{
    const QByteArray utf8 = message.toUtf8();
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",utf8.constData());
}

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double staticComponentValue(const RConditionComponent &component,
                            const QString &description)
{
    if (component.size() != 1)
    {
        fail(QString("Static TRUSS1 mode requires one constant value for %1.").arg(description));
    }
    const double value = component.getValue(0);
    if (!std::isfinite(value))
    {
        fail(QString("Static TRUSS1 mode requires a finite value for %1.").arg(description));
    }
    return value;
}

void validateEnvironmentConditions(const RElementGroup &group, uint groupId,
                                   bool allowThermalPiezo)
{
    for (uint conditionPosition=0;
         conditionPosition<group.getNEnvironmentConditions();
         conditionPosition++)
    {
        const REnvironmentCondition &condition = group.getEnvironmentCondition(conditionPosition);
        switch (condition.getType())
        {
            case R_ENVIRONMENT_CONDITION_G_ACCELERATION:
            {
                for (uint componentPosition=0;componentPosition<condition.size();componentPosition++)
                {
                    const RConditionComponent &component = condition.getComponent(componentPosition);
                    const RVariableType type = component.getType();
                    if (type != R_VARIABLE_G_ACCELERATION_X &&
                        type != R_VARIABLE_G_ACCELERATION_Y &&
                        type != R_VARIABLE_G_ACCELERATION_Z)
                    {
                        fail(QString("Static TRUSS1 mode does not support gravity component '%1' on group %2.")
                                 .arg(RVariable::getName(type)).arg(groupId));
                    }
                    if (!component.getEnabled())
                    {
                        continue;
                    }
                    const double value = staticComponentValue(
                        component,
                        QString("gravity component '%1' on group %2")
                            .arg(RVariable::getName(type)).arg(groupId));
                    if (value != 0.0)
                    {
                        fail(QString("Static TRUSS1 mode requires enabled gravity components to be zero on group %1.")
                                 .arg(groupId));
                    }
                }
                break;
            }
            case R_ENVIRONMENT_CONDITION_TEMPERATURE:
            {
                for (uint componentPosition=0;componentPosition<condition.size();componentPosition++)
                {
                    const RConditionComponent &component = condition.getComponent(componentPosition);
                    if (component.getType() != R_VARIABLE_TEMPERATURE)
                    {
                        fail(QString("Static TRUSS1 mode does not support temperature component '%1' on group %2.")
                                 .arg(RVariable::getName(component.getType())).arg(groupId));
                    }
                    if (!component.getEnabled())
                    {
                        continue;
                    }
                    const double value = staticComponentValue(
                        component,QString("temperature on group %1").arg(groupId));
                    const double tolerance = 1.0e-12*referenceTemperature;
                    if (value < 0.0 ||
                        (!allowThermalPiezo && std::fabs(value-referenceTemperature) > tolerance))
                    {
                        fail(QString("Static TRUSS1 mode is isothermal at 293.15 K; group %1 has a different temperature.")
                                 .arg(groupId));
                    }
                }
                break;
            }
            default:
                fail(QString("Static TRUSS1 mode does not support environment condition '%1' on group %2.")
                         .arg(REnvironmentCondition::getName(condition.getType())).arg(groupId));
        }
    }
}

std::vector<uint> uniqueGroupNodes(const RModel &model,
                                  const RElementGroup &group)
{
    std::set<uint> uniqueNodes;
    for (uint i=0;i<group.size();i++)
    {
        const uint elementId = group.get(i);
        if (elementId >= model.getNElements())
        {
            fail(QString("Static TRUSS1 boundary group references invalid element %1.").arg(elementId));
        }
        const RElement &element = model.getElement(elementId);
        for (uint localNode=0;localNode<element.size();localNode++)
        {
            const uint nodeId = element.getNodeId(localNode);
            if (nodeId >= model.getNNodes())
            {
                fail(QString("Element %1 in static TRUSS1 boundary group references invalid node %2.")
                         .arg(elementId).arg(nodeId));
            }
            uniqueNodes.insert(nodeId);
        }
    }
    return std::vector<uint>(uniqueNodes.begin(),uniqueNodes.end());
}

void prescribe(RStructuralStaticTrussModelData &data,
               uint node,
               uint component,
               double value,
               bool exactDeclarations)
{
    const std::size_t dof = std::size_t(3)*node+component;
    if (dof >= data.prescribedMask.size())
    {
        fail("Static TRUSS1 displacement condition addresses an invalid node component.");
    }
    if (data.prescribedMask[dof])
    {
        const double previous = data.prescribedDisplacements[dof];
        const double tolerance = 1.0e-12*std::max(1.0,std::max(std::fabs(previous),std::fabs(value)));
        if (exactDeclarations ? previous != value : std::fabs(previous-value) > tolerance)
        {
            fail(QString("Contradictory static TRUSS1 displacement prescriptions at node %1, component %2.")
                     .arg(node).arg(component));
        }
        return;
    }
    data.prescribedMask[dof] = true;
    data.prescribedDisplacements[dof] = value;
}

void addDeadLoad(RStructuralStaticTrussModelData &data,
                 uint node,
                 uint component,
                 double value)
{
    const std::size_t dof = std::size_t(3)*node+component;
    if (dof >= data.deadNodalLoads.size())
    {
        fail("Static TRUSS1 point force addresses an invalid node component.");
    }
    data.deadNodalLoads[dof] += value;
    if (!std::isfinite(data.deadNodalLoads[dof]))
    {
        fail("Combined static TRUSS1 point forces are non-finite.");
    }
}

void applyBoundaryCondition(REntityGroupType groupType,
                            const RBoundaryCondition &condition,
                            const std::vector<uint> &nodes,
                            RStructuralStaticTrussModelData &data,
                            uint groupId,
                            bool exactDeclarations)
{
    const RBoundaryConditionType type = condition.getType();
    const bool displacement = type == R_BOUNDARY_CONDITION_DISPLACEMENT;
    const bool force = type == R_BOUNDARY_CONDITION_FORCE;
    if (!displacement && !force)
    {
        fail(QString("Static TRUSS1 mode does not support boundary condition '%1' on group %2.")
                 .arg(RBoundaryCondition::getName(type)).arg(groupId));
    }
    if (groupType != R_ENTITY_GROUP_POINT && groupType != R_ENTITY_GROUP_LINE)
    {
        fail(QString("Static TRUSS1 mode supports displacement conditions only on point/line groups; group %1 is unsupported.")
                 .arg(groupId));
    }
    if (force && groupType != R_ENTITY_GROUP_POINT)
    {
        fail(QString("Static TRUSS1 mode supports dead nodal FORCE only on point groups; line group %1 is unsupported.")
                 .arg(groupId));
    }
    if (condition.getHasLocalDirection())
    {
        fail(QString("Static TRUSS1 mode does not support local-frame boundary conditions on group %1.").arg(groupId));
    }
    if (nodes.empty())
    {
        fail(QString("Static TRUSS1 boundary condition on empty group %1 would have no effect.").arg(groupId));
    }

    for (uint i=0;i<condition.size();i++)
    {
        const RConditionComponent &component = condition.getComponent(i);
        if (!component.getEnabled())
        {
            continue;
        }
        uint axis = RConstants::eod;
        if (displacement)
        {
            if (component.getType() == R_VARIABLE_DISPLACEMENT_X) axis = 0;
            else if (component.getType() == R_VARIABLE_DISPLACEMENT_Y) axis = 1;
            else if (component.getType() == R_VARIABLE_DISPLACEMENT_Z) axis = 2;
            else
                fail(QString("Static TRUSS1 mode does not support displacement component '%1' on group %2.")
                         .arg(RVariable::getName(component.getType())).arg(groupId));
        }
        else
        {
            if (component.getType() == R_VARIABLE_FORCE_X) axis = 0;
            else if (component.getType() == R_VARIABLE_FORCE_Y) axis = 1;
            else if (component.getType() == R_VARIABLE_FORCE_Z) axis = 2;
            else
                fail(QString("Static TRUSS1 mode does not support force component '%1' on group %2.")
                         .arg(RVariable::getName(component.getType())).arg(groupId));
        }

        const double value = staticComponentValue(
            component,
            QString("%1 component on group %2")
                .arg(RVariable::getName(component.getType())).arg(groupId));
        if (displacement)
        {
            for (uint node : nodes)
            {
                prescribe(data,node,axis,value,exactDeclarations);
            }
        }
        else
        {
            const double nodalValue = value/double(nodes.size());
            if (!std::isfinite(nodalValue))
            {
                fail(QString("Static TRUSS1 point force on group %1 cannot be divided over its nodes.").arg(groupId));
            }
            for (uint node : nodes)
            {
                addDeadLoad(data,node,axis,nodalValue);
            }
        }
    }
}

void setVariableRange(RVariable &variable, const std::vector<double> &values)
{
    if (values.empty())
    {
        variable.getVariableData().setMinMaxDisplayValue(0.0,0.0);
        return;
    }
    const auto bounds = std::minmax_element(values.begin(),values.end());
    variable.getVariableData().setMinMaxDisplayValue(*bounds.first,*bounds.second);
}

void storeScalar(RModel &model,
                 RVariableType type,
                 const std::vector<double> &values)
{
    uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        position = model.addVariable(type);
    }
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    variable.resize(1,uint(values.size()));
    for (uint i=0;i<values.size();i++)
    {
        variable.setValue(0,i,values[i]);
    }
    setVariableRange(variable,values);
}

void storeVector(RModel &model,
                 RVariableType type,
                 uint nNodes,
                 const std::vector<double> &values,
                 bool magnitudeRange)
{
    uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        position = model.addVariable(type);
    }
    RVariable &variable = model.getVariable(position);
    variable.setApplyType(R_VARIABLE_APPLY_NODE);
    variable.resize(3,nNodes);
    std::vector<double> rangeValues;
    rangeValues.reserve(nNodes);
    for (uint node=0;node<nNodes;node++)
    {
        const double x = values[3*node];
        const double y = values[3*node+1];
        const double z = values[3*node+2];
        variable.setValue(0,node,x);
        variable.setValue(1,node,y);
        variable.setValue(2,node,z);
        rangeValues.push_back(magnitudeRange ? std::hypot(std::hypot(x,y),z)
                                             : std::max(std::fabs(x),std::max(std::fabs(y),std::fabs(z))));
    }
    setVariableRange(variable,rangeValues);
}

}

RStructuralStaticTrussModelData RStructuralStaticTrussModel::build(
    const RModel &model,bool allowStructuralDynamics,bool allowThermalPiezo)
{
    return buildWithPolicy(model,allowStructuralDynamics,allowThermalPiezo,ExtractionPolicy::Legacy);
}

RStructuralStaticTrussModelData RStructuralStaticTrussModel::buildMixedDynamics(const RModel &model)
{
    return buildWithPolicy(model,true,false,ExtractionPolicy::NativeMixedDynamics);
}

RStructuralStaticTrussModelData RStructuralStaticTrussModel::buildWithPolicy(
    const RModel &model,bool allowStructuralDynamics,bool allowThermalPiezo,ExtractionPolicy policy)
{
    const bool mixed = policy == ExtractionPolicy::NativeMixedDynamics;
    const RProblemTypeMask taskMask = model.getProblemTaskTree().getProblemTypeMask();
    if (taskMask != R_PROBLEM_STRESS &&
        !(allowThermalPiezo && taskMask == (R_PROBLEM_HEAT | R_PROBLEM_STRESS)))
    {
        fail("Static TRUSS1 mode requires an R_PROBLEM_STRESS-only task.");
    }
    if (model.getTimeSolver().getEnabled() && !allowStructuralDynamics)
    {
        fail("Static TRUSS1 mode does not support the legacy time solver.");
    }
    const RProblemSetup &problemSetup = model.getProblemSetup();
    const auto &mixedSetup = problemSetup.getStructuralTrussDynamicsSetup();
    if (!mixed && (mixedSetup.enabled || problemSetup.getStructuralTrussDynamicsRestartData().available))
        fail("Legacy TRUSS1 extraction cannot consume enabled native mixed setup or native restart data.");
    if (mixed)
    {
        if (!mixedSetup.enabled) fail("Native mixed TRUSS1 extraction requires enabled native setup.");
        mixedSetup.validate(model);
        if (!problemSetup.getStructuralDynamicsSetup().getModalDampingRatios().empty())
            fail("Native mixed TRUSS1 extraction does not accept modal damping records.");
        if (problemSetup.getStructuralCableSetup().enabled || problemSetup.getStructuralPiezoSetup().enabled)
            fail("Native mixed TRUSS1 extraction rejects enabled legacy cable/piezo modes.");
        if (problemSetup.getRestart() && !problemSetup.getStructuralTrussDynamicsRestartData().available)
            fail("Native mixed TRUSS1 restart requires its dedicated accepted-state record.");
        const auto &thermal = problemSetup.getStructuralThermalControlSetup();
        if (thermal.enabled || !thermal.channels.empty() || !thermal.sparseObservationNodeIds.empty() ||
            thermal.thermalSubstepSeconds != 0.0 || thermal.thermalSubstepsPerMechanicalMacrostep != 1 ||
            thermal.mechanicalMacrostepsPerControllerSample != 1)
            fail("Native mixed TRUSS1 extraction rejects thermal scheduling.");
    }
    if (problemSetup.getRestart() && !mixed)
    {
        fail("Static TRUSS1 mode does not support restart.");
    }
    if (problemSetup.getStructuralDynamicsSetup().getEnabled() && !allowStructuralDynamics)
    {
        fail("Static TRUSS1 mode does not support second-order dynamics.");
    }
    if (problemSetup.getStructuralDynamicsRestartData().available)
    {
        fail("Static TRUSS1 mode does not support structural restart data.");
    }
    const RStructuralControlSetup &control = problemSetup.getStructuralControlSetup();
    if (control.enabled || !control.controllerId.isEmpty() || !control.channelConditioning.empty() ||
        (mixed && (control.deterministicSeed != 0 || control.controlProtocolVersion != 1)))
    {
        fail("Static TRUSS1 mode does not support structural control setup.");
    }
    const RStructuralInteractionSetup &interaction = problemSetup.getStructuralInteractionSetup();
    if (!interaction.getLinearConstraints().empty() || !interaction.getCouplings().empty() ||
        !interaction.getSensors().empty() || !interaction.getActuators().empty() ||
        (!mixed && !interaction.getPrescribedMotions().empty()))
    {
        fail("Static TRUSS1 mode does not support MPCs, couplings, sensors, actuators or prescribed motion.");
    }
    if (mixed)
        for (const auto &motion : interaction.getPrescribedMotions())
        {
            motion.validate();
            if (uint(motion.target.dof) > 2 || motion.target.nodeId >= model.getNNodes())
                fail(QString("Native mixed motion '%1' requires a valid translational node target.").arg(motion.id));
        }
    if (model.findVariable(R_VARIABLE_TEMPERATURE) != RConstants::eod && !allowThermalPiezo)
    {
        fail("Static TRUSS1 mode is isothermal at 293.15 K and does not accept a temperature field.");
    }
    if (model.getNNodes() == 0)
    {
        fail("Static TRUSS1 mode requires at least one node.");
    }

    RStructuralStaticTrussModelData data;
    data.nodes.resize(model.getNNodes());
    for (uint node=0;node<model.getNNodes();node++)
    {
        const RNode &rNode = model.getNode(node);
        data.nodes[node] = {{rNode.getX(),rNode.getY(),rNode.getZ()}};
        for (double coordinate : data.nodes[node])
        {
            if (!std::isfinite(coordinate))
            {
                fail(QString("Static TRUSS1 node %1 has a non-finite reference coordinate.").arg(node));
            }
        }
    }

    const uint nElements = model.getNElements();
    std::vector<uint> lineOwnerCount(nElements,0);
    std::vector<double> areaByElement(nElements,0.0);
    std::vector<double> modulusByElement(nElements,0.0);
    std::set<uint> effectiveActiveIds;
    if (mixed)
        for (const auto &member : mixedSetup.activeMembers) effectiveActiveIds.insert(member.physicalElementId);
    for (uint lineId=0;lineId<model.getNLines();lineId++)
    {
        const RLine &line = model.getLine(lineId);
        if (line.empty())
        {
            continue;
        }
        const double area = line.getCrossArea();
        if (!finitePositive(area))
        {
            fail(QString("Static TRUSS1 line group %1 requires a finite positive cross-sectional area.").arg(lineId));
        }
        bool requiresModulus = !mixed;
        for (uint member=0;member<line.size() && !requiresModulus;++member)
        {
            const uint elementId = line.get(member);
            if (elementId >= nElements)
                fail(QString("Native mixed line group %1 references invalid element %2.").arg(lineId).arg(elementId));
            requiresModulus = effectiveActiveIds.count(model.getElementPhysicalId(elementId)) == 0;
        }
        double youngModulus = 0.0;
        if (requiresModulus)
        {
            const RMaterial &material = line.getMaterial();
            const uint modulusPosition = material.findPosition(RMaterialProperty::ModulusOfElasiticity);
            if (modulusPosition >= material.size())
                fail(QString("Static TRUSS1 line group %1 has no Young's modulus.").arg(lineId));
            const RMaterialProperty &modulusProperty = material.get(modulusPosition);
            if (modulusProperty.size() == 0)
                fail(QString("Static TRUSS1 line group %1 has an empty Young's modulus table.").arg(lineId));
            youngModulus = modulusProperty.get(referenceTemperature);
            if (!finitePositive(youngModulus))
                fail(QString("Static TRUSS1 line group %1 requires a finite positive Young's modulus at 293.15 K.").arg(lineId));
            if (!finitePositive(area*youngModulus))
                fail(QString("Static TRUSS1 line group %1 has invalid axial rigidity EA.").arg(lineId));
        }
        for (uint member=0;member<line.size();member++)
        {
            const uint elementId = line.get(member);
            if (elementId >= nElements)
            {
                fail(QString("Static TRUSS1 line group %1 references invalid element %2.").arg(lineId).arg(elementId));
            }
            const RElement &element = model.getElement(elementId);
            if (element.getType() != R_ELEMENT_TRUSS1)
            {
                fail(QString("Every static truss line member must be TRUSS1; element %1 is '%2'.")
                         .arg(elementId).arg(RElement::getName(element.getType())));
            }
            lineOwnerCount[elementId]++;
            areaByElement[elementId] = area;
            modulusByElement[elementId] = mixed && effectiveActiveIds.count(model.getElementPhysicalId(elementId))
                ? 0.0 : youngModulus;
        }
    }

    uint nonPointElementCount = 0;
    for (uint elementId=0;elementId<nElements;elementId++)
    {
        const RElement &element = model.getElement(elementId);
        if (element.getType() == R_ELEMENT_POINT)
        {
            if (lineOwnerCount[elementId] != 0)
            {
                fail(QString("Point element %1 cannot be owned by a static truss line group.").arg(elementId));
            }
            continue;
        }
        nonPointElementCount++;
        if (element.getType() != R_ELEMENT_TRUSS1)
        {
            fail(QString("Static TRUSS1 mode only supports TRUSS1 nonpoint elements; element %1 is '%2'.")
                     .arg(elementId).arg(RElement::getName(element.getType())));
        }
        if (element.size() != 2)
        {
            fail(QString("Static TRUSS1 element %1 must have two nodes.").arg(elementId));
        }
        const uint node1 = element.getNodeId(0);
        const uint node2 = element.getNodeId(1);
        if (node1 >= model.getNNodes() || node2 >= model.getNNodes() || node1 == node2)
        {
            fail(QString("Static TRUSS1 element %1 has invalid node connectivity.").arg(elementId));
        }
        if (lineOwnerCount[elementId] != 1)
        {
            fail(QString("Static TRUSS1 element %1 must be uniquely owned by one line group.").arg(elementId));
        }
        const bool effectiveActive = mixed && effectiveActiveIds.count(model.getElementPhysicalId(elementId));
        if (!finitePositive(areaByElement[elementId]) || (!effectiveActive && !finitePositive(modulusByElement[elementId])))
        {
            fail(QString("Static TRUSS1 element %1 has no usable area/material assignment.").arg(elementId));
        }

        const auto &p1 = data.nodes[node1];
        const auto &p2 = data.nodes[node2];
        const double dx = p2[0]-p1[0];
        const double dy = p2[1]-p1[1];
        const double dz = p2[2]-p1[2];
        const double length = std::hypot(std::hypot(dx,dy),dz);
        if (!finitePositive(length))
        {
            fail(QString("Static TRUSS1 element %1 has invalid reference length.").arg(elementId));
        }

        RStructuralStaticTrussMember member;
        member.elementId = elementId;
        member.node1 = node1;
        member.node2 = node2;
        member.area = areaByElement[elementId];
        member.youngModulus = modulusByElement[elementId];
        member.referenceLength = length;
        data.members.push_back(member);
    }
    if (nonPointElementCount == 0)
    {
        fail("Static TRUSS1 mode requires at least one nonpoint member.");
    }

    const std::size_t nDofs = std::size_t(3)*model.getNNodes();
    data.prescribedMask.assign(nDofs,false);
    data.prescribedDisplacements.assign(nDofs,0.0);
    data.deadNodalLoads.assign(nDofs,0.0);
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        const RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group)
        {
            continue;
        }
        validateEnvironmentConditions(*group,groupId,allowThermalPiezo);
        if (group->getNInitialConditions() != 0 && !mixed)
        {
            fail(QString("Static TRUSS1 mode does not support initial conditions on group %1.").arg(groupId));
        }
        if (mixed)
            for (uint ci=0;ci<group->getNInitialConditions();++ci)
            {
                const auto type = group->getInitialCondition(ci).getType();
                if (type != R_INITIAL_CONDITION_DISPLACEMENT && type != R_INITIAL_CONDITION_VELOCITY)
                    fail(QString("Native mixed TRUSS1 extraction only accepts displacement/velocity initial conditions on physical group %1.").arg(group->getPhysicalId()));
            }
        if (group->getNBoundaryConditions() == 0)
        {
            continue;
        }
        const REntityGroupType groupType = model.getEntityGroupType(groupId);
        const std::vector<uint> nodes = uniqueGroupNodes(model,*group);
        for (uint bcPosition=0;bcPosition<group->getNBoundaryConditions();bcPosition++)
        {
            const RBoundaryCondition &condition = group->getBoundaryCondition(bcPosition);
            if (allowThermalPiezo &&
                (condition.getType() == R_BOUNDARY_CONDITION_TEMPERATURE ||
                 condition.getType() == R_BOUNDARY_CONDITION_HEAT ||
                 condition.getType() == R_BOUNDARY_CONDITION_HEAT_RATE_AREA ||
                 condition.getType() == R_BOUNDARY_CONDITION_HEAT_RATE_VOLUME))
                continue;
            applyBoundaryCondition(groupType,condition,
                                   nodes,data,groupId,mixed);
        }
    }

    return data;
}

void RStructuralStaticTrussModel::clearResults(RModel &model)
{
    const RVariableType outputTypes[] =
    {
        R_VARIABLE_DISPLACEMENT,
        R_VARIABLE_ROTATION,
        R_VARIABLE_FORCE,
        R_VARIABLE_MOMENT,
        R_VARIABLE_STRESS_X,
        R_VARIABLE_STRESS_Y,
        R_VARIABLE_STRESS_Z,
        R_VARIABLE_STRESS_YZ,
        R_VARIABLE_STRESS_XZ,
        R_VARIABLE_STRESS_XY,
        R_VARIABLE_STRESS_VON_MISES,
        R_VARIABLE_STRESS_NORMAL,
        R_VARIABLE_STRESS_SHEAR,
        R_VARIABLE_ELECTRIC_CHARGE
    };
    for (RVariableType type : outputTypes)
    {
        const uint position = model.findVariable(type);
        if (position != RConstants::eod)
        {
            model.removeVariable(position);
        }
    }
}

void RStructuralStaticTrussModel::publishMechanical(
    RModel &model,
    const std::vector<double> &displacements,
    const std::vector<double> &internalForces,
    const std::vector<double> &memberAxialForces,
    bool allowStructuralDynamics)
{
    publishWithPolicy(model,displacements,internalForces,memberAxialForces,
                      allowStructuralDynamics,ExtractionPolicy::Legacy);
}

void RStructuralStaticTrussModel::publishMechanicalMixed(
    RModel &model,const std::vector<double> &displacements,
    const std::vector<double> &internalForces,const std::vector<double> &memberAxialForces)
{
    publishWithPolicy(model,displacements,internalForces,memberAxialForces,
                      true,ExtractionPolicy::NativeMixedDynamics);
}

void RStructuralStaticTrussModel::publishWithPolicy(
    RModel &model,const std::vector<double> &displacements,
    const std::vector<double> &internalForces,const std::vector<double> &memberAxialForces,
    bool allowStructuralDynamics,ExtractionPolicy policy)
{
    const bool mixed = policy == ExtractionPolicy::NativeMixedDynamics;
    const bool thermalPiezo = !mixed && model.getProblemSetup().getStructuralPiezoSetup().enabled &&
        model.getProblemTaskTree().getProblemTypeMask() ==
            (R_PROBLEM_HEAT | R_PROBLEM_STRESS);
    const RStructuralStaticTrussModelData data = buildWithPolicy(
        model,allowStructuralDynamics,thermalPiezo,policy);
    // Native q/v/a and prescribed kinematics are certified by the adapter's
    // driver before this output-only helper runs. A motion owns absolute q;
    // data.prescribedDisplacements is only the source held-DOF declaration,
    // not a baseline to add to, or compare with, a moving accepted q.
    const std::size_t expectedDofs = std::size_t(3)*model.getNNodes();
    if (displacements.size() != expectedDofs || internalForces.size() != expectedDofs ||
        memberAxialForces.size() != data.members.size())
    {
        fail("Converged static TRUSS1 response dimensions do not match the current model.");
    }
    for (double value : displacements)
    {
        if (!std::isfinite(value)) fail("Converged static TRUSS1 displacement contains a non-finite value.");
    }
    for (double value : internalForces)
    {
        if (!std::isfinite(value)) fail("Converged static TRUSS1 internal force contains a non-finite value.");
    }
    for (double value : memberAxialForces)
    {
        if (!std::isfinite(value)) fail("Converged static TRUSS1 member axial force contains a non-finite value.");
    }
    for (uint node=0;node<model.getNNodes();node++)
    {
        const double magnitude = std::hypot(
            std::hypot(displacements[3*node],displacements[3*node+1]),
            displacements[3*node+2]);
        if (!std::isfinite(magnitude))
        {
            fail("Converged static TRUSS1 displacement magnitude is non-finite.");
        }
    }

    std::vector<double> normalStress(model.getNElements(),0.0);
    std::vector<double> vonMisesStress(model.getNElements(),0.0);
    std::vector<double> localStressX(model.getNElements(),0.0);
    for (std::size_t memberPosition=0;memberPosition<data.members.size();memberPosition++)
    {
        const RStructuralStaticTrussMember &member = data.members[memberPosition];
        const double stress = memberAxialForces[memberPosition]/member.area;
        if (!std::isfinite(stress))
        {
            fail(QString("Static TRUSS1 stress recovery is non-finite for element %1.").arg(member.elementId));
        }
        normalStress[member.elementId] = stress;
        vonMisesStress[member.elementId] = std::fabs(stress);
        localStressX[member.elementId] = stress;
    }

    // All response dimensions, finite values, magnitudes and recovered stresses
    // have passed validation before the first candidate-model mutation.
    if (mixed)
    {
        const uint potential = model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL);
        if (potential != RConstants::eod) model.removeVariable(potential);
    }
    storeVector(model,R_VARIABLE_DISPLACEMENT,model.getNNodes(),displacements,true);
    storeVector(model,R_VARIABLE_FORCE,model.getNNodes(),internalForces,false);
    storeScalar(model,R_VARIABLE_STRESS_NORMAL,normalStress);
    storeScalar(model,R_VARIABLE_STRESS_VON_MISES,vonMisesStress);
    storeScalar(model,R_VARIABLE_STRESS_X,localStressX);
}
