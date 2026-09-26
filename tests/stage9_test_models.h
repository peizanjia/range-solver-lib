#ifndef STAGE9_TEST_MODELS_H
#define STAGE9_TEST_MODELS_H

#include <array>

#include <rml_boundary_condition.h>
#include <rml_element_group.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_interaction.h>
#include <rml_variable.h>

namespace Stage9TestModels
{
inline constexpr double YoungModulus = 2.0e5;
inline constexpr double Density = 1.0e3;
inline constexpr double Area = 1.0e-2;
inline constexpr double Length = 2.0;
inline constexpr double CableTension = 100.0;

inline void setComponent(RCondition &condition,
                         RVariableType type,
                         bool enabled,
                         double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled)
    {
        component.add(0.0,value);
    }
}

inline void addPoint(RModel &model, uint node, uint pointGroup)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,pointGroup);
}

inline RBoundaryCondition displacement(const std::array<bool,3> &enabled,
                                       const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {
        R_VARIABLE_DISPLACEMENT_X,
        R_VARIABLE_DISPLACEMENT_Y,
        R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;++axis)
    {
        setComponent(condition,types[axis],enabled[axis],values[axis]);
    }
    return condition;
}

inline void assignMaterial(RLine &line, bool includeDensity = true)
{
    line.setCrossArea(Area);
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,YoungModulus);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.0);
    material.add(poisson);
    if (includeDensity)
    {
        RMaterialProperty density(RMaterialProperty::Density);
        density.add(293.15,Density);
        material.add(density);
    }
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);
    line.setMaterial(material);
}

inline void setZeroGravity(RModel &model)
{
    const RVariableType types[3] = {
        R_VARIABLE_G_ACCELERATION_X,
        R_VARIABLE_G_ACCELERATION_Y,
        R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();++groupId)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group)
        {
            continue;
        }
        REnvironmentCondition &gravity =
            group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for (RVariableType type : types)
        {
            setComponent(gravity,type,true,0.0);
        }
    }
}

inline void addCatalog(RModel &model, uint tipNode,
                       const RR3Vector &direction = RR3Vector(1.0,0.0,0.0))
{
    RStructuralInteractionSetup &interaction =
        model.getProblemSetup().getStructuralInteractionSetup();
    RStructuralActuator actuator;
    actuator.id = QStringLiteral("tip-force-x");
    actuator.target.nodeId = tipNode;
    actuator.variable = R_VARIABLE_FORCE_X;
    actuator.direction = direction;
    actuator.unit = RVariable::getUnits(actuator.variable);
    actuator.minimum = -1.0e6;
    actuator.maximum = 1.0e6;
    interaction.getActuators().push_back(actuator);

    const RStructuralSensorQuantity quantities[3] = {
        RStructuralSensorQuantity::Displacement,
        RStructuralSensorQuantity::Velocity,
        RStructuralSensorQuantity::Acceleration
    };
    const QString ids[3] = {
        QStringLiteral("tip-displacement-x"),
        QStringLiteral("tip-velocity-x"),
        QStringLiteral("tip-acceleration-x")
    };
    for (uint i=0;i<3;++i)
    {
        RStructuralSensor sensor;
        sensor.id = ids[i];
        sensor.target.nodeId = tipNode;
        sensor.quantity = quantities[i];
        sensor.direction = direction;
        sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
        interaction.getSensors().push_back(sensor);
    }
}

//! A fixed-root axial truss. The optional line constraints remove the two
//! transverse mechanisms while leaving the tip's axial translation free.
inline RModel makeTruss(uint nElements = 1,
                        double tipLoad = 0.0,
                        double rootDisplacement = 0.0,
                        bool constrainTransverse = true,
                        bool includeDensity = true)
{
    if (nElements == 0)
    {
        nElements = 1;
    }
    RModel model;
    model.setNNodes(nElements+1);
    for (uint node=0;node<=nElements;++node)
    {
        model.getNode(node).set(Length*double(node)/double(nElements),0.0,0.0);
        model.setNodePhysicalId(node,7101+17*node);
    }
    for (uint elementId=0;elementId<nElements;++elementId)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,elementId);
        element.setNodeId(1,elementId+1);
        model.addElement(element,true,0);
        model.setElementPhysicalId(elementId,8101+19*elementId);
    }
    assignMaterial(model.getLine(0),includeDensity);
    if (constrainTransverse)
    {
        RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
        setComponent(transverse,R_VARIABLE_DISPLACEMENT_X,false);
        setComponent(transverse,R_VARIABLE_DISPLACEMENT_Y,true,0.0);
        setComponent(transverse,R_VARIABLE_DISPLACEMENT_Z,true,0.0);
        model.getLine(0).addBoundaryCondition(transverse);
    }

    addPoint(model,0,0);
    addPoint(model,nElements,1);
    model.getPoint(0).addBoundaryCondition(
        displacement({{true,true,true}},{{rootDisplacement,0.0,0.0}}));
    if (tipLoad != 0.0)
    {
        RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
        setComponent(force,R_VARIABLE_FORCE_X,true,tipLoad);
        setComponent(force,R_VARIABLE_FORCE_Y,false);
        setComponent(force,R_VARIABLE_FORCE_Z,false);
        model.getPoint(1).addBoundaryCondition(force);
    }

    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RMatrixSolverConf &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setOutputFrequency(0);
    linear.setSolverCvgValue(1.0e-13);
    linear.setNInnerIterations(200);
    linear.setNOuterIterations(2000);
    setZeroGravity(model);
    addCatalog(model,nElements);
    return model;
}

//! Rotate the axial truss to global Y and make global X/Z its held directions.
//! Its one free equation is then a useful check of exported local-frame axes.
inline void rotateTrussToY(RModel &model)
{
    for (uint node=0;node<model.getNNodes();++node)
    {
        const double y = model.getNode(node).getX();
        model.getNode(node).set(0.0,y,0.0);
    }
    RBoundaryCondition &transverse = model.getLine(0).getBoundaryCondition(0);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_X,true,0.0);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Y,false);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Z,true,0.0);
    RStructuralInteractionSetup &interaction =
        model.getProblemSetup().getStructuralInteractionSetup();
    interaction.getActuators()[0].variable = R_VARIABLE_FORCE_Y;
    interaction.getActuators()[0].id = QStringLiteral("tip-force-y");
    interaction.getActuators()[0].direction = RR3Vector(0.0,1.0,0.0);
    for (RStructuralSensor &sensor : interaction.getSensors())
    {
        sensor.direction = RR3Vector(0.0,1.0,0.0);
    }
    interaction.getSensors()[0].id = QStringLiteral("tip-displacement-y");
    interaction.getSensors()[1].id = QStringLiteral("tip-velocity-y");
    interaction.getSensors()[2].id = QStringLiteral("tip-acceleration-y");
}

inline RModel makeFrameCantilever(double tipForceY = 0.0,
                                  double tipMomentZ = 0.0)
{
    constexpr double iy = 2.0e-5;
    constexpr double iz = 5.0e-5;
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(Length,0.0,0.0);
    model.setNodePhysicalId(0,9203);
    model.setNodePhysicalId(1,9277);

    RElement beam(R_ELEMENT_BEAM1);
    beam.setNodeId(0,0);
    beam.setNodeId(1,1);
    model.addElement(beam,true,0);
    RLine &line = model.getLine(0);
    assignMaterial(line);
    line.setSecondMomentY(iy);
    line.setSecondMomentZ(iz);
    line.setTorsionalConstant(iy+iz);
    line.setSectionModulusY(4.0e-4);
    line.setSectionModulusZ(5.0e-4);
    line.setTorsionalSectionModulus(6.0e-4);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));

    addPoint(model,0,0);
    addPoint(model,1,1);
    model.getPoint(0).addBoundaryCondition(
        displacement({{true,true,true}},{{0.0,0.0,0.0}}));
    RBoundaryCondition rotation(R_BOUNDARY_CONDITION_ROTATION);
    setComponent(rotation,R_VARIABLE_ROTATION_X,true,0.0);
    setComponent(rotation,R_VARIABLE_ROTATION_Y,true,0.0);
    setComponent(rotation,R_VARIABLE_ROTATION_Z,true,0.0);
    model.getPoint(0).addBoundaryCondition(rotation);
    if (tipForceY != 0.0)
    {
        RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
        setComponent(force,R_VARIABLE_FORCE_X,false);
        setComponent(force,R_VARIABLE_FORCE_Y,true,tipForceY);
        setComponent(force,R_VARIABLE_FORCE_Z,false);
        model.getPoint(1).addBoundaryCondition(force);
    }
    if (tipMomentZ != 0.0)
    {
        RBoundaryCondition moment(R_BOUNDARY_CONDITION_MOMENT);
        setComponent(moment,R_VARIABLE_MOMENT_X,false);
        setComponent(moment,R_VARIABLE_MOMENT_Y,false);
        setComponent(moment,R_VARIABLE_MOMENT_Z,true,tipMomentZ);
        model.getPoint(1).addBoundaryCondition(moment);
    }

    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RMatrixSolverConf &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setOutputFrequency(0);
    linear.setSolverCvgValue(1.0e-13);
    linear.setNInnerIterations(200);
    linear.setNOuterIterations(2000);
    setZeroGravity(model);
    addCatalog(model,1);
    return model;
}

inline RModel makeTautCable(double tipLoad = CableTension)
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(Length,0.0,0.0);
    model.setNodePhysicalId(0,6301);
    model.setNodePhysicalId(1,6397);

    RElement cable(R_ELEMENT_TRUSS1);
    cable.setNodeId(0,0);
    cable.setNodeId(1,1);
    model.addElement(cable,true,0);
    assignMaterial(model.getLine(0));
    addPoint(model,0,0);
    addPoint(model,1,1);
    model.getPoint(0).addBoundaryCondition(
        displacement({{true,true,true}},{{0.0,0.0,0.0}}));
    if (tipLoad != 0.0)
    {
        RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
        setComponent(force,R_VARIABLE_FORCE_X,true,tipLoad);
        setComponent(force,R_VARIABLE_FORCE_Y,false);
        setComponent(force,R_VARIABLE_FORCE_Z,false);
        model.getPoint(1).addBoundaryCondition(force);
    }

    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RMatrixSolverConf &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setOutputFrequency(0);
    linear.setSolverCvgValue(1.0e-13);
    linear.setNInnerIterations(200);
    linear.setNOuterIterations(2000);
    setZeroGravity(model);
    RStructuralCableSetup &cableSetup = model.getProblemSetup().getStructuralCableSetup();
    cableSetup.enabled = true;
    cableSetup.members.push_back({0,Length,CableTension});
    cableSetup.bindToModel(model);
    addCatalog(model,1);
    return model;
}
}

#endif // STAGE9_TEST_MODELS_H
