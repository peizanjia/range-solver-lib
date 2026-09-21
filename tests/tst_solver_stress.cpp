#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "rml_model.h"
#include "rsolver.h"

namespace
{

// Steel with a zero Poisson ratio, so that an axially loaded bar behaves
// exactly one dimensionally and can be compared against a closed form.
constexpr double elasticityModulus = 2.0e11; // [Pa]
constexpr double poissonRatio = 0.0;
constexpr double density = 8000.0;           // [kg/m^3]
constexpr double thermalExpansion = 0.0;     // [1/K]

constexpr double barLength = 2.0;            // [m]
constexpr double barSide = 0.1;              // [m]
constexpr double barArea = barSide * barSide;
constexpr uint   nSlices = 30;
constexpr double gravityAcceleration = -9.80665; // [m/s^2]

uint nodeId(uint i, uint j, uint k)
{
    return 4*i + 2*j + k;
}

void addTetrahedron(RModel &model, uint n1, uint n2, uint n3, uint n4)
{
    // Keep a positive signed volume, otherwise the Jacobian determinant and
    // with it the element stiffness comes out negative.
    RR3Vector a(model.getNode(n1).toVector());
    RR3Vector b(model.getNode(n2).toVector());
    RR3Vector c(model.getNode(n3).toVector());
    RR3Vector d(model.getNode(n4).toVector());

    RR3Vector ab,ac,ad,cross;
    RR3Vector::subtract(b,a,ab);
    RR3Vector::subtract(c,a,ac);
    RR3Vector::subtract(d,a,ad);
    RR3Vector::cross(ab,ac,cross);

    if (RR3Vector::dot(cross,ad) < 0.0)
    {
        std::swap(n3,n4);
    }

    RElement element(R_ELEMENT_TETRA1);
    element.setNodeId(0,n1);
    element.setNodeId(1,n2);
    element.setNodeId(2,n3);
    element.setNodeId(3,n4);
    model.addElement(element,true,0);
}

//! Build a fixed-free bar of tetrahedra along the x axis.
//!
//! The lateral degrees of freedom of the whole volume are removed with a
//! Displacement boundary condition whose X component is switched off, which
//! makes the bar a clean one dimensional problem and at the same time exercises
//! the per-component switch of an optional boundary condition.
RModel buildBar(double tipForce,
                bool lateralRestraint = true,
                RElementType volumeElementType = R_ELEMENT_TETRA1)
{
    RModel model;

    model.setNNodes(4*(nSlices+1));
    for (uint i=0;i<=nSlices;i++)
    {
        double x = barLength * double(i) / double(nSlices);
        for (uint j=0;j<2;j++)
        {
            for (uint k=0;k<2;k++)
            {
                model.getNode(nodeId(i,j,k)).set(x,double(j)*barSide,double(k)*barSide);
            }
        }
    }

    // Use the same structured nodes for every first-order solid so the element
    // formulations can be compared against the same closed-form bar solution.
    for (uint i=0;i<nSlices;i++)
    {
        uint v000 = nodeId(i,0,0);
        uint v001 = nodeId(i,0,1);
        uint v010 = nodeId(i,1,0);
        uint v011 = nodeId(i,1,1);
        uint v100 = nodeId(i+1,0,0);
        uint v101 = nodeId(i+1,0,1);
        uint v110 = nodeId(i+1,1,0);
        uint v111 = nodeId(i+1,1,1);

        if (volumeElementType == R_ELEMENT_TETRA1)
        {
            // Six tetrahedra per slice - a Kuhn decomposition of the brick.
            addTetrahedron(model,v000,v100,v110,v111);
            addTetrahedron(model,v000,v110,v010,v111);
            addTetrahedron(model,v000,v010,v011,v111);
            addTetrahedron(model,v000,v011,v001,v111);
            addTetrahedron(model,v000,v001,v101,v111);
            addTetrahedron(model,v000,v101,v100,v111);
        }
        else if (volumeElementType == R_ELEMENT_WEDGE1)
        {
            RElement lower(R_ELEMENT_WEDGE1);
            lower.setNodeId(0,v000);
            lower.setNodeId(1,v010);
            lower.setNodeId(2,v001);
            lower.setNodeId(3,v100);
            lower.setNodeId(4,v110);
            lower.setNodeId(5,v101);
            model.addElement(lower,true,0);

            RElement upper(R_ELEMENT_WEDGE1);
            upper.setNodeId(0,v011);
            upper.setNodeId(1,v001);
            upper.setNodeId(2,v010);
            upper.setNodeId(3,v111);
            upper.setNodeId(4,v101);
            upper.setNodeId(5,v110);
            model.addElement(upper,true,0);
        }
        else if (volumeElementType == R_ELEMENT_HEXA1)
        {
            RElement element(R_ELEMENT_HEXA1);
            element.setNodeId(0,v000);
            element.setNodeId(1,v100);
            element.setNodeId(2,v110);
            element.setNodeId(3,v010);
            element.setNodeId(4,v001);
            element.setNodeId(5,v101);
            element.setNodeId(6,v111);
            element.setNodeId(7,v011);
            model.addElement(element,true,0);
        }
        else
        {
            throw std::invalid_argument("Unsupported volume element in the bar test model");
        }
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Steel");

    RMaterialProperty modulusProperty(RMaterialProperty::ModulusOfElasiticity);
    modulusProperty.add(293.15,elasticityModulus);
    material.add(modulusProperty);

    RMaterialProperty poissonProperty(RMaterialProperty::PoissonRatio);
    poissonProperty.add(293.15,poissonRatio);
    material.add(poissonProperty);

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);

    RMaterialProperty expansionProperty(RMaterialProperty::ThermalExpansionCoefficient);
    expansionProperty.add(293.15,thermalExpansion);
    material.add(expansionProperty);

    RVolume &volume = model.getVolume(0);
    volume.setName("Bar");
    volume.setMaterial(material);

    // Lateral restraint - X switched off, Y and Z held at zero. It is left out
    // when a node of the model carries a local frame, because a globally phrased
    // constraint is read in that local frame and would over-constrain the node.
    if (lateralRestraint)
    {
        RBoundaryCondition lateralBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
        lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
        lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
        lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
        volume.addBoundaryCondition(lateralBc);
    }

    // Fixed end - the four nodes at x = 0.
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            RElement element(R_ELEMENT_POINT);
            element.setNodeId(0,nodeId(0,j,k));
            model.addElement(element,true,0);
        }
    }

    RPoint &fixedPoint = model.getPoint(0);
    fixedPoint.setName("Fixed end");
    fixedPoint.setVolume(0.0);

    RBoundaryCondition fixedBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    fixedPoint.addBoundaryCondition(fixedBc);

    // Loaded end - the four nodes at x = L. The diagonal nodes of the two-wedge
    // section have twice the tributary area of the other corners, so use four
    // point entities there to apply the consistent nodal force distribution.
    const double wedgeWeights[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            const uint corner = 2*j+k;
            const uint pointGroup = volumeElementType == R_ELEMENT_WEDGE1 ? 1+corner : 1;
            RElement element(R_ELEMENT_POINT);
            element.setNodeId(0,nodeId(nSlices,j,k));
            model.addElement(element,true,pointGroup);
        }
    }

    const uint nLoadedGroups = volumeElementType == R_ELEMENT_WEDGE1 ? 4 : 1;
    for (uint i=0;i<nLoadedGroups;i++)
    {
        RPoint &loadedPoint = model.getPoint(1+i);
        loadedPoint.setName("Free end");
        loadedPoint.setVolume(0.0);

        if (tipForce != 0.0)
        {
            const double groupForce = volumeElementType == R_ELEMENT_WEDGE1
                                    ? tipForce*wedgeWeights[i]
                                    : tipForce;
            RBoundaryCondition forceBc(R_BOUNDARY_CONDITION_FORCE);
            forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,groupForce);
            forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
            forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
            loadedPoint.addBoundaryCondition(forceBc);
        }
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    model.getTimeSolver().setEnabled(false);

    return model;
}

//! Build the same fixed-free bar out of truss elements.
//!
//! A truss has no lateral stiffness and no lateral mass, so the y and z degrees
//! of freedom are removed with a Displacement boundary condition on the line
//! entity that has its X component switched off.
RModel buildTruss(double tipForce)
{
    RModel model;

    model.setNNodes(nSlices+1);
    for (uint i=0;i<=nSlices;i++)
    {
        model.getNode(i).set(barLength*double(i)/double(nSlices),0.0,0.0);
    }

    for (uint i=0;i<nSlices;i++)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,i);
        element.setNodeId(1,i+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Steel");

    RMaterialProperty modulusProperty(RMaterialProperty::ModulusOfElasiticity);
    modulusProperty.add(293.15,elasticityModulus);
    material.add(modulusProperty);

    RMaterialProperty poissonProperty(RMaterialProperty::PoissonRatio);
    poissonProperty.add(293.15,poissonRatio);
    material.add(poissonProperty);

    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);

    RMaterialProperty expansionProperty(RMaterialProperty::ThermalExpansionCoefficient);
    expansionProperty.add(293.15,thermalExpansion);
    material.add(expansionProperty);

    RLine &line = model.getLine(0);
    line.setName("Bar");
    line.setCrossArea(barArea);
    line.setMaterial(material);

    RBoundaryCondition lateralBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    line.addBoundaryCondition(lateralBc);

    RElement fixedElement(R_ELEMENT_POINT);
    fixedElement.setNodeId(0,0);
    model.addElement(fixedElement,true,0);

    RPoint &fixedPoint = model.getPoint(0);
    fixedPoint.setName("Fixed end");
    fixedPoint.setVolume(0.0);

    RBoundaryCondition fixedBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    fixedPoint.addBoundaryCondition(fixedBc);

    RElement loadedElement(R_ELEMENT_POINT);
    loadedElement.setNodeId(0,nSlices);
    model.addElement(loadedElement,true,1);

    RPoint &loadedPoint = model.getPoint(1);
    loadedPoint.setName("Free end");
    loadedPoint.setVolume(0.0);

    if (tipForce != 0.0)
    {
        RBoundaryCondition forceBc(R_BOUNDARY_CONDITION_FORCE);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,tipForce);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
        loadedPoint.addBoundaryCondition(forceBc);
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);

    model.getTimeSolver().setEnabled(false);

    return model;
}

RModel buildQuadraticTruss(double tipForce)
{
    constexpr uint elements = 4;
    RModel model;
    model.setNNodes(2*elements+1);
    for (uint i=0;i<=2*elements;i++)
    {
        model.getNode(i).set(barLength*double(i)/double(2*elements),0.0,0.0);
    }
    for (uint i=0;i<elements;i++)
    {
        RElement element(R_ELEMENT_TRUSS2);
        element.setNodeId(0,2*i);
        element.setNodeId(1,2*i+2);
        element.setNodeId(2,2*i+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,elasticityModulus); material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,poissonRatio); material.add(poisson);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density); material.add(densityProperty);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,thermalExpansion); material.add(expansion);
    RLine &line = model.getLine(0);
    line.setCrossArea(barArea);
    line.setMaterial(material);

    RBoundaryCondition lateral(R_BOUNDARY_CONDITION_DISPLACEMENT);
    lateral.getComponent(lateral.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    lateral.getComponent(lateral.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    lateral.getComponent(lateral.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    line.addBoundaryCondition(lateral);

    RElement fixed(R_ELEMENT_POINT); fixed.setNodeId(0,0); model.addElement(fixed,true,0);
    RBoundaryCondition fixedBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixedBc);

    RElement loaded(R_ELEMENT_POINT); loaded.setNodeId(0,2*elements); model.addElement(loaded,true,1);
    RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,tipForce);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
    model.getPoint(1).addBoundaryCondition(force);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

RModel buildQuadraticSurfacePatch(RElementType type, double axialStrain)
{
    RModel model;
    if (type == R_ELEMENT_TRI2)
    {
        model.setNNodes(6);
        const RNode nodes[6] = {RNode(0.0,0.0,0.0),RNode(1.0,0.0,0.0),RNode(0.0,1.0,0.0),
                                RNode(0.5,0.0,0.0),RNode(0.5,0.5,0.0),RNode(0.0,0.5,0.0)};
        for (uint i=0;i<6;i++) model.getNode(i) = nodes[i];
    }
    else
    {
        model.setNNodes(9);
        const RNode nodes[9] = {RNode(0.0,0.0,0.0),RNode(1.0,0.0,0.0),RNode(1.0,1.0,0.0),
                                RNode(0.0,1.0,0.0),RNode(0.5,0.0,0.0),RNode(1.0,0.5,0.0),
                                RNode(0.5,1.0,0.0),RNode(0.0,0.5,0.0),RNode(0.5,0.5,0.0)};
        for (uint i=0;i<9;i++) model.getNode(i) = nodes[i];
    }
    RElement element(type);
    for (uint i=0;i<element.size();i++) element.setNodeId(i,i);
    model.addElement(element,true,0);

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,elasticityModulus); material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.25); material.add(poisson);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density); material.add(densityProperty);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0); material.add(expansion);
    RSurface &surface = model.getSurface(0);
    surface.setThickness(0.1);
    surface.setStructuralFormulation(RElementFormulation::PlaneStress);
    surface.setMaterial(material);

    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    surface.addBoundaryCondition(transverse);

    const std::vector<uint> left = type == R_ELEMENT_TRI2
                                 ? std::vector<uint>{0,2,5}
                                 : std::vector<uint>{0,3,7};
    const std::vector<uint> right = type == R_ELEMENT_TRI2
                                  ? std::vector<uint>{1}
                                  : std::vector<uint>{1,2,5};
    const std::vector<uint> middle = type == R_ELEMENT_TRI2
                                   ? std::vector<uint>{3,4}
                                   : std::vector<uint>{4,6,8};
    for (uint node : left)
    {
        RElement point(R_ELEMENT_POINT); point.setNodeId(0,node); model.addElement(point,true,0);
    }
    for (uint node : right)
    {
        RElement point(R_ELEMENT_POINT); point.setNodeId(0,node); model.addElement(point,true,1);
    }
    for (uint node : middle)
    {
        RElement point(R_ELEMENT_POINT); point.setNodeId(0,node); model.addElement(point,true,2);
    }
    RBoundaryCondition leftBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    leftBc.getComponent(leftBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    leftBc.getComponent(leftBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).setEnabled(false);
    leftBc.getComponent(leftBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).setEnabled(false);
    model.getPoint(0).addBoundaryCondition(leftBc);
    RBoundaryCondition rightBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    rightBc.getComponent(rightBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,axialStrain);
    rightBc.getComponent(rightBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).setEnabled(false);
    rightBc.getComponent(rightBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).setEnabled(false);
    model.getPoint(1).addBoundaryCondition(rightBc);
    RBoundaryCondition middleBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    middleBc.getComponent(middleBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.5*axialStrain);
    middleBc.getComponent(middleBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).setEnabled(false);
    middleBc.getComponent(middleBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).setEnabled(false);
    model.getPoint(2).addBoundaryCondition(middleBc);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

//! Build a clamped-free 3D Euler-Bernoulli frame along global X.
RModel buildBeam(double tipForceY, double tipForceZ)
{
    constexpr uint beamSlices = 6;
    constexpr double iy = 2.0e-5;
    constexpr double iz = 5.0e-5;

    RModel model;
    model.setNNodes(beamSlices+1);
    for (uint i=0;i<=beamSlices;i++)
    {
        model.getNode(i).set(barLength*double(i)/double(beamSlices),0.0,0.0);
    }
    for (uint i=0;i<beamSlices;i++)
    {
        RElement element(R_ELEMENT_BEAM1);
        element.setNodeId(0,i);
        element.setNodeId(1,i+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulusProperty(RMaterialProperty::ModulusOfElasiticity);
    modulusProperty.add(293.15,elasticityModulus);
    material.add(modulusProperty);
    RMaterialProperty poissonProperty(RMaterialProperty::PoissonRatio);
    poissonProperty.add(293.15,0.3);
    material.add(poissonProperty);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);
    RMaterialProperty expansionProperty(RMaterialProperty::ThermalExpansionCoefficient);
    expansionProperty.add(293.15,thermalExpansion);
    material.add(expansionProperty);

    RLine &line = model.getLine(0);
    line.setName("Frame");
    line.setCrossArea(barArea);
    line.setSecondMomentY(iy);
    line.setSecondMomentZ(iz);
    line.setTorsionalConstant(iy+iz);
    line.setSectionModulusY(4.0e-4);
    line.setSectionModulusZ(5.0e-4);
    line.setTorsionalSectionModulus(6.0e-4);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));
    line.setMaterial(material);
    REnvironmentCondition gravity(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_X)).add(0.0,0.0);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y)).add(0.0,0.0);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Z)).add(0.0,0.0);
    line.addEnvironmentCondition(gravity);
    RElement fixedElement(R_ELEMENT_POINT);
    fixedElement.setNodeId(0,0);
    model.addElement(fixedElement,true,0);
    RBoundaryCondition fixed(R_BOUNDARY_CONDITION_DISPLACEMENT);
    fixed.getComponent(fixed.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    fixed.getComponent(fixed.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    fixed.getComponent(fixed.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixed);

    RElement loadedElement(R_ELEMENT_POINT);
    loadedElement.setNodeId(0,beamSlices);
    model.addElement(loadedElement,true,1);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    load.getComponent(load.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,0.0);
    load.getComponent(load.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,tipForceY);
    load.getComponent(load.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,tipForceZ);
    model.getPoint(1).addBoundaryCondition(load);

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(1000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(10000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-13);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

RModel buildQuadraticBeam(double tipForceY)
{
    RModel model;
    model.setNNodes(3);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(0.5*barLength,0.0,0.0);
    model.getNode(2).set(barLength,0.0,0.0);
    RElement beam(R_ELEMENT_BEAM2);
    beam.setNodeId(0,0);
    beam.setNodeId(1,2);
    beam.setNodeId(2,1);
    model.addElement(beam,true,0);

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,elasticityModulus); material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.3); material.add(poisson);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density); material.add(densityProperty);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,thermalExpansion); material.add(expansion);
    RLine &line = model.getLine(0);
    line.setCrossArea(barArea);
    line.setSecondMomentY(2.0e-5);
    line.setSecondMomentZ(5.0e-5);
    line.setTorsionalConstant(7.0e-5);
    line.setSectionModulusY(4.0e-4);
    line.setSectionModulusZ(5.0e-4);
    line.setTorsionalSectionModulus(6.0e-4);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));
    line.setMaterial(material);
    REnvironmentCondition quadraticGravity(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    quadraticGravity.getComponent(quadraticGravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_X)).add(0.0,0.0);
    quadraticGravity.getComponent(quadraticGravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y)).add(0.0,0.0);
    quadraticGravity.getComponent(quadraticGravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Z)).add(0.0,0.0);
    line.addEnvironmentCondition(quadraticGravity);

    RElement fixed(R_ELEMENT_POINT); fixed.setNodeId(0,0); model.addElement(fixed,true,0);
    RBoundaryCondition displacement(R_BOUNDARY_CONDITION_DISPLACEMENT);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(displacement);
    RBoundaryCondition rotation(R_BOUNDARY_CONDITION_ROTATION);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(rotation);

    RElement loaded(R_ELEMENT_POINT); loaded.setNodeId(0,2); model.addElement(loaded,true,1);
    RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,0.0);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,tipForceY);
    force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
    model.getPoint(1).addBoundaryCondition(force);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

//! Build a membrane strip in the XY plane. Its transverse in-plane strain and
//! out-of-plane displacement are restrained so the two constitutive contracts
//! have distinct closed-form axial stiffness and transverse stress values.
RModel buildSurfaceStrip(double tipForce,
                         RElementType elementType,
                         RElementFormulation formulation,
                         double expansion = thermalExpansion,
                         double environmentTemperature = 293.15)
{
    constexpr uint surfaceSlices = 12;
    constexpr double width = 0.1;
    constexpr double thickness = 0.02;
    constexpr double surfacePoissonRatio = 0.25;

    RModel model;
    model.setNNodes(2*(surfaceSlices+1));
    for (uint i=0;i<=surfaceSlices;i++)
    {
        const double x = barLength*double(i)/double(surfaceSlices);
        model.getNode(2*i+0).set(x,0.0,0.0);
        model.getNode(2*i+1).set(x,width,0.0);
    }

    for (uint i=0;i<surfaceSlices;i++)
    {
        const uint n00 = 2*i;
        const uint n01 = n00+1;
        const uint n10 = 2*(i+1);
        const uint n11 = n10+1;
        if (elementType == R_ELEMENT_TRI1)
        {
            RElement lower(R_ELEMENT_TRI1);
            lower.setNodeId(0,n00);
            lower.setNodeId(1,n10);
            lower.setNodeId(2,n11);
            model.addElement(lower,true,0);

            RElement upper(R_ELEMENT_TRI1);
            // Keep the local x axis parallel to the strip axis so the
            // element-local stress components have the same closed form.
            upper.setNodeId(0,n11);
            upper.setNodeId(1,n01);
            upper.setNodeId(2,n00);
            model.addElement(upper,true,0);
        }
        else if (elementType == R_ELEMENT_QUAD1)
        {
            RElement element(R_ELEMENT_QUAD1);
            element.setNodeId(0,n00);
            element.setNodeId(1,n10);
            element.setNodeId(2,n11);
            element.setNodeId(3,n01);
            model.addElement(element,true,0);
        }
        else
        {
            throw std::invalid_argument("Unsupported surface element in strip model");
        }
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Membrane material");
    RMaterialProperty modulusProperty(RMaterialProperty::ModulusOfElasiticity);
    modulusProperty.add(293.15,elasticityModulus);
    material.add(modulusProperty);
    RMaterialProperty poissonProperty(RMaterialProperty::PoissonRatio);
    poissonProperty.add(293.15,surfacePoissonRatio);
    material.add(poissonProperty);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);
    RMaterialProperty expansionProperty(RMaterialProperty::ThermalExpansionCoefficient);
    expansionProperty.add(293.15,expansion);
    material.add(expansionProperty);

    RSurface &surface = model.getSurface(0);
    surface.setName("Membrane strip");
    surface.setThickness(thickness);
    surface.setStructuralFormulation(formulation);
    surface.setMaterial(material);

    REnvironmentCondition temperature(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    temperature.getComponent(temperature.findComponentPosition(R_VARIABLE_TEMPERATURE))
               .add(0.0,environmentTemperature);
    surface.addEnvironmentCondition(temperature);

    RBoundaryCondition lateralBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    lateralBc.getComponent(lateralBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    surface.addBoundaryCondition(lateralBc);

    for (uint j=0;j<2;j++)
    {
        RElement fixed(R_ELEMENT_POINT);
        fixed.setNodeId(0,j);
        model.addElement(fixed,true,0);

        RElement loaded(R_ELEMENT_POINT);
        loaded.setNodeId(0,2*surfaceSlices+j);
        model.addElement(loaded,true,1);
    }

    RBoundaryCondition fixedBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    fixedBc.getComponent(fixedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixedBc);

    if (tipForce != 0.0)
    {
        RBoundaryCondition forceBc(R_BOUNDARY_CONDITION_FORCE);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_X)).add(0.0,tipForce);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
        forceBc.getComponent(forceBc.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
        model.getPoint(1).addBoundaryCondition(forceBc);
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

//! Build an annular meridian in the global XY plane. The axisymmetric solver
//! contract uses X as radius, Y as the revolution axis and Z as circumferential.
RModel buildAxisymmetricAnnulus(uint radialSlices,
                                RElementType elementType,
                                bool restrainAxialStrain,
                                double expansion = thermalExpansion,
                                double environmentTemperature = 293.15,
                                double innerPressure = 0.0,
                                double outerPressure = 0.0)
{
    constexpr double innerRadius = 1.0;
    constexpr double outerRadius = 2.0;
    constexpr double height = 1.0;
    constexpr double axisymmetricPoissonRatio = 0.25;
    constexpr double fullAngle = 2.0*3.14159265358979323846;

    RModel model;
    model.setNNodes(2*(radialSlices+1));
    for (uint i=0;i<=radialSlices;i++)
    {
        const double radius = innerRadius
                            + (outerRadius-innerRadius)*double(i)/double(radialSlices);
        model.getNode(2*i+0).set(radius,0.0,0.0);
        model.getNode(2*i+1).set(radius,height,0.0);
    }

    for (uint i=0;i<radialSlices;i++)
    {
        const uint n00 = 2*i;
        const uint n01 = n00+1;
        const uint n10 = 2*(i+1);
        const uint n11 = n10+1;
        if (elementType == R_ELEMENT_QUAD1)
        {
            RElement element(R_ELEMENT_QUAD1);
            element.setNodeId(0,n00);
            element.setNodeId(1,n10);
            element.setNodeId(2,n11);
            element.setNodeId(3,n01);
            model.addElement(element,true,0);
        }
        else if (elementType == R_ELEMENT_TRI1)
        {
            RElement lower(R_ELEMENT_TRI1);
            lower.setNodeId(0,n00);
            lower.setNodeId(1,n10);
            lower.setNodeId(2,n11);
            model.addElement(lower,true,0);

            RElement upper(R_ELEMENT_TRI1);
            upper.setNodeId(0,n11);
            upper.setNodeId(1,n01);
            upper.setNodeId(2,n00);
            model.addElement(upper,true,0);
        }
        else
        {
            throw std::invalid_argument("Unsupported axisymmetric surface element");
        }
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Axisymmetric material");
    RMaterialProperty modulusProperty(RMaterialProperty::ModulusOfElasiticity);
    modulusProperty.add(293.15,elasticityModulus);
    material.add(modulusProperty);
    RMaterialProperty poissonProperty(RMaterialProperty::PoissonRatio);
    poissonProperty.add(293.15,axisymmetricPoissonRatio);
    material.add(poissonProperty);
    RMaterialProperty densityProperty(RMaterialProperty::Density);
    densityProperty.add(293.15,density);
    material.add(densityProperty);
    RMaterialProperty expansionProperty(RMaterialProperty::ThermalExpansionCoefficient);
    expansionProperty.add(293.15,expansion);
    material.add(expansionProperty);

    RSurface &surface = model.getSurface(0);
    surface.setName("Axisymmetric annulus");
    surface.setThickness(0.0);
    surface.setStructuralFormulation(RElementFormulation::Axisymmetric);
    surface.setMaterial(material);

    REnvironmentCondition temperature(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    temperature.getComponent(temperature.findComponentPosition(R_VARIABLE_TEMPERATURE))
               .add(0.0,environmentTemperature);
    surface.addEnvironmentCondition(temperature);

    RBoundaryCondition meridianBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    meridianBc.getComponent(meridianBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    if (restrainAxialStrain)
    {
        meridianBc.getComponent(meridianBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    }
    else
    {
        meridianBc.getComponent(meridianBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).setEnabled(false);
    }
    meridianBc.getComponent(meridianBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    surface.addBoundaryCondition(meridianBc);

    if (!restrainAxialStrain)
    {
        for (uint i=0;i<=radialSlices;i++)
        {
            RElement point(R_ELEMENT_POINT);
            point.setNodeId(0,2*i);
            model.addElement(point,true,0);
        }
        RBoundaryCondition bottomBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
        bottomBc.getComponent(bottomBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
        bottomBc.getComponent(bottomBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
        bottomBc.getComponent(bottomBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).setEnabled(false);
        model.getPoint(0).addBoundaryCondition(bottomBc);
    }

    if (innerPressure != 0.0 || outerPressure != 0.0)
    {
        for (uint axial=0;axial<2;axial++)
        {
            RElement innerPoint(R_ELEMENT_POINT);
            innerPoint.setNodeId(0,axial);
            model.addElement(innerPoint,true,0);

            RElement outerPoint(R_ELEMENT_POINT);
            outerPoint.setNodeId(0,2*radialSlices+axial);
            model.addElement(outerPoint,true,1);
        }
        if (innerPressure != 0.0)
        {
            RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_X))
                 .add(0.0,innerPressure*fullAngle*innerRadius*height);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
            model.getPoint(0).addBoundaryCondition(force);
        }
        if (outerPressure != 0.0)
        {
            RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_X))
                 .add(0.0,-outerPressure*fullAngle*outerRadius*height);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Y)).add(0.0,0.0);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_Z)).add(0.0,0.0);
            model.getPoint(1).addBoundaryCondition(force);
        }
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(500);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(5000);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-14);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getTimeSolver().setEnabled(false);
    return model;
}

//! Build a fixed-free bar whose free end face carries a Roller displacement
//! oriented transversely to the bar.
//!
//! The end face normal points along x, so the boundary condition would restrain
//! the axial direction if the local frame were taken from the geometry. The
//! entered direction points along y, which restrains nothing that is not
//! already zero. The two cases give opposite answers, which is what makes this
//! a usable check of whether the entered direction is honoured.
RModel buildBarWithRoller(double tipForce, bool useEnteredDirection, bool atRoot, const RR3Vector &direction, bool lateralRestraint)
{
    RModel model = buildBar(tipForce,lateralRestraint);

    // Two triangles covering one end face.
    uint slice = atRoot ? 0 : nSlices;
    uint n00 = nodeId(slice,0,0);
    uint n01 = nodeId(slice,0,1);
    uint n10 = nodeId(slice,1,0);
    uint n11 = nodeId(slice,1,1);

    RElement t1(R_ELEMENT_TRI1);
    t1.setNodeId(0,n00);
    t1.setNodeId(1,n10);
    t1.setNodeId(2,n11);
    model.addElement(t1,true,0);

    RElement t2(R_ELEMENT_TRI1);
    t2.setNodeId(0,n00);
    t2.setNodeId(1,n11);
    t2.setNodeId(2,n01);
    model.addElement(t2,true,0);

    RSurface &surface = model.getSurface(0);
    surface.setName("End face");
    surface.setThickness(0.0);

    RBoundaryCondition rollerBc(R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER);
    rollerBc.getComponent(rollerBc.findComponentPosition(R_VARIABLE_DISPLACEMENT)).add(0.0,0.0);
    rollerBc.setExplicitLocalDirection(useEnteredDirection);
    rollerBc.setLocalDirection(direction);
    surface.addBoundaryCondition(rollerBc);

    return model;
}

//! Build a bar loaded by its own weight, with a roller on a mid span cross
//! section whose local direction is tilted away from every global axis.
//!
//! The rotated nodes keep two free degrees of freedom and sit on elements which
//! carry a body force, which is the only situation in which the rotation of the
//! element load vector is visible at all.
RModel buildBarUnderGravity()
{
    RModel model = buildBar(0.0,false);

    REnvironmentCondition gravity(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_X)).add(0.0,0.0);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y)).add(0.0,0.0);
    gravity.getComponent(gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Z)).add(0.0,gravityAcceleration);
    model.getVolume(0).addEnvironmentCondition(gravity);

    return model;
}

RModel buildBarUnderGravityWithTiltedRoller()
{
    RModel model = buildBarUnderGravity();

    // Two triangles across the bar at mid span.
    uint slice = nSlices/2;
    uint n00 = nodeId(slice,0,0);
    uint n01 = nodeId(slice,0,1);
    uint n10 = nodeId(slice,1,0);
    uint n11 = nodeId(slice,1,1);

    RElement t1(R_ELEMENT_TRI1);
    t1.setNodeId(0,n00);
    t1.setNodeId(1,n10);
    t1.setNodeId(2,n11);
    model.addElement(t1,true,0);

    RElement t2(R_ELEMENT_TRI1);
    t2.setNodeId(0,n00);
    t2.setNodeId(1,n11);
    t2.setNodeId(2,n01);
    model.addElement(t2,true,0);

    RSurface &surface = model.getSurface(0);
    surface.setName("Mid span");
    surface.setThickness(0.0);

    // A tilted direction, so that the rotation matrix is genuinely
    // non-symmetric - an axis aligned one would hide a transposed transform.
    RBoundaryCondition rollerBc(R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER);
    rollerBc.getComponent(rollerBc.findComponentPosition(R_VARIABLE_DISPLACEMENT)).add(0.0,0.0);
    rollerBc.setExplicitLocalDirection(true);
    rollerBc.setLocalDirection(RR3Vector(1.0,1.0,1.0));
    surface.addBoundaryCondition(rollerBc);

    return model;
}

//! Run the solver, turning a solver error into a readable test failure.
bool runSolver(RModel &model, QString &errorMessage)
{
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        errorMessage = error.getMessage();
        return false;
    }
    catch (const std::exception &error)
    {
        errorMessage = QString(error.what());
        return false;
    }
    return true;
}

double nodeValue(const RModel &model, RVariableType variableType, uint vectorPosition, uint nodeID)
{
    uint position = model.findVariable(variableType);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(vectorPosition,nodeID);
}

double elementValue(const RModel &model, RVariableType variableType, uint elementID)
{
    uint position = model.findVariable(variableType);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(0,elementID);
}

}

class TestSolverStress : public QObject
{
    Q_OBJECT

private slots:

    void staticAxialBar();
    void staticAdditionalVolumeElements_data();
    void staticAdditionalVolumeElements();
    void vonMisesMatchesComponents();
    void modalAxialBar_data();
    void modalAxialBar();
    void staticTrussBar();
    void staticQuadraticTrussBar();
    void staticBeamBending_data();
    void staticBeamBending();
    void staticQuadraticBeamBending();
    void staticBeamRotatedFrame();
    void staticBeamEndMoment();
    void staticBeamTorsion();
    void prescribedBeamTipRotation();
    void modalBeamBending_data();
    void modalBeamBending();
    void modalTrussBar_data();
    void modalTrussBar();
    void modalSurfaceMassMatrices_data();
    void modalSurfaceMassMatrices();
    void staticSurfaceFormulations_data();
    void staticSurfaceFormulations();
    void quadraticSurfacePatch_data();
    void quadraticSurfacePatch();
    void thermalSurfaceExpansion_data();
    void thermalSurfaceExpansion();
    void axisymmetricThermalExpansion_data();
    void axisymmetricThermalExpansion();
    void axisymmetricCylinderConverges();
    void axisymmetricGeometryIsValidated();
    void surfaceConstraintUsesEnteredDirection();
    void surfaceConstraintFallsBackToTheNormal();
    void explicitDirectionMatchingTheNormalAgrees();
    void reactionsBalanceGravityWithTiltedRoller();
    void prescribedDisplacementIsAppliedInLocalFrame();
    void contradictingConstraintsAreReported();
    void linearConstraintCouplesTrussNodes();
    void contradictoryLinearConstraintsAreReported();
};

//! A fixed-free bar under an axial tip force stretches by F*L/(E*A) and carries
//! a uniform axial stress F/A.
void TestSolverStress::staticAxialBar()
{
    const double tipForce = 1000.0; // [N]

    RModel model = buildBar(tipForce);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    // Switching the X component off must leave the axial direction free.
    double expectedDisplacement = tipForce * barLength / (elasticityModulus * barArea);
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nodeId(nSlices,j,k));
            QVERIFY2(std::fabs(tipDisplacement-expectedDisplacement) < 0.02*expectedDisplacement,
                     qPrintable(QString("Tip displacement = %1 m, expected %2 m").arg(tipDisplacement).arg(expectedDisplacement)));
        }
    }

    // ... while Y and Z stay held at zero everywhere.
    for (uint i=0;i<model.getNNodes();i++)
    {
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,i)) < 1.0e-14);
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,2,i)) < 1.0e-14);
    }

    // Uniform axial stress F/A, checked away from the loaded end.
    double expectedStress = tipForce / barArea;
    uint nChecked = 0;
    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_VOLUME(model.getElement(i).getType()))
        {
            continue;
        }
        RR3Vector center;
        model.getElement(i).findCenter(model.getNodes(),center[0],center[1],center[2]);
        if (center[0] < 0.2*barLength || center[0] > 0.8*barLength)
        {
            continue;
        }

        double stressX = elementValue(model,R_VARIABLE_STRESS_X,i);
        QVERIFY2(std::fabs(stressX-expectedStress) < 0.05*expectedStress,
                 qPrintable(QString("Element %1: stress X = %2 Pa, expected %3 Pa").arg(i).arg(stressX).arg(expectedStress)));

        // Uniaxial tension - von Mises equals the axial stress.
        double vonMises = elementValue(model,R_VARIABLE_STRESS_VON_MISES,i);
        QVERIFY2(std::fabs(vonMises-expectedStress) < 0.05*expectedStress,
                 qPrintable(QString("Element %1: von Mises = %2 Pa, expected %3 Pa").arg(i).arg(vonMises).arg(expectedStress)));

        nChecked++;
    }
    QVERIFY2(nChecked > 0,"No interior element was checked");

    // The reaction at the fixed end balances the applied force.
    double reaction = 0.0;
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            reaction += nodeValue(model,R_VARIABLE_FORCE,0,nodeId(0,j,k));
        }
    }
    QVERIFY2(std::fabs(reaction+tipForce) < 0.02*tipForce,
             qPrintable(QString("Reaction = %1 N, expected %2 N").arg(reaction).arg(-tipForce)));
}

//! WEDGE6 and HEX8 reuse the generic three-dimensional assembly path. This
//! test protects the complete path from interpolation and quadrature through
//! global assembly, result recovery and reactions.
void TestSolverStress::staticAdditionalVolumeElements_data()
{
    QTest::addColumn<int>("elementType");
    QTest::newRow("wedge6") << int(R_ELEMENT_WEDGE1);
    QTest::newRow("hex8") << int(R_ELEMENT_HEXA1);
}

void TestSolverStress::staticAdditionalVolumeElements()
{
    QFETCH(int,elementType);
    const RElementType type = static_cast<RElementType>(elementType);
    const double tipForce = 1000.0; // [N]

    RModel model = buildBar(tipForce,true,type);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedDisplacement = tipForce*barLength/(elasticityModulus*barArea);
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            const double displacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nodeId(nSlices,j,k));
            QVERIFY2(std::fabs(displacement-expectedDisplacement) < 0.01*expectedDisplacement,
                     qPrintable(QString("Tip displacement = %1 m, expected %2 m")
                                .arg(displacement).arg(expectedDisplacement)));
        }
    }

    const double expectedStress = tipForce/barArea;
    uint nChecked = 0;
    for (uint i=0;i<model.getNElements();i++)
    {
        if (model.getElement(i).getType() != type)
        {
            continue;
        }
        const double stressX = elementValue(model,R_VARIABLE_STRESS_X,i);
        QVERIFY2(std::fabs(stressX-expectedStress) < 0.01*expectedStress,
                 qPrintable(QString("Element %1: stress X = %2 Pa, expected %3 Pa")
                            .arg(i).arg(stressX).arg(expectedStress)));
        nChecked++;
    }
    QVERIFY(nChecked > 0);

    double reaction = 0.0;
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            reaction += nodeValue(model,R_VARIABLE_FORCE,0,nodeId(0,j,k));
        }
    }
    QVERIFY2(std::fabs(reaction+tipForce) < 0.01*tipForce,
             qPrintable(QString("Reaction = %1 N, expected %2 N").arg(reaction).arg(-tipForce)));
}

void TestSolverStress::staticSurfaceFormulations_data()
{
    QTest::addColumn<int>("elementType");
    QTest::addColumn<int>("formulation");
    QTest::newRow("tri3-plane-stress")
        << int(R_ELEMENT_TRI1) << int(RElementFormulation::PlaneStress);
    QTest::newRow("tri3-plane-strain")
        << int(R_ELEMENT_TRI1) << int(RElementFormulation::PlaneStrain);
    QTest::newRow("quad4-plane-stress")
        << int(R_ELEMENT_QUAD1) << int(RElementFormulation::PlaneStress);
    QTest::newRow("quad4-plane-strain")
        << int(R_ELEMENT_QUAD1) << int(RElementFormulation::PlaneStrain);
}

//! TRI3 and QUAD4 use the formulation stored on their surface group. With
//! epsilon-y restrained, both the axial stiffness and sigma-y distinguish the
//! plane-stress and plane-strain constitutive matrices unambiguously.
void TestSolverStress::staticSurfaceFormulations()
{
    QFETCH(int,elementType);
    QFETCH(int,formulation);
    constexpr double tipForce = 1000.0;
    constexpr double width = 0.1;
    constexpr double thickness = 0.02;
    constexpr double v = 0.25;

    const RElementFormulation selected = static_cast<RElementFormulation>(formulation);
    RModel model = buildSurfaceStrip(tipForce,static_cast<RElementType>(elementType),selected);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double d00 = selected == RElementFormulation::PlaneStress
                     ? elasticityModulus/(1.0-v*v)
                     : elasticityModulus*(1.0-v)/((1.0+v)*(1.0-2.0*v));
    const double expectedDisplacement = tipForce*barLength/(d00*width*thickness);
    const uint tip0 = model.getNNodes()-2;
    const uint tip1 = model.getNNodes()-1;
    const double tipDisplacement = 0.5*(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,tip0)
                                      + nodeValue(model,R_VARIABLE_DISPLACEMENT,0,tip1));
    QVERIFY2(std::fabs(tipDisplacement-expectedDisplacement) < 0.02*expectedDisplacement,
             qPrintable(QString("Tip displacement = %1 m, expected %2 m")
                        .arg(tipDisplacement).arg(expectedDisplacement)));

    const double expectedStressX = tipForce/(width*thickness);
    const double expectedStressY = expectedStressX
                                 * (selected == RElementFormulation::PlaneStress ? v : v/(1.0-v));
    uint checked = 0;
    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_SURFACE(model.getElement(i).getType()))
        {
            continue;
        }
        const double stressX = elementValue(model,R_VARIABLE_STRESS_X,i);
        const double stressY = elementValue(model,R_VARIABLE_STRESS_Y,i);
        QVERIFY2(std::fabs(stressX-expectedStressX) < 0.03*expectedStressX,
                 qPrintable(QString("Element %1: stress X = %2 Pa, expected %3 Pa")
                            .arg(i).arg(stressX).arg(expectedStressX)));
        QVERIFY2(std::fabs(stressY-expectedStressY) < 0.03*expectedStressX,
                 qPrintable(QString("Element %1: stress Y = %2 Pa, expected %3 Pa")
                            .arg(i).arg(stressY).arg(expectedStressY)));
        checked++;
    }
    QVERIFY(checked > 0);

    const double reaction = nodeValue(model,R_VARIABLE_FORCE,0,0)
                          + nodeValue(model,R_VARIABLE_FORCE,0,1);
    QVERIFY2(std::fabs(reaction+tipForce) < 0.02*tipForce,
             qPrintable(QString("Reaction = %1 N, expected %2 N").arg(reaction).arg(-tipForce)));
}

void TestSolverStress::axisymmetricThermalExpansion_data()
{
    QTest::addColumn<int>("elementType");
    QTest::newRow("tri3") << int(R_ELEMENT_TRI1);
    QTest::newRow("quad4") << int(R_ELEMENT_QUAD1);
}

//! Uniform thermal expansion is an exact axisymmetric patch: radial, axial and
//! hoop strains are all alpha*dT and all four stress components must vanish.
void TestSolverStress::axisymmetricThermalExpansion()
{
    QFETCH(int,elementType);
    constexpr uint radialSlices = 6;
    constexpr double alpha = 1.2e-5;
    constexpr double deltaT = 50.0;
    RModel model = buildAxisymmetricAnnulus(radialSlices,
                                             static_cast<RElementType>(elementType),
                                             false,alpha,293.15-deltaT);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    for (uint i=0;i<2*(radialSlices+1);i++)
    {
        const double expectedRadial = alpha*deltaT*model.getNode(i).getX();
        const double expectedAxial = alpha*deltaT*model.getNode(i).getY();
        const double radial = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,i);
        const double axial = nodeValue(model,R_VARIABLE_DISPLACEMENT,1,i);
        QVERIFY2(std::fabs(radial-expectedRadial) < 1.0e-5*alpha*deltaT*2.0,
                 qPrintable(QString("Node %1 radial displacement = %2, expected %3")
                            .arg(i).arg(radial).arg(expectedRadial)));
        QVERIFY2(std::fabs(axial-expectedAxial) < 1.0e-5*alpha*deltaT,
                 qPrintable(QString("Node %1 axial displacement = %2, expected %3")
                            .arg(i).arg(axial).arg(expectedAxial)));
    }

    const double stressScale = elasticityModulus*alpha*deltaT;
    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_SURFACE(model.getElement(i).getType()))
        {
            continue;
        }
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_X,i)) < 1.0e-5*stressScale);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_Y,i)) < 1.0e-5*stressScale);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_Z,i)) < 1.0e-5*stressScale);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_XY,i)) < 1.0e-5*stressScale);
    }

    const double forceTolerance = 1.0e-6*stressScale*3.14159265358979323846*(4.0-1.0);
    for (uint i=0;i<2*(radialSlices+1);i++)
    {
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,0,i)) < forceTolerance);
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,1,i)) < forceTolerance);
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,2,i)) < forceTolerance);
    }
}

//! A long thick cylinder under internal pressure has the Lame solution. Mesh
//! refinement must reduce both radial-displacement and recovered-stress error.
void TestSolverStress::axisymmetricCylinderConverges()
{
    constexpr double innerRadius = 1.0;
    constexpr double outerRadius = 2.0;
    constexpr double pressure = 1.0e6;
    constexpr double v = 0.25;
    const double A = pressure*innerRadius*innerRadius
                   / (outerRadius*outerRadius-innerRadius*innerRadius);
    const double B = pressure*innerRadius*innerRadius*outerRadius*outerRadius
                   / (outerRadius*outerRadius-innerRadius*innerRadius);

    double previousDisplacementError = std::numeric_limits<double>::infinity();
    double previousStressError = std::numeric_limits<double>::infinity();
    for (uint radialSlices : {2u,4u,8u,16u})
    {
        RModel model = buildAxisymmetricAnnulus(radialSlices,R_ELEMENT_QUAD1,true,
                                                 0.0,293.15,pressure,0.0);
        model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
        QString solverError;
        QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

        double displacementError = 0.0;
        for (uint i=0;i<=radialSlices;i++)
        {
            const double radius = model.getNode(2*i).getX();
            const double expected = (1.0+v)/elasticityModulus
                                  * ((1.0-2.0*v)*A*radius+B/radius);
            const double value = 0.5*(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,2*i)
                                      + nodeValue(model,R_VARIABLE_DISPLACEMENT,0,2*i+1));
            displacementError = std::max(displacementError,std::fabs(value-expected));
        }

        double stressError = 0.0;
        for (uint i=0;i<radialSlices;i++)
        {
            const double radius = 0.5*(model.getNode(2*i).getX()+model.getNode(2*(i+1)).getX());
            const double expectedRadial = A-B/(radius*radius);
            const double expectedHoop = A+B/(radius*radius);
            const double expectedAxial = 2.0*v*A;
            stressError = std::max(stressError,
                                   std::fabs(elementValue(model,R_VARIABLE_STRESS_X,i)-expectedRadial));
            stressError = std::max(stressError,
                                   std::fabs(elementValue(model,R_VARIABLE_STRESS_Z,i)-expectedHoop));
            stressError = std::max(stressError,
                                   std::fabs(elementValue(model,R_VARIABLE_STRESS_Y,i)-expectedAxial));
        }

        QVERIFY2(displacementError < previousDisplacementError,
                 qPrintable(QString("Displacement error did not decrease at %1 slices: %2 >= %3")
                            .arg(radialSlices).arg(displacementError).arg(previousDisplacementError)));
        QVERIFY2(stressError < previousStressError,
                 qPrintable(QString("Stress error did not decrease at %1 slices: %2 >= %3")
                            .arg(radialSlices).arg(stressError).arg(previousStressError)));
        previousDisplacementError = displacementError;
        previousStressError = stressError;
    }

    QVERIFY(previousDisplacementError < 2.0e-3*pressure/elasticityModulus);
    QVERIFY(previousStressError < 0.03*pressure);
}

void TestSolverStress::axisymmetricGeometryIsValidated()
{
    RModel model = buildAxisymmetricAnnulus(2,R_ELEMENT_QUAD1,true);
    model.getNode(0).setZ(0.1);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(!runSolver(model,solverError),"Out-of-plane axisymmetric geometry was accepted");
    QVERIFY2(solverError.contains("meridian plane"),qPrintable(solverError));
}

void TestSolverStress::thermalSurfaceExpansion_data()
{
    QTest::addColumn<int>("formulation");
    QTest::newRow("plane-stress") << int(RElementFormulation::PlaneStress);
    QTest::newRow("plane-strain") << int(RElementFormulation::PlaneStrain);
}

//! A free axial edge with restrained transverse strain expands according to
//! the selected 2-D constitutive law. This protects both the equivalent nodal
//! thermal force and the stress-recovery subtraction from thickness leakage.
void TestSolverStress::thermalSurfaceExpansion()
{
    QFETCH(int,formulation);
    constexpr double alpha = 1.2e-5;
    constexpr double deltaT = 50.0;
    constexpr double v = 0.25;
    const RElementFormulation selected = static_cast<RElementFormulation>(formulation);
    RModel model = buildSurfaceStrip(0.0,R_ELEMENT_QUAD1,selected,alpha,293.15-deltaT);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedStrain = selected == RElementFormulation::PlaneStress
                                ? alpha*deltaT*(1.0+v)
                                : alpha*deltaT/(1.0-v);
    const double expectedTip = expectedStrain*barLength;
    const uint tip0 = model.getNNodes()-2;
    const uint tip1 = model.getNNodes()-1;
    const double tip = 0.5*(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,tip0)
                           + nodeValue(model,R_VARIABLE_DISPLACEMENT,0,tip1));
    QVERIFY2(std::fabs(tip-expectedTip) < 0.02*expectedTip,
             qPrintable(QString("Thermal tip displacement = %1 m, expected %2 m")
                        .arg(tip).arg(expectedTip)));

    const double expectedStressY = selected == RElementFormulation::PlaneStress
                                 ? -elasticityModulus*alpha*deltaT
                                 : -elasticityModulus*alpha*deltaT/(1.0-v*v);
    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_SURFACE(model.getElement(i).getType()))
        {
            continue;
        }
        const double stressX = elementValue(model,R_VARIABLE_STRESS_X,i);
        const double stressY = elementValue(model,R_VARIABLE_STRESS_Y,i);
        QVERIFY2(std::fabs(stressX) < 0.02*std::fabs(expectedStressY),
                 qPrintable(QString("Element %1: free thermal stress X = %2 Pa").arg(i).arg(stressX)));
        QVERIFY2(std::fabs(stressY-expectedStressY) < 0.02*std::fabs(expectedStressY),
                 qPrintable(QString("Element %1: thermal stress Y = %2 Pa, expected %3 Pa")
                            .arg(i).arg(stressY).arg(expectedStressY)));
    }
}

//! The reported von Mises stress has to be the classical invariant of the
//! stored stress components, which combines the normal and the shear part in
//! quadrature rather than by adding them.
void TestSolverStress::vonMisesMatchesComponents()
{
    RModel model = buildBar(1000.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    uint nShearElements = 0;

    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_VOLUME(model.getElement(i).getType()))
        {
            continue;
        }

        double sx = elementValue(model,R_VARIABLE_STRESS_X,i);
        double sy = elementValue(model,R_VARIABLE_STRESS_Y,i);
        double sz = elementValue(model,R_VARIABLE_STRESS_Z,i);
        double tyz = elementValue(model,R_VARIABLE_STRESS_YZ,i);
        double txz = elementValue(model,R_VARIABLE_STRESS_XZ,i);
        double txy = elementValue(model,R_VARIABLE_STRESS_XY,i);

        double expected = std::sqrt(sx*sx + sy*sy + sz*sz - (sx*sy + sy*sz + sz*sx)
                                    + 3.0*(tyz*tyz + txz*txz + txy*txy));
        double vonMises = elementValue(model,R_VARIABLE_STRESS_VON_MISES,i);

        QVERIFY2(std::fabs(vonMises-expected) < 1.0e-6*std::max(1.0,std::fabs(expected)),
                 qPrintable(QString("Element %1: von Mises = %2 Pa, invariant of the components = %3 Pa")
                            .arg(i).arg(vonMises).arg(expected)));

        if (std::fabs(tyz) + std::fabs(txz) + std::fabs(txy) > 1.0e-6*std::fabs(sx))
        {
            nShearElements++;
        }
    }

    // The check only discriminates the two formulas where a shear component is
    // present, so make sure at least some elements carry one.
    QVERIFY2(nShearElements > 0,"No element carried a shear stress - the check would not discriminate");
}

//! The fundamental axial mode of a fixed-free bar is at c/(4*L) with
//! c = sqrt(E/rho). The modal setup has to report it in Hz.
void TestSolverStress::modalAxialBar_data()
{
    QTest::addColumn<int>("massMatrixType");
    QTest::addColumn<int>("expectedBias");
    QTest::newRow("consistent-mass") << int(RMassMatrixType::Consistent) << 1;
    QTest::newRow("lumped-mass") << int(RMassMatrixType::Lumped) << -1;
}

void TestSolverStress::modalAxialBar()
{
    QFETCH(int,massMatrixType);
    QFETCH(int,expectedBias);
    RModel model = buildBar(0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.getProblemSetup().setMassMatrixType(
        static_cast<RMassMatrixType>(massMatrixType));

    RModalSetup &modalSetup = model.getProblemSetup().getModalSetup();
    modalSetup.setMethod(R_MODAL_MULTIPLE_MODES);
    modalSetup.setNModesToExtract(4);
    modalSetup.setNIterations(300);
    modalSetup.setConvergenceValue(1.0e-12);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    // Modes are processed from the highest index down, so when the run is over
    // the modal setup holds the fundamental one. For a fixed-free bar that is
    // the first axial mode at c/(4*L) with c = sqrt(E/rho).
    double soundSpeed = std::sqrt(elasticityModulus/density);
    double expectedFrequency = soundSpeed / (4.0 * barLength);
    double frequency = modalSetup.getFrequency();

    QVERIFY2(std::isfinite(frequency) && frequency > 0.0,
             qPrintable(QString("Frequency = %1 Hz is not a usable value").arg(frequency)));
    QVERIFY2(std::fabs(frequency-expectedFrequency) < 0.02*expectedFrequency,
             qPrintable(QString("Fundamental frequency = %1 Hz, expected %2 Hz").arg(frequency).arg(expectedFrequency)));

    // Every extracted mode is stored as a displacement field.
    uint displacementPosition = model.findVariable(R_VARIABLE_DISPLACEMENT);
    QVERIFY(displacementPosition != RConstants::eod);

    double maxDisplacement = 0.0;
    for (uint i=0;i<model.getNNodes();i++)
    {
        maxDisplacement = std::max(maxDisplacement,std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,i)));
    }
    QVERIFY2(maxDisplacement > 0.0,"The mode shape is identically zero");
}

//! The same bar built from truss elements has to give the same closed form
//! answer - the axial stiffness must carry the integration weight, and the
//! reported axial stress must be a stress and not an axial force.
void TestSolverStress::staticTrussBar()
{
    const double tipForce = 1000.0; // [N]

    RModel model = buildTruss(tipForce);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double expectedDisplacement = tipForce * barLength / (elasticityModulus * barArea);
    double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nSlices);
    QVERIFY2(std::fabs(tipDisplacement-expectedDisplacement) < 0.01*expectedDisplacement,
             qPrintable(QString("Tip displacement = %1 m, expected %2 m").arg(tipDisplacement).arg(expectedDisplacement)));

    // The displacement grows linearly along a uniformly loaded bar.
    for (uint i=0;i<=nSlices;i++)
    {
        double expected = expectedDisplacement * double(i) / double(nSlices);
        double value = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,i);
        QVERIFY2(std::fabs(value-expected) < 0.01*expectedDisplacement,
                 qPrintable(QString("Node %1: displacement = %2 m, expected %3 m").arg(i).arg(value).arg(expected)));
    }

    double expectedStress = tipForce / barArea;
    for (uint i=0;i<nSlices;i++)
    {
        double normalStress = elementValue(model,R_VARIABLE_STRESS_NORMAL,i);
        QVERIFY2(std::fabs(normalStress-expectedStress) < 0.01*expectedStress,
                 qPrintable(QString("Element %1: normal stress = %2 Pa, expected %3 Pa").arg(i).arg(normalStress).arg(expectedStress)));
    }

    double reaction = nodeValue(model,R_VARIABLE_FORCE,0,0);
    QVERIFY2(std::fabs(reaction+tipForce) < 0.01*tipForce,
             qPrintable(QString("Reaction = %1 N, expected %2 N").arg(reaction).arg(-tipForce)));
}

//! A truss bar has the same axial modes as the solid one.
void TestSolverStress::modalTrussBar_data()
{
    QTest::addColumn<int>("massMatrixType");
    QTest::addColumn<int>("expectedBias");
    QTest::newRow("consistent-mass") << int(RMassMatrixType::Consistent) << 1;
    QTest::newRow("lumped-mass") << int(RMassMatrixType::Lumped) << -1;
}

void TestSolverStress::staticQuadraticTrussBar()
{
    constexpr double tipForce = 1000.0;
    RModel model = buildQuadraticTruss(tipForce);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedDisplacement = tipForce*barLength/(elasticityModulus*barArea);
    const uint tipNode = model.getNNodes()-1;
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,tipNode)-expectedDisplacement)
            < 1.0e-5*expectedDisplacement);
    for (uint i=0;i<model.getNNodes();i++)
    {
        const double expected = expectedDisplacement*model.getNode(i).getX()/barLength;
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,i)-expected)
                < 1.0e-5*expectedDisplacement);
    }
    const double expectedStress = tipForce/barArea;
    for (uint i=0;i<4;i++)
    {
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,i)-expectedStress)
                < 1.0e-5*expectedStress);
    }
}

void TestSolverStress::quadraticSurfacePatch_data()
{
    QTest::addColumn<int>("elementType");
    QTest::newRow("tri6") << int(R_ELEMENT_TRI2);
    QTest::newRow("quad9") << int(R_ELEMENT_QUAD2);
}

void TestSolverStress::quadraticSurfacePatch()
{
    QFETCH(int,elementType);
    constexpr double strain = 1.0e-4;
    RModel model = buildQuadraticSurfacePatch(static_cast<RElementType>(elementType),strain);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    constexpr double v = 0.25;
    const double expectedStressX = elasticityModulus*strain/(1.0-v*v);
    const double expectedStressY = v*expectedStressX;
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_X,0)-expectedStressX)
            < 1.0e-8*expectedStressX);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_Y,0)-expectedStressY)
            < 1.0e-8*expectedStressX);
    for (uint node=0;node<model.getNNodes();node++)
    {
        const double expected = strain*model.getNode(node).getX();
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,node)-expected)
                < 1.0e-8*strain);
    }
}

//! A tip-loaded Euler-Bernoulli cantilever follows the exact cubic deflection
//! curve. Distinct Iy and Iz values verify the section-axis convention.
void TestSolverStress::staticBeamBending_data()
{
    QTest::addColumn<int>("component");
    QTest::addColumn<double>("tipForce");
    QTest::addColumn<double>("secondMoment");
    QTest::newRow("local-y-bending-about-z") << 1 << 1200.0 << 5.0e-5;
    QTest::newRow("local-z-bending-about-y") << 2 << -900.0 << 2.0e-5;
}

void TestSolverStress::staticBeamBending()
{
    QFETCH(int,component);
    QFETCH(double,tipForce);
    QFETCH(double,secondMoment);
    RModel model = buildBeam(component == 1 ? tipForce : 0.0,
                             component == 2 ? tipForce : 0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const uint beamSlices = 6;
    const double expectedTip = tipForce*barLength*barLength*barLength
                             /(3.0*elasticityModulus*secondMoment);
    const double computedTip = nodeValue(model,R_VARIABLE_DISPLACEMENT,uint(component),beamSlices);
    QVERIFY2(std::fabs(computedTip-expectedTip) < 1.0e-5*std::fabs(expectedTip),
             qPrintable(QString("Beam tip displacement = %1 m, expected %2 m")
                        .arg(computedTip).arg(expectedTip)));

    for (uint i=0;i<=beamSlices;i++)
    {
        const double x = barLength*double(i)/double(beamSlices);
        const double expected = tipForce*x*x*(3.0*barLength-x)
                              /(6.0*elasticityModulus*secondMoment);
        const double value = nodeValue(model,R_VARIABLE_DISPLACEMENT,uint(component),i);
        QVERIFY2(std::fabs(value-expected) < 1.0e-5*std::fabs(expectedTip),
                 qPrintable(QString("Beam node %1 displacement = %2 m, expected %3 m")
                            .arg(i).arg(value).arg(expected)));
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,i)) < 1.0e-12);
        const uint other = component == 1 ? 2u : 1u;
        QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,other,i)) < 1.0e-12);
    }

    const double reaction = nodeValue(model,R_VARIABLE_FORCE,uint(component),0);
    QVERIFY2(std::fabs(reaction+tipForce) < 1.0e-5*std::fabs(tipForce),
             qPrintable(QString("Beam root reaction = %1 N, expected %2 N")
                        .arg(reaction).arg(-tipForce)));
}

void TestSolverStress::staticQuadraticBeamBending()
{
    constexpr double tipForce = 1200.0;
    constexpr double iz = 5.0e-5;
    RModel model = buildQuadraticBeam(tipForce);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedTip = tipForce*std::pow(barLength,3)/(3.0*elasticityModulus*iz);
    const double expectedMiddle = tipForce*std::pow(0.5*barLength,2)*(3.0*barLength-0.5*barLength)
                                /(6.0*elasticityModulus*iz);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,2)-expectedTip)
            < 1.0e-5*expectedTip);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,1)-expectedMiddle)
            < 1.0e-5*expectedTip);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,1,0)+tipForce) < 1.0e-5*tipForce);
    const double expectedStress = tipForce*barLength/model.getLine(0).getSectionModulusZ();
    const double recoveredStress = elementValue(model,R_VARIABLE_STRESS_NORMAL,0);
    QVERIFY2(std::fabs(recoveredStress-expectedStress) < 1.0e-5*expectedStress,
             qPrintable(QString("BEAM2 normal stress = %1 Pa, expected %2 Pa; root moment=%3, ux=%4, uy=%5, uz=%6")
                        .arg(recoveredStress).arg(expectedStress)
                        .arg(nodeValue(model,R_VARIABLE_MOMENT,2,0))
                        .arg(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,2))
                        .arg(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,2))
                        .arg(nodeValue(model,R_VARIABLE_DISPLACEMENT,2,2))));
}

//! A beam aligned with global Z exercises the local-to-global frame transform
//! while retaining the same local-Y bending reference solution.
void TestSolverStress::staticBeamRotatedFrame()
{
    constexpr uint beamSlices = 6;
    constexpr double tipForce = 1200.0;
    constexpr double iz = 5.0e-5;
    RModel model = buildBeam(tipForce,0.0);
    for (uint i=0;i<=beamSlices;i++)
    {
        model.getNode(i).set(0.0,0.0,barLength*double(i)/double(beamSlices));
    }
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedTip = tipForce*std::pow(barLength,3)
                             /(3.0*elasticityModulus*iz);
    const double computedTip = nodeValue(model,R_VARIABLE_DISPLACEMENT,1,beamSlices);
    QVERIFY2(std::fabs(computedTip-expectedTip) < 1.0e-5*expectedTip,
             qPrintable(QString("Rotated beam tip displacement = %1 m, expected %2 m")
                        .arg(computedTip).arg(expectedTip)));
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,beamSlices)) < 1.0e-12);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,2,beamSlices)) < 1.0e-12);
    QVERIFY2(std::fabs(nodeValue(model,R_VARIABLE_FORCE,1,0)+tipForce) < 1.0e-5*tipForce,
             "Rotated beam root reaction does not balance the tip load");
}

//! A point moment must enter the rotational equation directly. The explicit
//! root Rotation condition replaces the historical translation-only clamp
//! fallback and the stored nodal moments must balance the applied end moment.
void TestSolverStress::staticBeamEndMoment()
{
    constexpr uint beamSlices = 6;
    constexpr double endMoment = 750.0;
    constexpr double iz = 5.0e-5;
    RModel model = buildBeam(0.0,0.0);

    RBoundaryCondition fixedRotation(R_BOUNDARY_CONDITION_ROTATION);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixedRotation);

    RBoundaryCondition moment(R_BOUNDARY_CONDITION_MOMENT);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_X)).add(0.0,0.0);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_Y)).add(0.0,0.0);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_Z)).add(0.0,endMoment);
    model.getPoint(1).addBoundaryCondition(moment);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedRotation = endMoment*barLength/(elasticityModulus*iz);
    const double expectedDisplacement = endMoment*barLength*barLength/(2.0*elasticityModulus*iz);
    QCOMPARE(nodeValue(model,R_VARIABLE_ROTATION,0,0),0.0);
    QCOMPARE(nodeValue(model,R_VARIABLE_ROTATION,1,0),0.0);
    QCOMPARE(nodeValue(model,R_VARIABLE_ROTATION,2,0),0.0);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_ROTATION,2,beamSlices)-expectedRotation)
            < 1.0e-5*expectedRotation);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,beamSlices)-expectedDisplacement)
            < 1.0e-5*expectedDisplacement);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_MOMENT,2,0)+endMoment)
            < 1.0e-5*endMoment);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_MOMENT,2,beamSlices)-endMoment)
            < 1.0e-5*endMoment);
    const double expectedBendingStress = endMoment/model.getLine(0).getSectionModulusZ();
    for (uint i=0;i<beamSlices;i++)
    {
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,i)-expectedBendingStress)
                < 1.0e-5*expectedBendingStress);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_VON_MISES,i)-expectedBendingStress)
                < 1.0e-5*expectedBendingStress);
    }
}

//! Pure torsion verifies the torsional section modulus and the line-element
//! stress-invariant convention used by the generic result variables.
void TestSolverStress::staticBeamTorsion()
{
    constexpr uint beamSlices = 6;
    constexpr double torque = 600.0;
    RModel model = buildBeam(0.0,0.0);

    RBoundaryCondition fixedRotation(R_BOUNDARY_CONDITION_ROTATION);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixedRotation);

    RBoundaryCondition moment(R_BOUNDARY_CONDITION_MOMENT);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_X)).add(0.0,torque);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_Y)).add(0.0,0.0);
    moment.getComponent(moment.findComponentPosition(R_VARIABLE_MOMENT_Z)).add(0.0,0.0);
    model.getPoint(1).addBoundaryCondition(moment);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const RLine &line = model.getLine(0);
    const double shearModulus = elasticityModulus/(2.0*(1.0+0.3));
    const double expectedRotation = torque*barLength/(shearModulus*line.getTorsionalConstant());
    const double expectedShearInvariant = std::sqrt(3.0)*torque/line.getTorsionalSectionModulus();
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_ROTATION,0,beamSlices)-expectedRotation)
            < 1.0e-5*expectedRotation);
    for (uint i=0;i<beamSlices;i++)
    {
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,i)) < 1.0e-5);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_SHEAR,i)-expectedShearInvariant)
                < 1.0e-5*expectedShearInvariant);
        QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_VON_MISES,i)-expectedShearInvariant)
                < 1.0e-5*expectedShearInvariant);
    }
}

//! Non-zero prescribed rotation exercises the elimination term, rather than
//! merely removing an equation. With a free transverse tip the exact solution
//! has v(L)=L*theta/2 and equal/opposite end moments E*I*theta/L.
void TestSolverStress::prescribedBeamTipRotation()
{
    constexpr uint beamSlices = 6;
    constexpr double prescribedRotation = 2.0e-3;
    constexpr double iz = 5.0e-5;
    RModel model = buildBeam(0.0,0.0);

    RBoundaryCondition fixedRotation(R_BOUNDARY_CONDITION_ROTATION);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    fixedRotation.getComponent(fixedRotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(fixedRotation);

    RBoundaryCondition tipRotation(R_BOUNDARY_CONDITION_ROTATION);
    tipRotation.getComponent(tipRotation.findComponentPosition(R_VARIABLE_ROTATION_X)).setEnabled(false);
    tipRotation.getComponent(tipRotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).setEnabled(false);
    tipRotation.getComponent(tipRotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,prescribedRotation);
    model.getPoint(1).addBoundaryCondition(tipRotation);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    const double expectedDisplacement = barLength*prescribedRotation/2.0;
    const double expectedMoment = elasticityModulus*iz*prescribedRotation/barLength;
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_ROTATION,2,beamSlices)-prescribedRotation)
            < 1.0e-10*prescribedRotation);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,beamSlices)-expectedDisplacement)
            < 1.0e-5*expectedDisplacement);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_MOMENT,2,0)+expectedMoment)
            < 1.0e-5*expectedMoment);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_MOMENT,2,beamSlices)-expectedMoment)
            < 1.0e-5*expectedMoment);
}

//! The first cantilever mode checks the BEAM1 mass matrix together with the
//! rotational DOFs. The weaker Iy axis must govern the fundamental frequency.
void TestSolverStress::modalBeamBending_data()
{
    QTest::addColumn<int>("massMatrixType");
    QTest::newRow("consistent-mass") << int(RMassMatrixType::Consistent);
    QTest::newRow("lumped-mass") << int(RMassMatrixType::Lumped);
}

void TestSolverStress::modalBeamBending()
{
    QFETCH(int,massMatrixType);
    RModel model = buildBeam(0.0,0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.getProblemSetup().setMassMatrixType(
        static_cast<RMassMatrixType>(massMatrixType));

    RModalSetup &modalSetup = model.getProblemSetup().getModalSetup();
    modalSetup.setMethod(R_MODAL_MULTIPLE_MODES);
    modalSetup.setNModesToExtract(4);
    modalSetup.setNIterations(400);
    modalSetup.setConvergenceValue(1.0e-11);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    constexpr double beta1 = 1.875104068711961;
    constexpr double weakSecondMoment = 2.0e-5;
    const double expectedFrequency = beta1*beta1/(2.0*std::acos(-1.0))
                                   * std::sqrt(elasticityModulus*weakSecondMoment
                                               /(density*barArea*std::pow(barLength,4)));
    const double frequency = modalSetup.getFrequency();
    QVERIFY2(std::isfinite(frequency) && frequency > 0.0,
             qPrintable(QString("Beam frequency = %1 Hz is not usable").arg(frequency)));
    QVERIFY2(std::fabs(frequency-expectedFrequency) < 0.04*expectedFrequency,
             qPrintable(QString("Beam fundamental frequency = %1 Hz, expected %2 Hz")
                        .arg(frequency).arg(expectedFrequency)));
}

void TestSolverStress::modalTrussBar()
{
    QFETCH(int,massMatrixType);
    QFETCH(int,expectedBias);
    RModel model = buildTruss(0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.getProblemSetup().setMassMatrixType(
        static_cast<RMassMatrixType>(massMatrixType));

    RModalSetup &modalSetup = model.getProblemSetup().getModalSetup();
    modalSetup.setMethod(R_MODAL_MULTIPLE_MODES);
    modalSetup.setNModesToExtract(3);
    modalSetup.setNIterations(300);
    modalSetup.setConvergenceValue(1.0e-12);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double soundSpeed = std::sqrt(elasticityModulus/density);
    double expectedFrequency = soundSpeed / (4.0 * barLength);
    double frequency = modalSetup.getFrequency();

    QVERIFY2(std::fabs(frequency-expectedFrequency) < 0.02*expectedFrequency,
             qPrintable(QString("Fundamental frequency = %1 Hz, expected %2 Hz").arg(frequency).arg(expectedFrequency)));
    QVERIFY2(expectedBias*(frequency-expectedFrequency) > 0.0,
             qPrintable(QString("Frequency = %1 Hz does not show the expected mass bias around %2 Hz")
                        .arg(frequency).arg(expectedFrequency)));
    QVERIFY2(expectedBias*(frequency-expectedFrequency) > 0.0,
             qPrintable(QString("Frequency = %1 Hz does not show the expected %2 mass bias around %3 Hz")
                        .arg(frequency)
                        .arg(expectedBias > 0 ? "positive consistent" : "negative lumped")
                        .arg(expectedFrequency)));
}

//! Membrane modes exercise the local two-dimensional mass matrix before it is
//! rotated into the three-dimensional model frame. Both supported first-order
//! topologies must retain the analytical fixed-free axial mode.
void TestSolverStress::modalSurfaceMassMatrices_data()
{
    QTest::addColumn<int>("elementType");
    QTest::addColumn<int>("massMatrixType");
    QTest::addColumn<int>("expectedBias");
    QTest::newRow("tri-consistent") << int(R_ELEMENT_TRI1) << int(RMassMatrixType::Consistent) << 1;
    QTest::newRow("tri-lumped") << int(R_ELEMENT_TRI1) << int(RMassMatrixType::Lumped) << -1;
    QTest::newRow("quad-consistent") << int(R_ELEMENT_QUAD1) << int(RMassMatrixType::Consistent) << 1;
    QTest::newRow("quad-lumped") << int(R_ELEMENT_QUAD1) << int(RMassMatrixType::Lumped) << -1;
}

void TestSolverStress::modalSurfaceMassMatrices()
{
    QFETCH(int,elementType);
    QFETCH(int,massMatrixType);
    QFETCH(int,expectedBias);

    RModel model = buildSurfaceStrip(0.0,
                                     static_cast<RElementType>(elementType),
                                     RElementFormulation::PlaneStress);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.getProblemSetup().setMassMatrixType(
        static_cast<RMassMatrixType>(massMatrixType));

    RModalSetup &modalSetup = model.getProblemSetup().getModalSetup();
    modalSetup.setMethod(R_MODAL_MULTIPLE_MODES);
    modalSetup.setNModesToExtract(3);
    modalSetup.setNIterations(300);
    modalSetup.setConvergenceValue(1.0e-12);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    constexpr double surfacePoissonRatio = 0.25;
    const double effectiveModulus = elasticityModulus
                                  / (1.0-surfacePoissonRatio*surfacePoissonRatio);
    const double expectedFrequency = std::sqrt(effectiveModulus/density)
                                   / (4.0*barLength);
    const double frequency = modalSetup.getFrequency();

    QVERIFY2(std::isfinite(frequency) && frequency > 0.0,
             qPrintable(QString("Frequency = %1 Hz is not a usable value").arg(frequency)));
    QVERIFY2(std::fabs(frequency-expectedFrequency) < 0.02*expectedFrequency,
             qPrintable(QString("Fundamental frequency = %1 Hz, expected %2 Hz")
                        .arg(frequency).arg(expectedFrequency)));
    QVERIFY2(expectedBias*(frequency-expectedFrequency) > 0.0,
             qPrintable(QString("Frequency = %1 Hz does not show the expected mass bias around %2 Hz")
                        .arg(frequency).arg(expectedFrequency)));
}

//! A roller on the free end face restrains the direction of its local frame.
//!
//! The face normal points along the bar, so falling back to the geometry clamps
//! the axial direction and nothing moves. Entering a transverse direction has
//! to restrain that direction instead and leave the bar free to stretch.
//!
//! The global lateral restraint has to be left out here, because a globally
//! phrased constraint on a node which carries a local frame is read in that
//! frame and would lock the tip completely. Without it this slender one element
//! thick bar is soft in bending, so only the presence or absence of the axial
//! restraint is compared, never the magnitude.
void TestSolverStress::surfaceConstraintUsesEnteredDirection()
{
    const double tipForce = 1000.0; // [N]
    double freeDisplacement = tipForce * barLength / (elasticityModulus * barArea);

    RModel model = buildBarWithRoller(tipForce,true,false,RR3Vector(0.0,1.0,0.0),false);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nodeId(nSlices,0,0));

    QVERIFY2(std::fabs(tipDisplacement) > 0.5*freeDisplacement,
             qPrintable(QString("Tip displacement = %1 m - the axial direction is still restrained, "
                                "so the entered local direction was ignored").arg(tipDisplacement)));
}

//! With the entered direction switched off, the local frame comes from the
//! element normals as before. The end face normal points along the bar, so the
//! same boundary condition now clamps the free end and nothing moves.
void TestSolverStress::surfaceConstraintFallsBackToTheNormal()
{
    const double tipForce = 1000.0; // [N]

    RModel model = buildBarWithRoller(tipForce,false,false,RR3Vector(0.0,1.0,0.0),false);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double freeDisplacement = tipForce * barLength / (elasticityModulus * barArea);
    double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nodeId(nSlices,0,0));

    QVERIFY2(std::fabs(tipDisplacement) < 0.02*freeDisplacement,
             qPrintable(QString("Tip displacement = %1 m, expected the end to be clamped by its normal")
                        .arg(tipDisplacement)));
}

//! Entering the direction the geometry would have produced anyway has to leave
//! the answer untouched, which shows the entered direction feeds the local
//! rotation machinery in the same way and not just differently.
//!
//! The roller sits on the root face, which is fully fixed already, and names
//! the direction its own normal would have given. It therefore adds nothing,
//! and the bar still has to stretch by exactly F*L/(E*A).
void TestSolverStress::explicitDirectionMatchingTheNormalAgrees()
{
    const double tipForce = 1000.0; // [N]

    RModel model = buildBarWithRoller(tipForce,true,true,RR3Vector(1.0,0.0,0.0),true);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    double expectedDisplacement = tipForce * barLength / (elasticityModulus * barArea);
    double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nodeId(nSlices,0,0));

    QVERIFY2(std::fabs(tipDisplacement-expectedDisplacement) < 0.02*expectedDisplacement,
             qPrintable(QString("Tip displacement = %1 m, expected %2 m").arg(tipDisplacement).arg(expectedDisplacement)));
}

//! A roller may only add a reaction along the direction it restrains.
//!
//! The nodal force is recovered from the global element stiffness, so at a free
//! degree of freedom it equals the load applied there - here the share of the
//! bar weight carried by that node, which does not depend on the constraints.
//! Adding the roller may therefore only change the nodal force along its own
//! local direction; the two components across it have to stay put. That fails
//! when the element load vector is rotated into the node local frame with the
//! wrong transform, because the solution then satisfies equilibrium with a
//! rotated load instead of the real one.
void TestSolverStress::reactionsBalanceGravityWithTiltedRoller()
{
    RModel withRoller = buildBarUnderGravityWithTiltedRoller();
    withRoller.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    RModel plain = buildBarUnderGravity();
    plain.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(withRoller,solverError),qPrintable(solverError));
    QVERIFY2(runSolver(plain,solverError),qPrintable(solverError));

    RR3Vector d(1.0,1.0,1.0);
    d.normalize();

    uint slice = nSlices/2;

    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            uint node = nodeId(slice,j,k);

            RR3Vector difference(0.0,0.0,0.0);
            double reference = 0.0;
            for (uint c=0;c<3;c++)
            {
                double a = nodeValue(withRoller,R_VARIABLE_FORCE,c,node);
                double b = nodeValue(plain,R_VARIABLE_FORCE,c,node);
                difference[c] = a - b;
                reference += b*b;
            }
            reference = std::sqrt(reference);

            double along = RR3Vector::dot(difference,d);
            double perpendicular = 0.0;
            for (uint c=0;c<3;c++)
            {
                double p = difference[c] - along*d[c];
                perpendicular += p*p;
            }
            perpendicular = std::sqrt(perpendicular);

            QVERIFY2(perpendicular < 0.01*reference,
                     qPrintable(QString("Node %1: the roller changed the nodal force by %2 N across the "
                                        "direction it restrains, against a nodal load of %3 N")
                                .arg(node).arg(perpendicular).arg(reference)));
        }
    }
}

//! Constraints coming from several entities have to combine on a shared node,
//! and a prescribed value has to survive the move into the node local frame.
//!
//! The whole volume carries a Displacement boundary condition holding Y and Z
//! at zero, and the free end face carries a roller along the bar axis with a
//! non-zero value. The tip node therefore collects three held directions from
//! two different entities, one of which prescribes a value. Reducing the node
//! to a single boundary condition would lose one of them; reading the global
//! Y and Z of the volume in the frame of the roller would hold the wrong
//! directions; and dropping the value would leave the bar unstretched.
void TestSolverStress::prescribedDisplacementIsAppliedInLocalFrame()
{
    const double prescribed = 1.0e-4; // [m]

    RModel model = buildBar(0.0,true);

    uint n00 = nodeId(nSlices,0,0);
    uint n01 = nodeId(nSlices,0,1);
    uint n10 = nodeId(nSlices,1,0);
    uint n11 = nodeId(nSlices,1,1);

    RElement t1(R_ELEMENT_TRI1);
    t1.setNodeId(0,n00);
    t1.setNodeId(1,n10);
    t1.setNodeId(2,n11);
    model.addElement(t1,true,0);

    RElement t2(R_ELEMENT_TRI1);
    t2.setNodeId(0,n00);
    t2.setNodeId(1,n11);
    t2.setNodeId(2,n01);
    model.addElement(t2,true,0);

    RSurface &surface = model.getSurface(0);
    surface.setName("Free end");
    surface.setThickness(0.0);

    RBoundaryCondition rollerBc(R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER);
    rollerBc.getComponent(rollerBc.findComponentPosition(R_VARIABLE_DISPLACEMENT)).add(0.0,prescribed);
    rollerBc.setExplicitLocalDirection(true);
    rollerBc.setLocalDirection(RR3Vector(1.0,0.0,0.0));
    surface.addBoundaryCondition(rollerBc);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));

    // The prescribed value arrives at the tip.
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            uint node = nodeId(nSlices,j,k);
            double tipDisplacement = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,node);
            QVERIFY2(std::fabs(tipDisplacement-prescribed) < 1.0e-3*prescribed,
                     qPrintable(QString("Node %1: displacement X = %2 m, expected the prescribed %3 m")
                                .arg(node).arg(tipDisplacement).arg(prescribed)));
        }
    }

    // The volume constraint is still read in global Y and Z.
    for (uint i=0;i<model.getNNodes();i++)
    {
        QVERIFY2(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,i)) < 1.0e-14,
                 qPrintable(QString("Node %1 moved in Y").arg(i)));
        QVERIFY2(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,2,i)) < 1.0e-14,
                 qPrintable(QString("Node %1 moved in Z").arg(i)));
    }

    // A uniform stretch of prescribed/L raises a uniform axial stress.
    double expectedStress = elasticityModulus * prescribed / barLength;
    uint nChecked = 0;
    for (uint i=0;i<model.getNElements();i++)
    {
        if (!R_ELEMENT_TYPE_IS_VOLUME(model.getElement(i).getType()))
        {
            continue;
        }
        RR3Vector center;
        model.getElement(i).findCenter(model.getNodes(),center[0],center[1],center[2]);
        if (center[0] < 0.2*barLength || center[0] > 0.8*barLength)
        {
            continue;
        }

        double stressX = elementValue(model,R_VARIABLE_STRESS_X,i);
        QVERIFY2(std::fabs(stressX-expectedStress) < 0.05*expectedStress,
                 qPrintable(QString("Element %1: stress X = %2 Pa, expected %3 Pa").arg(i).arg(stressX).arg(expectedStress)));
        nChecked++;
    }
    QVERIFY2(nChecked > 0,"No interior element was checked");
}

//! Two entities meeting at a node and prescribing different displacements in
//! the same direction cannot both be satisfied. That has to be reported rather
//! than silently resolved in favour of whichever was read last.
void TestSolverStress::contradictingConstraintsAreReported()
{
    RModel model = buildBar(0.0,false);

    // A second point entity on the nodes the first one already holds at zero,
    // this time asking them to move.
    for (uint j=0;j<2;j++)
    {
        for (uint k=0;k<2;k++)
        {
            RElement element(R_ELEMENT_POINT);
            element.setNodeId(0,nodeId(0,j,k));
            model.addElement(element,true,2);
        }
    }

    RPoint &movedPoint = model.getPoint(2);
    movedPoint.setName("Contradicting end");
    movedPoint.setVolume(0.0);

    RBoundaryCondition movedBc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    movedBc.getComponent(movedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,1.0e-4);
    movedBc.getComponent(movedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).setEnabled(false);
    movedBc.getComponent(movedBc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).setEnabled(false);
    movedPoint.addBoundaryCondition(movedBc);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));

    QString solverError;
    QVERIFY2(!runSolver(model,solverError),"The contradicting constraints were accepted");
    QVERIFY2(solverError.contains("contradict"),
             qPrintable(QString("Unexpected error message: %1").arg(solverError)));
}

void TestSolverStress::linearConstraintCouplesTrussNodes()
{
    constexpr double tipForce = 1000.0;
    RModel model = buildTruss(tipForce);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RMatrixSolverConf &gmres = model.getMatrixSolverConf(RMatrixSolverConf::GMRES);
    gmres.setNInnerIterations(50);
    gmres.setNOuterIterations(100);
    gmres.setSolverCvgValue(1.0e-10);
    gmres.setOutputFrequency(0);

    RLinearConstraint constraint;
    constraint.id = "tie-last-segment";
    constraint.name = "Equal tip displacement";
    constraint.terms = {{{nSlices,RDofType::UX},1.0},{{nSlices-1,RDofType::UX},-1.0}};
    model.getProblemSetup().getStructuralInteractionSetup().getLinearConstraints().push_back(constraint);

    QString solverError;
    QVERIFY2(runSolver(model,solverError),qPrintable(solverError));
    const double expected = tipForce*(barLength-barLength/double(nSlices))/(elasticityModulus*barArea);
    const double tip = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nSlices);
    const double previous = nodeValue(model,R_VARIABLE_DISPLACEMENT,0,nSlices-1);
    QVERIFY2(std::fabs(tip-previous) < 1.0e-8*expected,
             qPrintable(QString("MPC residual = %1").arg(tip-previous)));
    QVERIFY2(std::fabs(tip-expected) < 2.0e-5*expected,
             qPrintable(QString("MPC tip = %1, expected %2").arg(tip).arg(expected)));
}

void TestSolverStress::contradictoryLinearConstraintsAreReported()
{
    RModel model = buildTruss(0.0);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RLinearConstraint first;
    first.id = "tip-zero";
    first.terms = {{{nSlices,RDofType::UX},1.0}};
    RLinearConstraint second = first;
    second.id = "tip-one";
    second.rhs = 1.0;
    RStructuralInteractionSetup &setup = model.getProblemSetup().getStructuralInteractionSetup();
    setup.getLinearConstraints().push_back(first);
    setup.getLinearConstraints().push_back(second);

    QString solverError;
    QVERIFY2(!runSolver(model,solverError),"Contradictory MPC equations were accepted");
    QVERIFY2(solverError.contains("contradict",Qt::CaseInsensitive),qPrintable(solverError));
}

QTEST_APPLESS_MAIN(TestSolverStress)

#include "tst_solver_stress.moc"
