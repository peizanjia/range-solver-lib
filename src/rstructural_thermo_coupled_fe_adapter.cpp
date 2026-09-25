#include "rstructural_thermo_coupled_fe_adapter.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <rml_element.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_material_property.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_environment_condition.h>
#include <rml_condition_component.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

namespace
{
constexpr double geometryTolerance = 1.0e-12;

[[noreturn]] void invalidInput(const char *message)
{
    throw std::invalid_argument(std::string("Thermo-mechanical RModel adapter: ")+message);
}

void require(bool condition,const char *message)
{
    if (!condition)
        invalidInput(message);
}

bool isSingleHeatThenStressTask(const RProblemTaskItem &root)
{
    return root.getProblemType() == R_PROBLEM_NONE && root.getNIterations() == 1U &&
           root.getNChildren() == 2U &&
           root.getChild(0).getProblemType() == R_PROBLEM_HEAT &&
           root.getChild(0).getNChildren() == 0U &&
           root.getChild(1).getProblemType() == R_PROBLEM_STRESS &&
           root.getChild(1).getNChildren() == 0U;
}

bool sampleEquals(const RConditionComponent &component,double expected)
{
    if (!component.getEnabled())
        return true;
    for (uint i=0;i<component.size();i++)
    {
        const double key = component.getKey(i);
        const double value = component.getValue(i);
        if (!std::isfinite(key) || !std::isfinite(value) ||
            std::abs(value-expected) > 1.0e-12*std::max({1.0,std::abs(expected),std::abs(value)}))
            return false;
    }
    return true;
}

bool neutralEnvironment(const REnvironmentCondition &condition,double referenceTemperatureKelvin)
{
    if (condition.getType() == R_ENVIRONMENT_CONDITION_G_ACCELERATION)
    {
        for (uint i=0;i<condition.size();i++)
        {
            const RConditionComponent &component = condition.getComponent(i);
            const RVariableType type = component.getType();
            if (type != R_VARIABLE_G_ACCELERATION_X && type != R_VARIABLE_G_ACCELERATION_Y &&
                type != R_VARIABLE_G_ACCELERATION_Z)
                return false;
            if (!sampleEquals(component,0.0))
                return false;
        }
        return true;
    }
    if (condition.getType() == R_ENVIRONMENT_CONDITION_TEMPERATURE)
    {
        for (uint i=0;i<condition.size();i++)
        {
            const RConditionComponent &component = condition.getComponent(i);
            if (component.getType() != R_VARIABLE_TEMPERATURE ||
                !sampleEquals(component,referenceTemperatureKelvin))
                return false;
        }
        return true;
    }
    for (uint i=0;i<condition.size();i++)
    {
        const RConditionComponent &component = condition.getComponent(i);
        if (component.getEnabled() && !component.empty())
            return false;
    }
    return true;
}

void rejectUnsupportedModelConditions(const RModel &model,double referenceTemperatureKelvin)
{
    for (uint groupIndex=0;groupIndex<model.getNElementGroups();groupIndex++)
    {
        const RElementGroup *group = model.getElementGroupPtr(groupIndex);
        require(group != nullptr,"RModel contains a null element group.");
        require(group->getNBoundaryConditions() == 0U,
                "Explicit kernel conditions are required; model boundary conditions are unsupported.");
        require(group->getNInitialConditions() == 0U,
                "Accepted nodal fields are the only initial state; model initial conditions are unsupported.");
        for (uint conditionIndex=0;conditionIndex<group->getNEnvironmentConditions();conditionIndex++)
            require(neutralEnvironment(group->getEnvironmentCondition(conditionIndex),
                                      referenceTemperatureKelvin),
                    "Model environment conditions must be neutral for the explicit 1-D kernel inputs.");
    }
}

RMaterialProperty propertyFor(const RMaterial &material,RMaterialProperty::Type type)
{
    const uint position = material.findPosition(type);
    require(position < material.size(),"TRUSS1 material is missing a required property.");
    const RMaterialProperty property = material.get(position);
    require(!property.empty(),"TRUSS1 material property has no temperature samples.");
    for (uint i=0;i<property.size();i++)
        require(std::isfinite(property.getKey(i)) && std::isfinite(property.getValue(i)),
                "TRUSS1 material property contains a non-finite sample.");
    return property;
}

double constantProperty(const RMaterial &material,RMaterialProperty::Type type)
{
    const RMaterialProperty property = propertyFor(material,type);
    const double value = property.getValue(0);
    for (uint i=1;i<property.size();i++)
        require(property.getValue(i) == value,
                "Thermo-coupled 1-D kernel requires constant density, heat capacity, and conductivity.");
    require(std::isfinite(value) && value > 0.0,
            "TRUSS1 thermal material property must be finite and positive.");
    return value;
}

struct ElementMaterialBinding
{
    unsigned int lineIndex = 0;
    double area = 0.0;
    double density = 0.0;
    double specificHeat = 0.0;
    double conductivity = 0.0;
    RMaterialProperty youngModulus;
    RMaterialProperty expansion;
};

std::vector<ElementMaterialBinding> elementMaterialBindings(const RModel &model)
{
    require(model.getNElements() > 0,"RModel contains no TRUSS1 elements.");
    require(model.getNLines() > 0,"TRUSS1 elements must belong to RLine groups.");

    std::vector<ElementMaterialBinding> bindings(model.getNElements());
    std::vector<unsigned int> membershipCount(model.getNElements(),0U);
    for (uint lineIndex=0;lineIndex<model.getNLines();lineIndex++)
    {
        const RLine &line = model.getLine(lineIndex);
        if (line.empty())
            continue;
        const double area = line.getCrossArea();
        require(std::isfinite(area) && area > 0.0,
                "Every populated RLine group requires a finite positive cross-sectional area.");
        const RMaterial &material = line.getMaterial();
        require(material.getState() == RMaterial::Solid,
                "Thermo-mechanical TRUSS1 adapter requires solid line materials.");

        const RMaterialProperty young = propertyFor(material,RMaterialProperty::ModulusOfElasiticity);
        const RMaterialProperty expansion = propertyFor(
            material,RMaterialProperty::ThermalExpansionCoefficient);
        const double density = constantProperty(material,RMaterialProperty::Density);
        const double specificHeat = constantProperty(material,RMaterialProperty::HeatCapacity);
        const double conductivity = constantProperty(material,RMaterialProperty::ThermalConductivity);

        for (uint position=0;position<line.size();position++)
        {
            const unsigned int elementIndex = line.get(position);
            require(elementIndex < model.getNElements(),"RLine group references an invalid element.");
            require(model.getElement(elementIndex).getType() == R_ELEMENT_TRUSS1,
                    "Only TRUSS1 elements are supported; mixed line formulations fail closed.");
            ++membershipCount[elementIndex];
            require(membershipCount[elementIndex] == 1U,
                    "Each TRUSS1 element must belong to exactly one RLine material group.");
            bindings[elementIndex] = {lineIndex,area,density,specificHeat,conductivity,young,expansion};
        }
    }

    for (uint elementIndex=0;elementIndex<model.getNElements();elementIndex++)
    {
        require(model.getElement(elementIndex).getType() == R_ELEMENT_TRUSS1,
                "The bounded thermo-mechanical adapter accepts TRUSS1 elements only.");
        require(membershipCount[elementIndex] == 1U,
                "Each TRUSS1 element must belong to exactly one RLine material group.");
    }
    return bindings;
}

void validateVariable(const RModel &model,RVariableType type,RVariableApplyType applyType,
                      uint nVectors,uint nValues,const char *message)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
        return;
    const RVariable &variable = model.getVariable(position);
    require(variable.getApplyType() == applyType && variable.getNVectors() == nVectors &&
                variable.getNValues() == nValues,message);
}

RVariable makeVariable(RVariableType type,RVariableApplyType applyType,uint vectors,uint values)
{
    RVariable variable(type,applyType);
    variable.resize(vectors,values,true);
    return variable;
}
}

RStructuralThermoCoupled1DResult RStructuralThermoCoupledFEAdapter::advance(
    RModel &model,
    const RStructuralThermoFEViscosityByElement &viscosityPaSecondsByElement,
    const RStructuralThermoCoupled1DConditions &conditions,
    const RStructuralThermoCoupled1DOptions &options,
    double referenceTemperatureKelvin)
{
    require(model.getNNodes() >= 2,"RModel requires at least two nodes.");
    require(isSingleHeatThenStressTask(model.getProblemTaskTree()),
            "RModel task tree must be exactly one Heat -> Stress iteration with two leaf tasks.");
    require(!model.getProblemSetup().getStructuralDynamicsSetup().getEnabled(),
            "The 1-D thermo-mechanical kernel is quasi-static and cannot advance a dynamics-enabled RModel.");
    require(std::isfinite(referenceTemperatureKelvin) && referenceTemperatureKelvin >= 0.0,
            "Reference temperature must be finite and nonnegative kelvin.");
    require(!model.getTimeSolver().getEnabled(),
            "This direct macrostep adapter requires the RModel time scheduler to be disabled.");
    const double acceptedTime = model.getTimeSolver().getComputedTime();
    require(std::isfinite(acceptedTime) && std::isfinite(options.dt) && options.dt > 0.0 &&
                std::isfinite(acceptedTime+options.dt),
            "Accepted RModel time plus the physical step interval must remain finite.");
    rejectUnsupportedModelConditions(model,referenceTemperatureKelvin);

    const auto bindings = elementMaterialBindings(model);
    require(viscosityPaSecondsByElement.size() == model.getNElements(),
            "Supply one explicit physical Kelvin-Voigt viscosity for every TRUSS1 element index.");
    for (auto iter=viscosityPaSecondsByElement.cbegin();iter!=viscosityPaSecondsByElement.cend();++iter)
        require(iter.key() < model.getNElements() && std::isfinite(iter.value()) && iter.value() >= 0.0,
                "Explicit TRUSS1 Kelvin-Voigt viscosities must be nonnegative finite Pa*s values keyed by valid element index.");

    const double referenceY = model.getNode(0).getY();
    const double referenceZ = model.getNode(0).getZ();
    require(std::isfinite(referenceY) && std::isfinite(referenceZ),
            "RModel coordinates must be finite.");
    for (uint node=0;node<model.getNNodes();node++)
    {
        const RNode &position = model.getNode(node);
        require(std::isfinite(position.getX()) && std::isfinite(position.getY()) &&
                    std::isfinite(position.getZ()),"RModel coordinates must be finite.");
        require(std::abs(position.getY()-referenceY) <= geometryTolerance &&
                    std::abs(position.getZ()-referenceZ) <= geometryTolerance,
                "Only x-aligned 1-D TRUSS1 geometry is supported.");
    }

    const uint temperaturePosition = model.findVariable(R_VARIABLE_TEMPERATURE);
    const uint displacementPosition = model.findVariable(R_VARIABLE_DISPLACEMENT);
    validateVariable(model,R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE,1U,model.getNNodes(),
                     "Existing temperature must be one nodal kelvin component.");
    validateVariable(model,R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3U,model.getNNodes(),
                     "Existing displacement must have three nodal components.");

    RStructuralThermoCoupled1DModel kernelModel;
    kernelModel.referenceTemperatureKelvin = referenceTemperatureKelvin;
    kernelModel.nodeIds.reserve(model.getNNodes());
    std::vector<RStructuralThermoCoupled1DNodeState> accepted;
    accepted.reserve(model.getNNodes());
    for (uint node=0;node<model.getNNodes();node++)
    {
        const QString id = QString::number(node);
        kernelModel.nodeIds.push_back(id);
        const double temperature = temperaturePosition == RConstants::eod
            ? referenceTemperatureKelvin
            : model.getVariable(temperaturePosition).getValue(0,node);
        const double axialDisplacement = displacementPosition == RConstants::eod
            ? 0.0 : model.getVariable(displacementPosition).getValue(0,node);
        if (displacementPosition != RConstants::eod)
        {
            const RVariable &displacement = model.getVariable(displacementPosition);
            require(std::abs(displacement.getValue(1,node)) <= geometryTolerance &&
                        std::abs(displacement.getValue(2,node)) <= geometryTolerance,
                    "Existing transverse displacement is outside the x-aligned adapter scope.");
        }
        require(std::isfinite(temperature) && temperature >= 0.0 &&
                    std::isfinite(axialDisplacement),
                "Initial nodal temperature or axial displacement is invalid.");
        accepted.push_back({id,temperature,axialDisplacement});
    }

    kernelModel.members.reserve(model.getNElements());
    std::vector<double> elementAreas(model.getNElements(),0.0);
    std::vector<bool> usedNodes(model.getNNodes(),false);
    for (uint elementIndex=0;elementIndex<model.getNElements();elementIndex++)
    {
        const RElement &element = model.getElement(elementIndex);
        require(element.size() == 2 && !element.hasDuplicateNodes(),
                "TRUSS1 element must contain two distinct node indices.");
        const unsigned int firstNode = element.getNodeId(0);
        const unsigned int secondNode = element.getNodeId(1);
        require(firstNode < model.getNNodes() && secondNode < model.getNNodes(),
                "TRUSS1 element references an invalid node.");
        const double firstX = model.getNode(firstNode).getX();
        const double secondX = model.getNode(secondNode).getX();
        const double length = std::abs(secondX-firstX);
        require(std::isfinite(length) && length > 0.0,
                "TRUSS1 member must have positive length along global X.");

        const unsigned int node0 = firstX < secondX ? firstNode : secondNode;
        const unsigned int node1 = firstX < secondX ? secondNode : firstNode;
        const ElementMaterialBinding &binding = bindings[elementIndex];
        RStructuralThermoCoupled1DMember member;
        member.elementId = QStringLiteral("element-")+QString::number(elementIndex);
        member.node0Id = QString::number(node0);
        member.node1Id = QString::number(node1);
        member.area = binding.area;
        member.length = length;
        member.conductivity = binding.conductivity;
        member.density = binding.density;
        member.specificHeat = binding.specificHeat;
        member.internalViscosity = viscosityPaSecondsByElement.value(elementIndex);
        member.youngModulusAtKelvin = [property=binding.youngModulus](double kelvin) {
            return property.get(kelvin);
        };
        member.expansionAtKelvin = [property=binding.expansion](double kelvin) {
            return property.get(kelvin);
        };
        kernelModel.members.push_back(std::move(member));
        elementAreas[elementIndex] = binding.area;
        usedNodes[firstNode] = true;
        usedNodes[secondNode] = true;
    }
    require(std::all_of(usedNodes.cbegin(),usedNodes.cend(),[](bool used) { return used; }),
            "Unused isolated nodes are outside the bounded thermo-mechanical mesh.");

    // Keep a full RModel candidate throughout the solve. The core and every
    // output write must succeed before replacing the caller's accepted model.
    RModel candidate(model);
    const RStructuralThermoCoupled1DResult result = RStructuralThermoCoupled1D::step(
        kernelModel,accepted,conditions,options);
    require(result.nextState.size() == model.getNNodes() &&
                result.memberAxialForcesNewtons.size() == model.getNElements(),
            "Certified thermo-mechanical result dimensions differ from the RModel mesh.");

    RVariable temperature = makeVariable(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE,1U,
                                         candidate.getNNodes());
    RVariable displacement = makeVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3U,
                                          candidate.getNNodes());
    for (uint node=0;node<candidate.getNNodes();node++)
    {
        require(result.nextState[node].nodeId == QString::number(node),
                "Certified thermo-mechanical node order differs from RModel.");
        temperature.setValue(0,node,result.nextState[node].temperatureKelvin);
        displacement.setValue(0,node,result.nextState[node].axialDisplacement);
        displacement.setValue(1,node,0.0);
        displacement.setValue(2,node,0.0);
    }
    candidate.addVariable(temperature);
    candidate.addVariable(displacement);

    RVariable stressX = makeVariable(R_VARIABLE_STRESS_X,R_VARIABLE_APPLY_ELEMENT,1U,
                                     candidate.getNElements());
    RVariable stressY = makeVariable(R_VARIABLE_STRESS_Y,R_VARIABLE_APPLY_ELEMENT,1U,
                                     candidate.getNElements());
    RVariable stressZ = makeVariable(R_VARIABLE_STRESS_Z,R_VARIABLE_APPLY_ELEMENT,1U,
                                     candidate.getNElements());
    RVariable stressXY = makeVariable(R_VARIABLE_STRESS_XY,R_VARIABLE_APPLY_ELEMENT,1U,
                                      candidate.getNElements());
    RVariable stressXZ = makeVariable(R_VARIABLE_STRESS_XZ,R_VARIABLE_APPLY_ELEMENT,1U,
                                      candidate.getNElements());
    RVariable stressYZ = makeVariable(R_VARIABLE_STRESS_YZ,R_VARIABLE_APPLY_ELEMENT,1U,
                                      candidate.getNElements());
    for (uint element=0;element<candidate.getNElements();element++)
    {
        const double sigmaX = result.memberAxialForcesNewtons[element]/elementAreas[element];
        require(std::isfinite(sigmaX),"Recovered axial stress is non-finite.");
        stressX.setValue(0,element,sigmaX);
    }
    candidate.addVariable(stressX);
    candidate.addVariable(stressY);
    candidate.addVariable(stressZ);
    candidate.addVariable(stressXY);
    candidate.addVariable(stressXZ);
    candidate.addVariable(stressYZ);

    RVariable vonMises = makeVariable(R_VARIABLE_STRESS_VON_MISES,R_VARIABLE_APPLY_ELEMENT,1U,
                                      candidate.getNElements());
    RVariable normal = makeVariable(R_VARIABLE_STRESS_NORMAL,R_VARIABLE_APPLY_ELEMENT,1U,
                                   candidate.getNElements());
    RVariable shear = makeVariable(R_VARIABLE_STRESS_SHEAR,R_VARIABLE_APPLY_ELEMENT,1U,
                                   candidate.getNElements());
    for (uint element=0;element<candidate.getNElements();element++)
    {
        const double sigmaX = stressX.getValue(0,element);
        vonMises.setValue(0,element,std::abs(sigmaX));
        normal.setValue(0,element,sigmaX);
        shear.setValue(0,element,0.0);
    }
    candidate.addVariable(vonMises);
    candidate.addVariable(normal);
    candidate.addVariable(shear);

    candidate.getTimeSolver().setComputedTime(acceptedTime+options.dt);

    model = candidate;
    return result;
}
