#include "rstructural_thermo_monolithic_1d_fe_adapter.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <rml_condition_component.h>
#include <rml_element.h>
#include <rml_environment_condition.h>
#include <rml_element_group.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_material_property.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

namespace
{
constexpr double geometryTolerance = 1.0e-12;
constexpr unsigned int maximumNodeCount = 64U;
constexpr unsigned int maximumMemberCount = 1024U;

[[noreturn]] void invalidInput(const char *message)
{
    throw std::invalid_argument(std::string("Monolithic thermoelastic RModel adapter: ")+message);
}

void require(bool condition,const char *message)
{
    if (!condition)
    {
        invalidInput(message);
    }
}

void requireFinite(double value,const char *message)
{
    require(std::isfinite(value),message);
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
    {
        return true;
    }
    for (uint i=0;i<component.size();i++)
    {
        const double key = component.getKey(i);
        const double value = component.getValue(i);
        if (!std::isfinite(key) || !std::isfinite(value) ||
            std::abs(value-expected) > 1.0e-12*
                std::max({1.0,std::abs(expected),std::abs(value)}))
        {
            return false;
        }
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
            {
                return false;
            }
            if (!sampleEquals(component,0.0))
            {
                return false;
            }
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
            {
                return false;
            }
        }
        return true;
    }
    for (uint i=0;i<condition.size();i++)
    {
        const RConditionComponent &component = condition.getComponent(i);
        if (component.getEnabled() && !component.empty())
        {
            return false;
        }
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
                "Supply boundary conditions through the explicit monolithic conditions arrays.");
        require(group->getNInitialConditions() == 0U,
                "Accepted nodal fields are the only supported initial state.");
        for (uint conditionIndex=0;conditionIndex<group->getNEnvironmentConditions();conditionIndex++)
        {
            require(neutralEnvironment(group->getEnvironmentCondition(conditionIndex),
                                       referenceTemperatureKelvin),
                    "RModel environment conditions must be neutral for explicit monolithic inputs.");
        }
    }
}

RMaterialProperty propertyFor(const RMaterial &material,RMaterialProperty::Type type)
{
    const uint position = material.findPosition(type);
    require(position < material.size(),"TRUSS1 material is missing a required property.");
    const RMaterialProperty property = material.get(position);
    require(!property.empty(),"TRUSS1 material property has no temperature samples.");
    for (uint i=0;i<property.size();i++)
    {
        require(std::isfinite(property.getKey(i)) && std::isfinite(property.getValue(i)),
                "TRUSS1 material property contains a non-finite sample.");
    }
    return property;
}

double constantProperty(const RMaterial &material,RMaterialProperty::Type type)
{
    const RMaterialProperty property = propertyFor(material,type);
    const double value = property.getValue(0);
    for (uint i=1;i<property.size();i++)
    {
        require(property.getValue(i) == value,
                "The monolithic 1-D kernel requires temperature-independent material properties.");
    }
    return value;
}

struct ElementMaterialBinding
{
    double area = 0.0;
    double youngModulus = 0.0;
    double expansion = 0.0;
    double density = 0.0;
    double specificHeat = 0.0;
    double conductivity = 0.0;
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
        {
            continue;
        }
        const double area = line.getCrossArea();
        require(std::isfinite(area) && area > 0.0,
                "Every populated RLine requires a finite positive cross-sectional area.");
        const RMaterial &material = line.getMaterial();
        require(material.getState() == RMaterial::Solid,
                "TRUSS1 members must use solid materials.");

        const double young = constantProperty(material,RMaterialProperty::ModulusOfElasiticity);
        const double expansion = constantProperty(
            material,RMaterialProperty::ThermalExpansionCoefficient);
        const double density = constantProperty(material,RMaterialProperty::Density);
        const double specificHeat = constantProperty(material,RMaterialProperty::HeatCapacity);
        const double conductivity = constantProperty(material,RMaterialProperty::ThermalConductivity);
        require(std::isfinite(young) && young > 0.0,
                "TRUSS1 Young modulus must be finite and positive.");
        require(std::isfinite(expansion),"TRUSS1 thermal expansion must be finite.");
        require(std::isfinite(density) && density > 0.0,
                "TRUSS1 density must be finite and positive.");
        require(std::isfinite(specificHeat) && specificHeat > 0.0,
                "TRUSS1 specific heat must be finite and positive.");
        require(std::isfinite(conductivity) && conductivity >= 0.0,
                "TRUSS1 conductivity must be finite and nonnegative.");

        for (uint position=0;position<line.size();position++)
        {
            const unsigned int element = line.get(position);
            require(element < model.getNElements(),"RLine references an invalid element.");
            require(model.getElement(element).getType() == R_ELEMENT_TRUSS1,
                    "Every RLine member must be TRUSS1.");
            ++membershipCount[element];
            require(membershipCount[element] == 1U,
                    "Each TRUSS1 element must belong to exactly one RLine material group.");
            bindings[element] = {area,young,expansion,density,specificHeat,conductivity};
        }
    }

    for (uint element=0;element<model.getNElements();element++)
    {
        require(model.getElement(element).getType() == R_ELEMENT_TRUSS1,
                "The bounded monolithic adapter accepts TRUSS1 elements only.");
        require(membershipCount[element] == 1U,
                "Each TRUSS1 element must belong to exactly one RLine material group.");
    }
    return bindings;
}

void validateVariable(const RModel &model,RVariableType type,RVariableApplyType applyType,
                      uint nVectors,uint nValues,const char *message)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return;
    }
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

void validateConditions(const RStructuralThermoMonolithic1DConditions &conditions,
                        uint nNodes)
{
    require(conditions.mechanicalForceNewtons.size() == nNodes &&
                conditions.displacementPrescribed.size() == nNodes &&
                conditions.prescribedDisplacementMetres.size() == nNodes &&
                conditions.heatPowerWatts.size() == nNodes &&
                conditions.temperaturePrescribed.size() == nNodes &&
                conditions.prescribedTemperatureKelvin.size() == nNodes,
            "Every explicit load, mask and prescribed-value array must match the RModel node count.");
    for (uint node=0;node<nNodes;node++)
    {
        requireFinite(conditions.mechanicalForceNewtons[node],"Mechanical forces must be finite.");
        requireFinite(conditions.prescribedDisplacementMetres[node],
                      "Prescribed displacements must be finite.");
        requireFinite(conditions.heatPowerWatts[node],"Heat powers must be finite.");
        requireFinite(conditions.prescribedTemperatureKelvin[node],
                      "Prescribed temperatures must be finite.");
    }
}

void validateResult(const RStructuralThermoMonolithic1DResult &result,
                    uint nNodes,uint nElements)
{
    require(result.nextState.temperatureKelvin.size() == nNodes &&
                result.nextState.axialDisplacementMetres.size() == nNodes,
            "Monolithic candidate state dimensions differ from the RModel node count.");
    require(result.memberAxialStrain.size() == nElements &&
                result.memberAxialStressPascals.size() == nElements &&
                result.memberAxialForceNewtons.size() == nElements &&
                result.memberDissipatedEnergyJoules.size() == nElements,
            "Monolithic candidate element diagnostics differ from the RModel element count.");
    for (uint node=0;node<nNodes;node++)
    {
        require(std::isfinite(result.nextState.temperatureKelvin[node]) &&
                    result.nextState.temperatureKelvin[node] > 0.0 &&
                    std::isfinite(result.nextState.axialDisplacementMetres[node]),
                "Monolithic candidate nodal state is invalid.");
    }
    for (uint element=0;element<nElements;element++)
    {
        require(std::isfinite(result.memberAxialStrain[element]) &&
                    std::isfinite(result.memberAxialStressPascals[element]) &&
                    std::isfinite(result.memberAxialForceNewtons[element]) &&
                    std::isfinite(result.memberDissipatedEnergyJoules[element]) &&
                    result.memberDissipatedEnergyJoules[element] >= 0.0,
                "Monolithic candidate element result is invalid.");
    }
    const double diagnostics[] = {
        result.mechanicalResidualNormNewtons,
        result.thermalResidualNormWatts,
        result.thermalEnergyChangeJoules,
        result.heaterEnergyJoules,
        result.reversibleThermoelasticEnergyJoules,
        result.viscousDissipationEnergyJoules,
        result.dissipationHeatEnergyJoules,
        result.thermalBoundaryEnergyJoules,
        result.thermalEnergyBalanceResidualJoules
    };
    for (double value : diagnostics)
    {
        require(std::isfinite(value),"Monolithic candidate energy or residual diagnostics are non-finite.");
    }
    require(result.mechanicalResidualNormNewtons >= 0.0 &&
                result.thermalResidualNormWatts >= 0.0 &&
                result.viscousDissipationEnergyJoules >= 0.0 &&
                result.dissipationHeatEnergyJoules >= 0.0,
            "Monolithic residual norms and dissipated energies must be nonnegative.");
}
}

RStructuralThermoMonolithic1DResult RStructuralThermoMonolithic1DFEAdapter::advance(
    RModel &model,
    const RStructuralThermoMonolithicFEViscosityByElement &viscosityPaSecondsByElement,
    const RStructuralThermoMonolithic1DConditions &conditions,
    const RStructuralThermoMonolithic1DOptions &options,
    double referenceTemperatureKelvin)
{
    require(model.getNNodes() >= 2U && model.getNNodes() <= maximumNodeCount,
            "RModel must contain 2 to 64 nodes.");
    require(model.getNElements() > 0U && model.getNElements() <= maximumMemberCount,
            "RModel must contain 1 to 1024 TRUSS1 elements.");
    require(isSingleHeatThenStressTask(model.getProblemTaskTree()),
            "RModel task tree must declare one Heat -> Stress iteration with two leaf tasks.");
    require(!model.getProblemSetup().getStructuralDynamicsSetup().getEnabled(),
            "The quasi-static monolithic adapter cannot advance a dynamics-enabled RModel.");
    require(std::isfinite(referenceTemperatureKelvin) && referenceTemperatureKelvin > 0.0,
            "Reference temperature must be finite positive absolute Kelvin.");
    require(!model.getTimeSolver().getEnabled(),
            "The direct monolithic adapter requires the RModel time scheduler to be disabled.");
    const double acceptedTime = model.getTimeSolver().getComputedTime();
    require(std::isfinite(acceptedTime) && std::isfinite(options.timeStepSeconds) &&
                options.timeStepSeconds > 0.0 &&
                std::isfinite(acceptedTime+options.timeStepSeconds),
            "Accepted time plus the monolithic time step must remain finite.");
    validateConditions(conditions,model.getNNodes());
    rejectUnsupportedModelConditions(model,referenceTemperatureKelvin);

    const std::vector<ElementMaterialBinding> bindings = elementMaterialBindings(model);
    require(viscosityPaSecondsByElement.size() == model.getNElements(),
            "Supply one explicit Pa*s Kelvin-Voigt coefficient for every current element index.");
    for (auto iter=viscosityPaSecondsByElement.cbegin();
         iter!=viscosityPaSecondsByElement.cend();++iter)
    {
        require(iter.key() < model.getNElements() && std::isfinite(iter.value()) &&
                    iter.value() >= 0.0,
                "Element-index Kelvin-Voigt coefficients must be finite nonnegative Pa*s.");
    }

    const RNode &firstPosition = model.getNode(0);
    const double referenceY = firstPosition.getY();
    const double referenceZ = firstPosition.getZ();
    require(std::isfinite(referenceY) && std::isfinite(referenceZ),
            "RModel reference coordinates must be finite.");
    RStructuralThermoMonolithic1DModel kernelModel;
    kernelModel.referenceTemperatureKelvin = referenceTemperatureKelvin;
    kernelModel.nodeCoordinatesMetres.reserve(model.getNNodes());
    std::vector<bool> usedNodes(model.getNNodes(),false);
    for (uint node=0;node<model.getNNodes();node++)
    {
        const RNode &position = model.getNode(node);
        require(std::isfinite(position.getX()) && std::isfinite(position.getY()) &&
                    std::isfinite(position.getZ()),"RModel coordinates must be finite.");
        require(std::abs(position.getY()-referenceY) <= geometryTolerance &&
                    std::abs(position.getZ()-referenceZ) <= geometryTolerance,
                "Only straight global-X TRUSS1 geometry is supported.");
        kernelModel.nodeCoordinatesMetres.push_back(position.getX());
    }

    const uint temperaturePosition = model.findVariable(R_VARIABLE_TEMPERATURE);
    const uint displacementPosition = model.findVariable(R_VARIABLE_DISPLACEMENT);
    validateVariable(model,R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE,1U,model.getNNodes(),
                     "Existing temperature must be one nodal kelvin component.");
    validateVariable(model,R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3U,model.getNNodes(),
                     "Existing displacement must have three nodal components.");

    RStructuralThermoMonolithic1DState accepted;
    accepted.temperatureKelvin.resize(model.getNNodes(),referenceTemperatureKelvin);
    accepted.axialDisplacementMetres.resize(model.getNNodes(),0.0);
    for (uint node=0;node<model.getNNodes();node++)
    {
        if (temperaturePosition != RConstants::eod)
        {
            accepted.temperatureKelvin[node] = model.getVariable(temperaturePosition).getValue(0,node);
        }
        if (displacementPosition != RConstants::eod)
        {
            const RVariable &displacement = model.getVariable(displacementPosition);
            const double transverseY = displacement.getValue(1,node);
            const double transverseZ = displacement.getValue(2,node);
            require(std::isfinite(transverseY) && std::isfinite(transverseZ) &&
                        std::abs(transverseY) <= geometryTolerance &&
                        std::abs(transverseZ) <= geometryTolerance,
                    "Existing transverse displacement is outside the 1-D adapter scope.");
            accepted.axialDisplacementMetres[node] = displacement.getValue(0,node);
        }
    }

    kernelModel.members.reserve(model.getNElements());
    for (uint elementIndex=0;elementIndex<model.getNElements();elementIndex++)
    {
        const RElement &element = model.getElement(elementIndex);
        require(element.size() == 2U && !element.hasDuplicateNodes(),
                "TRUSS1 elements must connect two distinct nodes.");
        const unsigned int firstNode = element.getNodeId(0);
        const unsigned int secondNode = element.getNodeId(1);
        require(firstNode < model.getNNodes() && secondNode < model.getNNodes(),
                "TRUSS1 element references an invalid node.");
        const double firstX = model.getNode(firstNode).getX();
        const double secondX = model.getNode(secondNode).getX();
        const unsigned int node0 = firstX < secondX ? firstNode : secondNode;
        const unsigned int node1 = firstX < secondX ? secondNode : firstNode;
        const ElementMaterialBinding &binding = bindings[elementIndex];
        RStructuralThermoMonolithic1DMember member;
        member.node0 = node0;
        member.node1 = node1;
        member.areaSquareMetres = binding.area;
        member.youngModulusPascals = binding.youngModulus;
        member.expansionPerKelvin = binding.expansion;
        member.densityKilogramsPerCubicMetre = binding.density;
        member.specificHeatJoulesPerKilogramKelvin = binding.specificHeat;
        member.conductivityWattsPerMetreKelvin = binding.conductivity;
        member.kelvinVoigtViscosityPascalSeconds = viscosityPaSecondsByElement.value(elementIndex);
        kernelModel.members.push_back(member);
        usedNodes[firstNode] = true;
        usedNodes[secondNode] = true;
    }
    require(std::all_of(usedNodes.cbegin(),usedNodes.cend(),[](bool used) { return used; }),
            "Unused isolated RModel nodes are outside the bounded monolithic adapter.");

    const RStructuralThermoMonolithic1DResult result = RStructuralThermoMonolithic1D::step(
        kernelModel,accepted,conditions,options);
    validateResult(result,model.getNNodes(),model.getNElements());

    // Build every result on a private model copy; no input model field changes
    // before all candidate variables and time metadata have passed validation.
    RModel candidate(model);
    RVariable temperature = makeVariable(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE,1U,
                                         candidate.getNNodes());
    RVariable displacement = makeVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3U,
                                          candidate.getNNodes());
    for (uint node=0;node<candidate.getNNodes();node++)
    {
        temperature.setValue(0,node,result.nextState.temperatureKelvin[node]);
        displacement.setValue(0,node,result.nextState.axialDisplacementMetres[node]);
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
    RVariable vonMises = makeVariable(R_VARIABLE_STRESS_VON_MISES,R_VARIABLE_APPLY_ELEMENT,1U,
                                      candidate.getNElements());
    RVariable normal = makeVariable(R_VARIABLE_STRESS_NORMAL,R_VARIABLE_APPLY_ELEMENT,1U,
                                    candidate.getNElements());
    RVariable shear = makeVariable(R_VARIABLE_STRESS_SHEAR,R_VARIABLE_APPLY_ELEMENT,1U,
                                   candidate.getNElements());
    for (uint element=0;element<candidate.getNElements();element++)
    {
        const double sigmaX = result.memberAxialStressPascals[element];
        stressX.setValue(0,element,sigmaX);
        stressY.setValue(0,element,0.0);
        stressZ.setValue(0,element,0.0);
        stressXY.setValue(0,element,0.0);
        stressXZ.setValue(0,element,0.0);
        stressYZ.setValue(0,element,0.0);
        vonMises.setValue(0,element,std::abs(sigmaX));
        normal.setValue(0,element,sigmaX);
        shear.setValue(0,element,0.0);
    }
    candidate.addVariable(stressX);
    candidate.addVariable(stressY);
    candidate.addVariable(stressZ);
    candidate.addVariable(stressXY);
    candidate.addVariable(stressXZ);
    candidate.addVariable(stressYZ);
    candidate.addVariable(vonMises);
    candidate.addVariable(normal);
    candidate.addVariable(shear);

    candidate.getTimeSolver().setComputedTime(acceptedTime+options.timeStepSeconds);
    require(std::isfinite(candidate.getTimeSolver().getComputedTime()),
            "Candidate computed time is non-finite.");
    const uint candidateTemperaturePosition = candidate.findVariable(R_VARIABLE_TEMPERATURE);
    const uint candidateDisplacementPosition = candidate.findVariable(R_VARIABLE_DISPLACEMENT);
    const uint candidateStressPosition = candidate.findVariable(R_VARIABLE_STRESS_X);
    require(candidateTemperaturePosition != RConstants::eod &&
                candidateDisplacementPosition != RConstants::eod &&
                candidateStressPosition != RConstants::eod,
            "Candidate result variables were not published to the RModel copy.");
    validateVariable(candidate,R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE,1U,
                     candidate.getNNodes(),"Candidate nodal temperature shape is invalid.");
    validateVariable(candidate,R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,3U,
                     candidate.getNNodes(),"Candidate nodal displacement shape is invalid.");
    validateVariable(candidate,R_VARIABLE_STRESS_X,R_VARIABLE_APPLY_ELEMENT,1U,
                     candidate.getNElements(),"Candidate element stress shape is invalid.");
    for (RVariableType stressType : {R_VARIABLE_STRESS_Y,R_VARIABLE_STRESS_Z,
                                     R_VARIABLE_STRESS_XY,R_VARIABLE_STRESS_XZ,
                                     R_VARIABLE_STRESS_YZ,R_VARIABLE_STRESS_VON_MISES,
                                     R_VARIABLE_STRESS_NORMAL,R_VARIABLE_STRESS_SHEAR})
    {
        require(candidate.findVariable(stressType) != RConstants::eod,
                "Candidate element stress output was not published to the RModel copy.");
        validateVariable(candidate,stressType,R_VARIABLE_APPLY_ELEMENT,1U,
                         candidate.getNElements(),"Candidate element stress output shape is invalid.");
    }

    for (uint node=0;node<candidate.getNNodes();node++)
    {
        require(candidate.getVariable(candidateTemperaturePosition).getValue(0,node) ==
                    result.nextState.temperatureKelvin[node] &&
                    candidate.getVariable(candidateDisplacementPosition).getValue(0,node) ==
                    result.nextState.axialDisplacementMetres[node],
                "Candidate nodal output differs from the certified monolithic result.");
    }
    for (uint element=0;element<candidate.getNElements();element++)
    {
        require(candidate.getVariable(candidateStressPosition).getValue(0,element) ==
                    result.memberAxialStressPascals[element],
                "Candidate element stress differs from the certified monolithic result.");
        for (RVariableType stressType : {R_VARIABLE_STRESS_X,R_VARIABLE_STRESS_Y,
                                         R_VARIABLE_STRESS_Z,R_VARIABLE_STRESS_XY,
                                         R_VARIABLE_STRESS_XZ,R_VARIABLE_STRESS_YZ,
                                         R_VARIABLE_STRESS_VON_MISES,R_VARIABLE_STRESS_NORMAL,
                                         R_VARIABLE_STRESS_SHEAR})
        {
            const uint position = candidate.findVariable(stressType);
            require(std::isfinite(candidate.getVariable(position).getValue(0,element)),
                    "Candidate element stress output contains a non-finite value.");
        }
    }

    model = candidate;
    return result;
}
