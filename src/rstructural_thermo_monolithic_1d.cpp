#include "rstructural_thermo_monolithic_1d.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
constexpr std::size_t maximumNodeCount = 64;
constexpr std::size_t maximumMemberCount = 1024;
constexpr unsigned int maximumNewtonIterations = 100;
constexpr unsigned int maximumLineSearchSteps = 64;

[[noreturn]] void invalidInput(const std::string &message)
{
    throw std::invalid_argument("Monolithic thermoelastic 1-D step: "+message);
}

[[noreturn]] void failedStep(const std::string &message)
{
    throw std::runtime_error("Monolithic thermoelastic 1-D step: "+message);
}

void requireFinite(double value, const char *description)
{
    if (!std::isfinite(value))
    {
        invalidInput(std::string(description)+" must be finite.");
    }
}

double norm(const std::vector<double> &values, const char *description)
{
    double value = 0.0;
    for (double entry : values)
    {
        if (!std::isfinite(entry))
        {
            failedStep(std::string(description)+" contains a non-finite entry.");
        }
        value = std::hypot(value,entry);
    }
    if (!std::isfinite(value))
    {
        failedStep(std::string(description)+" norm is non-finite.");
    }
    return value;
}

double freeNorm(const std::vector<double> &residual,
                const std::vector<bool> &prescribed,
                const char *description)
{
    double value = 0.0;
    for (std::size_t i=0;i<residual.size();i++)
    {
        if (!prescribed[i])
        {
            if (!std::isfinite(residual[i]))
            {
                failedStep(std::string(description)+" contains a non-finite free residual.");
            }
            value = std::hypot(value,residual[i]);
        }
    }
    if (!std::isfinite(value))
    {
        failedStep(std::string(description)+" residual norm is non-finite.");
    }
    return value;
}

struct IndexedMember
{
    std::size_t node0 = 0;
    std::size_t node1 = 0;
    double length = 0.0;
};

struct ValidatedInput
{
    std::vector<IndexedMember> members;
    std::vector<std::size_t> freeUnknowns;
};

ValidatedInput validate(const RStructuralThermoMonolithic1DModel &model,
                         const RStructuralThermoMonolithic1DState &acceptedState,
                         const RStructuralThermoMonolithic1DConditions &conditions,
                         const RStructuralThermoMonolithic1DOptions &options)
{
    const std::size_t nNodes = model.nodeCoordinatesMetres.size();
    if (nNodes < 2 || nNodes > maximumNodeCount)
    {
        invalidInput("the TRUSS1 mesh must contain between 2 and 64 nodes.");
    }
    if (model.members.empty() || model.members.size() > maximumMemberCount)
    {
        invalidInput("the mesh must contain between 1 and 1024 TRUSS1 members.");
    }
    requireFinite(model.referenceTemperatureKelvin,"reference temperature");
    if (model.referenceTemperatureKelvin <= 0.0)
    {
        invalidInput("reference temperature must be positive absolute Kelvin.");
    }
    if (!std::isfinite(options.timeStepSeconds) || options.timeStepSeconds <= 0.0)
    {
        invalidInput("time step must be finite and positive.");
    }
    if (options.maximumNewtonIterations == 0 ||
        options.maximumNewtonIterations > maximumNewtonIterations ||
        options.maximumLineSearchSteps == 0 ||
        options.maximumLineSearchSteps > maximumLineSearchSteps)
    {
        invalidInput("Newton steps must be in [1,100] and line-search trials in [1,64].");
    }
    const double tolerances[] = {
        options.mechanicalAbsoluteResidualToleranceNewtons,
        options.mechanicalRelativeResidualTolerance,
        options.thermalAbsoluteResidualToleranceWatts,
        options.thermalRelativeResidualTolerance
    };
    for (double tolerance : tolerances)
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0)
        {
            invalidInput("residual tolerances must be finite and nonnegative.");
        }
    }
    if ((options.mechanicalAbsoluteResidualToleranceNewtons == 0.0 &&
         options.mechanicalRelativeResidualTolerance == 0.0) ||
        (options.thermalAbsoluteResidualToleranceWatts == 0.0 &&
         options.thermalRelativeResidualTolerance == 0.0))
    {
        invalidInput("each field requires a positive absolute or relative residual tolerance.");
    }

    auto requireNodeVectorSize = [nNodes](std::size_t size, const char *description) {
        if (size != nNodes)
        {
            invalidInput(std::string(description)+" must contain one value per node.");
        }
    };
    requireNodeVectorSize(acceptedState.temperatureKelvin.size(),"accepted temperatures");
    requireNodeVectorSize(acceptedState.axialDisplacementMetres.size(),"accepted displacements");
    requireNodeVectorSize(conditions.mechanicalForceNewtons.size(),"mechanical forces");
    requireNodeVectorSize(conditions.displacementPrescribed.size(),"displacement mask");
    requireNodeVectorSize(conditions.prescribedDisplacementMetres.size(),"prescribed displacements");
    requireNodeVectorSize(conditions.heatPowerWatts.size(),"heat powers");
    requireNodeVectorSize(conditions.temperaturePrescribed.size(),"temperature mask");
    requireNodeVectorSize(conditions.prescribedTemperatureKelvin.size(),"prescribed temperatures");

    ValidatedInput validated;
    validated.members.reserve(model.members.size());
    validated.freeUnknowns.reserve(2*nNodes);
    std::vector<bool> nodeConnected(nNodes,false);
    std::vector<std::vector<std::size_t>> adjacency(nNodes);
    for (std::size_t i=0;i<nNodes;i++)
    {
        const double coordinate = model.nodeCoordinatesMetres[i];
        const double temperature = acceptedState.temperatureKelvin[i];
        const double displacement = acceptedState.axialDisplacementMetres[i];
        requireFinite(coordinate,"node coordinate");
        requireFinite(temperature,"accepted temperature");
        requireFinite(displacement,"accepted displacement");
        requireFinite(conditions.mechanicalForceNewtons[i],"mechanical force");
        requireFinite(conditions.prescribedDisplacementMetres[i],"prescribed displacement");
        requireFinite(conditions.heatPowerWatts[i],"heat power");
        requireFinite(conditions.prescribedTemperatureKelvin[i],"prescribed temperature");
        if (temperature <= 0.0)
        {
            invalidInput("accepted temperatures must be positive absolute Kelvin.");
        }
        if (conditions.temperaturePrescribed[i] &&
            conditions.prescribedTemperatureKelvin[i] <= 0.0)
        {
            invalidInput("prescribed temperatures must be positive absolute Kelvin.");
        }
        if (!conditions.temperaturePrescribed[i])
        {
            validated.freeUnknowns.push_back(i);
        }
    }
    for (std::size_t i=0;i<nNodes;i++)
    {
        if (!conditions.displacementPrescribed[i])
        {
            validated.freeUnknowns.push_back(nNodes+i);
        }
    }

    for (const RStructuralThermoMonolithic1DMember &member : model.members)
    {
        if (member.node0 >= nNodes || member.node1 >= nNodes || member.node0 == member.node1)
        {
            invalidInput("each member must connect two distinct valid node indices.");
        }
        const double length = model.nodeCoordinatesMetres[member.node1]-
                              model.nodeCoordinatesMetres[member.node0];
        if (!std::isfinite(length) || length <= 0.0)
        {
            invalidInput("member connectivity must run in increasing global x with positive length.");
        }
        nodeConnected[member.node0] = true;
        nodeConnected[member.node1] = true;
        adjacency[member.node0].push_back(member.node1);
        adjacency[member.node1].push_back(member.node0);
        const double positiveProperties[] = {
            member.areaSquareMetres,
            member.youngModulusPascals,
            member.densityKilogramsPerCubicMetre,
            member.specificHeatJoulesPerKilogramKelvin
        };
        for (double property : positiveProperties)
        {
            if (!std::isfinite(property) || property <= 0.0)
            {
                invalidInput("area, Young modulus, density and specific heat must be finite and positive.");
            }
        }
        requireFinite(member.expansionPerKelvin,"thermal expansion coefficient");
        if (!std::isfinite(member.conductivityWattsPerMetreKelvin) ||
            member.conductivityWattsPerMetreKelvin < 0.0)
        {
            invalidInput("thermal conductivity must be finite and nonnegative.");
        }
        if (!std::isfinite(member.kelvinVoigtViscosityPascalSeconds) ||
            member.kelvinVoigtViscosityPascalSeconds < 0.0)
        {
            invalidInput("Kelvin-Voigt viscosity must be finite and nonnegative.");
        }
        validated.members.push_back({member.node0,member.node1,length});
    }
    if (std::find(nodeConnected.cbegin(),nodeConnected.cend(),false) != nodeConnected.cend())
    {
        invalidInput("every node must belong to at least one TRUSS1 member.");
    }
    std::vector<bool> visited(nNodes,false);
    for (std::size_t root=0;root<nNodes;root++)
    {
        if (visited[root])
        {
            continue;
        }
        std::vector<std::size_t> pending{root};
        bool anchored = false;
        visited[root] = true;
        while (!pending.empty())
        {
            const std::size_t node = pending.back();
            pending.pop_back();
            anchored = anchored || conditions.displacementPrescribed[node];
            for (std::size_t adjacent : adjacency[node])
            {
                if (!visited[adjacent])
                {
                    visited[adjacent] = true;
                    pending.push_back(adjacent);
                }
            }
        }
        if (!anchored)
        {
            invalidInput("each connected rod component requires at least one prescribed displacement to remove rigid translation.");
        }
    }
    return validated;
}

struct Evaluation
{
    std::vector<double> thermalResidual;
    std::vector<double> mechanicalResidual;
    std::vector<double> thermalStorageTerms;
    std::vector<double> thermalConductionTerms;
    std::vector<double> thermalHeaterTerms;
    std::vector<double> thermalReversibleTerms;
    std::vector<double> thermalDissipationTerms;
    std::vector<double> internalForce;
    std::vector<double> jacobian;
    std::vector<double> memberStrain;
    std::vector<double> memberStress;
    std::vector<double> memberAxialForce;
    std::vector<double> memberDissipatedEnergy;
    double mechanicalResidualNorm = 0.0;
    double thermalResidualNorm = 0.0;
    double mechanicalTarget = 0.0;
    double thermalTarget = 0.0;
    double thermalEnergyChange = 0.0;
    double heaterEnergy = 0.0;
    double reversibleEnergy = 0.0;
    double viscousDissipationEnergy = 0.0;
    double dissipationHeatEnergy = 0.0;
    double thermalBoundaryEnergy = 0.0;
    double thermalEnergyBalanceResidual = 0.0;
};

void addJacobian(Evaluation &evaluation, std::size_t dimension,
                 std::size_t row, std::size_t column, double value)
{
    const std::size_t index = row*dimension+column;
    evaluation.jacobian[index] += value;
    if (!std::isfinite(evaluation.jacobian[index]))
    {
        failedStep("coupled Jacobian coefficient overflowed.");
    }
}

Evaluation evaluate(const RStructuralThermoMonolithic1DModel &model,
                    const RStructuralThermoMonolithic1DState &acceptedState,
                    const RStructuralThermoMonolithic1DConditions &conditions,
                    const RStructuralThermoMonolithic1DOptions &options,
                    const ValidatedInput &validated,
                    const RStructuralThermoMonolithic1DState &candidate,
                    bool assembleJacobian)
{
    const std::size_t nNodes = model.nodeCoordinatesMetres.size();
    const std::size_t dimension = 2*nNodes;
    Evaluation evaluation;
    evaluation.thermalResidual.resize(nNodes,0.0);
    evaluation.mechanicalResidual.resize(nNodes,0.0);
    evaluation.thermalStorageTerms.resize(nNodes,0.0);
    evaluation.thermalConductionTerms.resize(nNodes,0.0);
    evaluation.thermalHeaterTerms.resize(nNodes,0.0);
    evaluation.thermalReversibleTerms.resize(nNodes,0.0);
    evaluation.thermalDissipationTerms.resize(nNodes,0.0);
    evaluation.internalForce.resize(nNodes,0.0);
    evaluation.memberStrain.resize(model.members.size(),0.0);
    evaluation.memberStress.resize(model.members.size(),0.0);
    evaluation.memberAxialForce.resize(model.members.size(),0.0);
    evaluation.memberDissipatedEnergy.resize(model.members.size(),0.0);
    if (assembleJacobian)
    {
        evaluation.jacobian.resize(dimension*dimension,0.0);
    }

    for (std::size_t i=0;i<nNodes;i++)
    {
        evaluation.mechanicalResidual[i] = -conditions.mechanicalForceNewtons[i];
        evaluation.thermalHeaterTerms[i] = -conditions.heatPowerWatts[i];
        evaluation.heaterEnergy += options.timeStepSeconds*conditions.heatPowerWatts[i];
    }

    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto &member = model.members[e];
        const auto &indexed = validated.members[e];
        const std::size_t i = indexed.node0;
        const std::size_t j = indexed.node1;
        const double length = indexed.length;
        const double volume = member.areaSquareMetres*length;
        const double volumetricCapacity = member.densityKilogramsPerCubicMetre*
            member.specificHeatJoulesPerKilogramKelvin;
        const double capacityScale = volumetricCapacity*volume/
            (6.0*options.timeStepSeconds);
        const double conductivityScale = member.conductivityWattsPerMetreKelvin*
            member.areaSquareMetres/length;
        const double temperatureIncrement0 = candidate.temperatureKelvin[i]-
                                               acceptedState.temperatureKelvin[i];
        const double temperatureIncrement1 = candidate.temperatureKelvin[j]-
                                               acceptedState.temperatureKelvin[j];
        const double storage0 = capacityScale*(2.0*temperatureIncrement0+
                                                 temperatureIncrement1);
        const double storage1 = capacityScale*(temperatureIncrement0+
                                                 2.0*temperatureIncrement1);
        const double conduction = conductivityScale*(candidate.temperatureKelvin[i]-
                                                       candidate.temperatureKelvin[j]);
        evaluation.thermalStorageTerms[i] += storage0;
        evaluation.thermalStorageTerms[j] += storage1;
        evaluation.thermalConductionTerms[i] += conduction;
        evaluation.thermalConductionTerms[j] -= conduction;
        evaluation.thermalResidual[i] += storage0+conduction;
        evaluation.thermalResidual[j] += storage1-conduction;

        if (assembleJacobian)
        {
            addJacobian(evaluation,dimension,i,j,capacityScale);
            addJacobian(evaluation,dimension,i,i,2.0*capacityScale+conductivityScale);
            addJacobian(evaluation,dimension,j,i,capacityScale);
            addJacobian(evaluation,dimension,j,j,2.0*capacityScale+conductivityScale);
            addJacobian(evaluation,dimension,i,j,-conductivityScale);
            addJacobian(evaluation,dimension,j,i,-conductivityScale);
        }

        const double displacementDifference = candidate.axialDisplacementMetres[j]-
                                              candidate.axialDisplacementMetres[i];
        const double acceptedDifference = acceptedState.axialDisplacementMetres[j]-
                                          acceptedState.axialDisplacementMetres[i];
        const double increment = displacementDifference-acceptedDifference;
        const double strain = displacementDifference/length;
        const double strainRate = increment/(length*options.timeStepSeconds);
        const double meanTemperature = 0.5*(candidate.temperatureKelvin[i]+
                                             candidate.temperatureKelvin[j]);
        const double stress = member.youngModulusPascals*(strain-
            member.expansionPerKelvin*(meanTemperature-model.referenceTemperatureKelvin))+
            member.kelvinVoigtViscosityPascalSeconds*strainRate;
        const double axialForce = member.areaSquareMetres*stress;
        evaluation.memberStrain[e] = strain;
        evaluation.memberStress[e] = stress;
        evaluation.memberAxialForce[e] = axialForce;
        evaluation.internalForce[i] -= axialForce;
        evaluation.internalForce[j] += axialForce;
        evaluation.mechanicalResidual[i] -= axialForce;
        evaluation.mechanicalResidual[j] += axialForce;

        if (assembleJacobian)
        {
            const double mechanicalStiffness = member.areaSquareMetres/
                length*(member.youngModulusPascals+
                         member.kelvinVoigtViscosityPascalSeconds/options.timeStepSeconds);
            const double thermalTangent = 0.5*member.areaSquareMetres*
                member.youngModulusPascals*member.expansionPerKelvin;
            addJacobian(evaluation,dimension,nNodes+i,nNodes+i,mechanicalStiffness);
            addJacobian(evaluation,dimension,nNodes+i,nNodes+j,-mechanicalStiffness);
            addJacobian(evaluation,dimension,nNodes+j,nNodes+i,-mechanicalStiffness);
            addJacobian(evaluation,dimension,nNodes+j,nNodes+j,mechanicalStiffness);
            addJacobian(evaluation,dimension,nNodes+i,i,thermalTangent);
            addJacobian(evaluation,dimension,nNodes+i,j,thermalTangent);
            addJacobian(evaluation,dimension,nNodes+j,i,-thermalTangent);
            addJacobian(evaluation,dimension,nNodes+j,j,-thermalTangent);
        }

        const double reversiblePower = model.referenceTemperatureKelvin*
            member.youngModulusPascals*member.expansionPerKelvin*
            member.areaSquareMetres*increment/(2.0*options.timeStepSeconds);
        evaluation.thermalReversibleTerms[i] += reversiblePower;
        evaluation.thermalReversibleTerms[j] += reversiblePower;
        evaluation.thermalResidual[i] += reversiblePower;
        evaluation.thermalResidual[j] += reversiblePower;
        evaluation.reversibleEnergy += 2.0*options.timeStepSeconds*reversiblePower;
        if (assembleJacobian)
        {
            const double reversibleTangent = model.referenceTemperatureKelvin*
                member.youngModulusPascals*member.expansionPerKelvin*
                member.areaSquareMetres/(2.0*options.timeStepSeconds);
            for (const std::size_t thermalNode : {i,j})
            {
                addJacobian(evaluation,dimension,thermalNode,nNodes+i,-reversibleTangent);
                addJacobian(evaluation,dimension,thermalNode,nNodes+j,reversibleTangent);
            }
        }

        const double dissipationPower = member.kelvinVoigtViscosityPascalSeconds*
            member.areaSquareMetres*length*strainRate*strainRate;
        const double dissipationEnergy = dissipationPower*options.timeStepSeconds;
        evaluation.memberDissipatedEnergy[e] = dissipationEnergy;
        evaluation.viscousDissipationEnergy += dissipationEnergy;
        const double nodalDissipationPower = 0.5*dissipationPower;
        const double appliedDissipationPower = options.heatKelvinVoigtDissipation
            ? nodalDissipationPower : 0.0;
        evaluation.thermalDissipationTerms[i] -= appliedDissipationPower;
        evaluation.thermalDissipationTerms[j] -= appliedDissipationPower;
        evaluation.thermalResidual[i] -= appliedDissipationPower;
        evaluation.thermalResidual[j] -= appliedDissipationPower;
        evaluation.dissipationHeatEnergy += options.timeStepSeconds*
            (2.0*appliedDissipationPower);
        if (assembleJacobian && options.heatKelvinVoigtDissipation &&
            member.kelvinVoigtViscosityPascalSeconds > 0.0)
        {
            const double dissipationTangent = member.kelvinVoigtViscosityPascalSeconds*
                member.areaSquareMetres*increment/
                (length*options.timeStepSeconds*options.timeStepSeconds);
            for (const std::size_t thermalNode : {i,j})
            {
                addJacobian(evaluation,dimension,thermalNode,nNodes+i,dissipationTangent);
                addJacobian(evaluation,dimension,thermalNode,nNodes+j,-dissipationTangent);
            }
        }

        const double averageTemperatureIncrement = 0.5*(temperatureIncrement0+
                                                          temperatureIncrement1);
        evaluation.thermalEnergyChange += volumetricCapacity*volume*
                                           averageTemperatureIncrement;
    }

    for (std::size_t i=0;i<nNodes;i++)
    {
        evaluation.thermalResidual[i] += evaluation.thermalHeaterTerms[i];
        if (conditions.temperaturePrescribed[i])
        {
            evaluation.thermalBoundaryEnergy += options.timeStepSeconds*
                evaluation.thermalResidual[i];
        }
    }

    evaluation.mechanicalResidualNorm = freeNorm(evaluation.mechanicalResidual,
        conditions.displacementPrescribed,"mechanical");
    evaluation.thermalResidualNorm = freeNorm(evaluation.thermalResidual,
        conditions.temperaturePrescribed,"thermal");
    const double mechanicalScale = norm(evaluation.internalForce,"internal force")+
        norm(conditions.mechanicalForceNewtons,"external force");
    const double thermalScale = norm(evaluation.thermalStorageTerms,"thermal storage")+
        norm(evaluation.thermalConductionTerms,"thermal conduction")+
        norm(evaluation.thermalHeaterTerms,"heater power")+
        norm(evaluation.thermalReversibleTerms,"reversible thermoelastic power")+
        norm(evaluation.thermalDissipationTerms,"dissipation power");
    evaluation.mechanicalTarget =
        options.mechanicalAbsoluteResidualToleranceNewtons+
        options.mechanicalRelativeResidualTolerance*mechanicalScale;
    evaluation.thermalTarget = options.thermalAbsoluteResidualToleranceWatts+
        options.thermalRelativeResidualTolerance*thermalScale;
    if (!std::isfinite(evaluation.mechanicalTarget) ||
        !std::isfinite(evaluation.thermalTarget))
    {
        failedStep("residual certification target is non-finite.");
    }
    evaluation.thermalEnergyBalanceResidual = evaluation.thermalEnergyChange+
        evaluation.reversibleEnergy-evaluation.heaterEnergy-
        evaluation.dissipationHeatEnergy-evaluation.thermalBoundaryEnergy;
    if (!std::isfinite(evaluation.thermalEnergyBalanceResidual) ||
        !std::isfinite(evaluation.heaterEnergy) ||
        !std::isfinite(evaluation.reversibleEnergy) ||
        !std::isfinite(evaluation.viscousDissipationEnergy) ||
        !std::isfinite(evaluation.dissipationHeatEnergy))
    {
        failedStep("thermal energy diagnostic overflowed.");
    }
    for (double residual : evaluation.thermalResidual)
    {
        if (!std::isfinite(residual))
        {
            failedStep("thermal residual contains a non-finite value.");
        }
    }
    for (double residual : evaluation.mechanicalResidual)
    {
        if (!std::isfinite(residual))
        {
            failedStep("mechanical residual contains a non-finite value.");
        }
    }
    return evaluation;
}

double normalizedResidual(double residual, double target)
{
    if (target > 0.0)
    {
        return residual/target;
    }
    return residual == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
}

double merit(const Evaluation &evaluation)
{
    return std::max(normalizedResidual(evaluation.mechanicalResidualNorm,
                                       evaluation.mechanicalTarget),
                    normalizedResidual(evaluation.thermalResidualNorm,
                                       evaluation.thermalTarget));
}

std::vector<double> solveEquilibratedDense(std::vector<double> matrix,
                                          std::vector<double> rhs)
{
    const std::size_t n = rhs.size();
    if (matrix.size() != n*n)
    {
        failedStep("coupled linear-system dimensions are inconsistent.");
    }
    if (n == 0)
    {
        return {};
    }

    // Equilibrate rows and columns before pivoting because displacement and
    // temperature equations have different units and coefficient scales.
    for (std::size_t row=0;row<n;row++)
    {
        double scale = 0.0;
        for (std::size_t column=0;column<n;column++)
        {
            scale = std::max(scale,std::fabs(matrix[row*n+column]));
        }
        if (!std::isfinite(scale))
        {
            failedStep("coupled Jacobian row scale is non-finite.");
        }
        if (scale == 0.0)
        {
            if (rhs[row] == 0.0)
            {
                for (std::size_t column=0;column<n;column++)
                {
                    if (matrix[column*n+row] != 0.0)
                    {
                        failedStep("coupled Jacobian contains an unconstrained zero row.");
                    }
                }
                matrix[row*n+row] = 1.0;
                rhs[row] = 0.0;
                continue;
            }
            failedStep("coupled Jacobian has a zero row with a nonzero residual.");
        }
        for (std::size_t column=0;column<n;column++)
        {
            matrix[row*n+column] /= scale;
        }
        rhs[row] /= scale;
    }

    std::vector<double> columnScales(n,0.0);
    for (std::size_t column=0;column<n;column++)
    {
        for (std::size_t row=0;row<n;row++)
        {
            columnScales[column] = std::max(columnScales[column],
                                             std::fabs(matrix[row*n+column]));
        }
        if (columnScales[column] == 0.0 || !std::isfinite(columnScales[column]))
        {
            failedStep("coupled Jacobian has an unconstrained column.");
        }
        for (std::size_t row=0;row<n;row++)
        {
            matrix[row*n+column] /= columnScales[column];
        }
    }

    const double pivotTolerance = 64.0*std::numeric_limits<double>::epsilon()*
                                  static_cast<double>(n);
    for (std::size_t column=0;column<n;column++)
    {
        std::size_t pivot = column;
        double pivotMagnitude = std::fabs(matrix[pivot*n+column]);
        for (std::size_t row=column+1;row<n;row++)
        {
            const double magnitude = std::fabs(matrix[row*n+column]);
            if (magnitude > pivotMagnitude)
            {
                pivot = row;
                pivotMagnitude = magnitude;
            }
        }
        if (!std::isfinite(pivotMagnitude) || pivotMagnitude <= pivotTolerance)
        {
            failedStep("coupled Jacobian is singular or numerically rank deficient.");
        }
        if (pivot != column)
        {
            for (std::size_t j=column;j<n;j++)
            {
                std::swap(matrix[column*n+j],matrix[pivot*n+j]);
            }
            std::swap(rhs[column],rhs[pivot]);
        }
        for (std::size_t row=column+1;row<n;row++)
        {
            const double factor = matrix[row*n+column]/matrix[column*n+column];
            matrix[row*n+column] = 0.0;
            for (std::size_t j=column+1;j<n;j++)
            {
                matrix[row*n+j] -= factor*matrix[column*n+j];
            }
            rhs[row] -= factor*rhs[column];
        }
    }

    std::vector<double> scaledSolution(n,0.0);
    for (std::size_t row=n;row-- > 0;)
    {
        double value = rhs[row];
        for (std::size_t column=row+1;column<n;column++)
        {
            value -= matrix[row*n+column]*scaledSolution[column];
        }
        scaledSolution[row] = value/matrix[row*n+row];
        if (!std::isfinite(scaledSolution[row]))
        {
            failedStep("coupled Newton correction is non-finite.");
        }
    }
    std::vector<double> solution(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        solution[i] = scaledSolution[i]/columnScales[i];
        if (!std::isfinite(solution[i]))
        {
            failedStep("coupled Newton correction is non-finite after scaling.");
        }
    }
    return solution;
}

bool allTemperaturesPositive(const RStructuralThermoMonolithic1DState &state)
{
    return std::all_of(state.temperatureKelvin.cbegin(),state.temperatureKelvin.cend(),
                       [](double value) { return std::isfinite(value) && value > 0.0; });
}

} // namespace

RStructuralThermoMonolithic1DResult RStructuralThermoMonolithic1D::step(
    const RStructuralThermoMonolithic1DModel &model,
    const RStructuralThermoMonolithic1DState &acceptedState,
    const RStructuralThermoMonolithic1DConditions &conditions,
    const RStructuralThermoMonolithic1DOptions &options)
{
    const ValidatedInput validated = validate(model,acceptedState,conditions,options);
    RStructuralThermoMonolithic1DState candidate = acceptedState;
    for (std::size_t i=0;i<model.nodeCoordinatesMetres.size();i++)
    {
        if (conditions.temperaturePrescribed[i])
        {
            candidate.temperatureKelvin[i] = conditions.prescribedTemperatureKelvin[i];
        }
        if (conditions.displacementPrescribed[i])
        {
            candidate.axialDisplacementMetres[i] = conditions.prescribedDisplacementMetres[i];
        }
    }
    if (!allTemperaturesPositive(candidate))
    {
        invalidInput("trial temperatures must be positive absolute Kelvin.");
    }

    unsigned int newtonIterations = 0;
    Evaluation current;
    bool converged = false;
    for (unsigned int iteration=0;iteration<=options.maximumNewtonIterations;iteration++)
    {
        current = evaluate(model,acceptedState,conditions,options,validated,candidate,true);
        if (merit(current) <= 1.0)
        {
            converged = true;
            break;
        }
        if (iteration == options.maximumNewtonIterations)
        {
            break;
        }

        const std::size_t fullDimension = 2*model.nodeCoordinatesMetres.size();
        const std::size_t reducedDimension = validated.freeUnknowns.size();
        std::vector<double> reducedMatrix(reducedDimension*reducedDimension,0.0);
        std::vector<double> reducedRightHandSide(reducedDimension,0.0);
        for (std::size_t row=0;row<reducedDimension;row++)
        {
            const std::size_t fullRow = validated.freeUnknowns[row];
            const double residual = fullRow < model.nodeCoordinatesMetres.size()
                ? current.thermalResidual[fullRow]
                : current.mechanicalResidual[fullRow-model.nodeCoordinatesMetres.size()];
            reducedRightHandSide[row] = -residual;
            for (std::size_t column=0;column<reducedDimension;column++)
            {
                reducedMatrix[row*reducedDimension+column] =
                    current.jacobian[fullRow*fullDimension+validated.freeUnknowns[column]];
            }
        }
        const std::vector<double> correction = solveEquilibratedDense(
            std::move(reducedMatrix),std::move(reducedRightHandSide));
        const double currentMerit = merit(current);
        bool acceptedLineSearch = false;
        double stepScale = 1.0;
        for (unsigned int search=0;search<options.maximumLineSearchSteps;search++)
        {
            RStructuralThermoMonolithic1DState trial = candidate;
            for (std::size_t i=0;i<reducedDimension;i++)
            {
                const std::size_t unknown = validated.freeUnknowns[i];
                if (unknown < model.nodeCoordinatesMetres.size())
                {
                    trial.temperatureKelvin[unknown] += stepScale*correction[i];
                }
                else
                {
                    trial.axialDisplacementMetres[
                        unknown-model.nodeCoordinatesMetres.size()] += stepScale*correction[i];
                }
            }
            if (allTemperaturesPositive(trial))
            {
                Evaluation trialEvaluation;
                try
                {
                    trialEvaluation = evaluate(model,acceptedState,conditions,
                        options,validated,trial,false);
                }
                catch (const std::runtime_error &)
                {
                    // A finite Newton correction can still produce an
                    // overflowing residual or energy diagnostic at this
                    // trial scale. Reject that trial and continue the bounded
                    // backtracking sequence; the next scale may be evaluable.
                    stepScale *= 0.5;
                    continue;
                }
                const double trialMerit = merit(trialEvaluation);
                if (trialMerit <= 1.0 ||
                    trialMerit < currentMerit*(1.0-1.0e-4*stepScale))
                {
                    candidate = std::move(trial);
                    acceptedLineSearch = true;
                    ++newtonIterations;
                    break;
                }
            }
            stepScale *= 0.5;
        }
        if (!acceptedLineSearch)
        {
            failedStep("Newton line search could not reduce the coupled residual while keeping temperature positive.");
        }
    }
    if (!converged)
    {
        failedStep("coupled Newton iteration did not converge within its configured limit.");
    }

    // Reassemble from the returned candidate rather than certifying from the
    // residual/Jacobian that generated the final Newton correction.
    const Evaluation certified = evaluate(model,acceptedState,conditions,options,
                                           validated,candidate,false);
    if (merit(certified) > 1.0 || !allTemperaturesPositive(candidate))
    {
        failedStep("fresh N/W residual or positive-temperature certification failed.");
    }

    RStructuralThermoMonolithic1DResult result;
    result.nextState = std::move(candidate);
    result.memberAxialStrain = certified.memberStrain;
    result.memberAxialStressPascals = certified.memberStress;
    result.memberAxialForceNewtons = certified.memberAxialForce;
    result.memberDissipatedEnergyJoules = certified.memberDissipatedEnergy;
    result.mechanicalResidualNormNewtons = certified.mechanicalResidualNorm;
    result.thermalResidualNormWatts = certified.thermalResidualNorm;
    result.thermalEnergyChangeJoules = certified.thermalEnergyChange;
    result.heaterEnergyJoules = certified.heaterEnergy;
    result.reversibleThermoelasticEnergyJoules = certified.reversibleEnergy;
    result.viscousDissipationEnergyJoules = certified.viscousDissipationEnergy;
    result.dissipationHeatEnergyJoules = certified.dissipationHeatEnergy;
    result.thermalBoundaryEnergyJoules = certified.thermalBoundaryEnergy;
    result.thermalEnergyBalanceResidualJoules = certified.thermalEnergyBalanceResidual;
    result.newtonIterations = newtonIterations;
    return result;
}
