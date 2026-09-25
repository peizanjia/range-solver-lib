#include "rstructural_thermo_coupled_1d.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include <QHash>
#include <QSet>

namespace
{
using SparseRow = std::map<std::size_t,double>;
using SparseRows = std::vector<SparseRow>;

struct DisjointSets
{
    explicit DisjointSets(std::size_t count) : parent(count),rank(count,0)
    {
        std::iota(parent.begin(),parent.end(),std::size_t(0));
    }

    std::size_t find(std::size_t value)
    {
        if (parent[value] != value)
            parent[value] = find(parent[value]);
        return parent[value];
    }

    void unite(std::size_t first,std::size_t second)
    {
        first = find(first);
        second = find(second);
        if (first == second)
            return;
        if (rank[first] < rank[second])
            std::swap(first,second);
        parent[second] = first;
        if (rank[first] == rank[second])
            ++rank[first];
    }

    std::vector<std::size_t> parent;
    std::vector<unsigned int> rank;
};

struct IndexedMember
{
    std::size_t node0 = 0;
    std::size_t node1 = 0;
};

struct ValidatedInput
{
    QHash<QString,std::size_t> nodeIndex;
    std::vector<IndexedMember> members;
    std::vector<double> oldTemperature;
    std::vector<double> oldDisplacement;
    std::vector<double> nodalCapacity;
    std::vector<double> heaterPower;
    std::vector<double> externalForce;
    std::vector<bool> displacementPrescribed;
    std::vector<double> prescribedDisplacement;
    std::vector<std::size_t> freeDofs;
    std::vector<std::size_t> globalToFree;
};

[[noreturn]] void invalidInput(const std::string &message)
{
    throw std::invalid_argument("Thermo-coupled 1-D step: " + message);
}

[[noreturn]] void failedStep(const std::string &message)
{
    throw std::runtime_error("Thermo-coupled 1-D step: " + message);
}

double stableNorm(const std::vector<double> &values,const char *description)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            failedStep(std::string(description) + " contains a non-finite value.");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        failedStep(std::string(description) + " norm is non-finite.");
    return norm;
}

double dot(const std::vector<double> &left,const std::vector<double> &right)
{
    if (left.size() != right.size())
        failedStep("internal linear-system vector sizes do not match.");
    double value = 0.0;
    for (std::size_t i=0;i<left.size();i++)
    {
        value += left[i]*right[i];
        if (!std::isfinite(value))
            failedStep("linear-system dot product is non-finite.");
    }
    return value;
}

void add(SparseRows &matrix,std::size_t row,std::size_t column,double value)
{
    if (!std::isfinite(value))
        failedStep("assembled matrix coefficient is non-finite.");
    matrix[row][column] += value;
    if (!std::isfinite(matrix[row][column]))
        failedStep("assembled matrix coefficient overflowed.");
}

std::vector<double> multiply(const SparseRows &matrix,const std::vector<double> &x)
{
    if (matrix.size() != x.size())
        failedStep("internal matrix-vector dimensions do not match.");
    std::vector<double> product(matrix.size(),0.0);
    for (std::size_t row=0;row<matrix.size();row++)
    {
        for (const auto &entry : matrix[row])
        {
            product[row] += entry.second*x[entry.first];
            if (!std::isfinite(product[row]))
                failedStep("linear-system matrix-vector product is non-finite.");
        }
    }
    return product;
}

std::vector<double> solveCg(const SparseRows &matrix,
                            const std::vector<double> &rhs,
                            unsigned int maxIterations,
                            double absoluteTolerance,
                            double relativeTolerance,
                            const char *systemName)
{
    const std::size_t n = rhs.size();
    if (matrix.size() != n)
        failedStep(std::string(systemName) + " matrix and right-hand side dimensions differ.");
    if (n == 0)
        return {};

    std::vector<double> diagonal(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        const auto position = matrix[i].find(i);
        if (position != matrix[i].end())
            diagonal[i] = position->second;
        if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0)
            failedStep(std::string(systemName) + " matrix has a nonpositive diagonal.");
    }

    std::vector<double> x(n,0.0);
    std::vector<double> residual(rhs);
    const double rhsNorm = stableNorm(rhs,"linear-system right-hand side");
    const double target = absoluteTolerance+relativeTolerance*rhsNorm;
    if (!std::isfinite(target) || target <= 0.0)
        failedStep(std::string(systemName) + " residual tolerance is invalid.");
    if (stableNorm(residual,"linear-system initial residual") <= target)
        return x;

    std::vector<double> preconditioned(n,0.0);
    std::vector<double> direction(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        preconditioned[i] = residual[i]/diagonal[i];
        direction[i] = preconditioned[i];
    }
    double residualPreconditioned = dot(residual,preconditioned);
    if (!std::isfinite(residualPreconditioned) || residualPreconditioned <= 0.0)
        failedStep(std::string(systemName) + " initial preconditioned residual is not positive.");

    bool converged = false;
    for (unsigned int iteration=0;iteration<maxIterations;iteration++)
    {
        const std::vector<double> product = multiply(matrix,direction);
        const double denominator = dot(direction,product);
        if (!std::isfinite(denominator) || denominator <= 0.0)
            failedStep(std::string(systemName) + " matrix is not positive definite on the free DOFs.");
        const double alpha = residualPreconditioned/denominator;
        if (!std::isfinite(alpha))
            failedStep(std::string(systemName) + " CG step length is non-finite.");
        for (std::size_t i=0;i<n;i++)
        {
            x[i] += alpha*direction[i];
            residual[i] -= alpha*product[i];
            if (!std::isfinite(x[i]) || !std::isfinite(residual[i]))
                failedStep(std::string(systemName) + " CG iterate is non-finite.");
        }
        if (stableNorm(residual,"linear-system CG residual") <= target)
        {
            converged = true;
            break;
        }
        for (std::size_t i=0;i<n;i++)
            preconditioned[i] = residual[i]/diagonal[i];
        const double nextResidualPreconditioned = dot(residual,preconditioned);
        if (!std::isfinite(nextResidualPreconditioned) || nextResidualPreconditioned <= 0.0)
            failedStep(std::string(systemName) + " CG residual lost positive definiteness.");
        const double beta = nextResidualPreconditioned/residualPreconditioned;
        if (!std::isfinite(beta))
            failedStep(std::string(systemName) + " CG direction is non-finite.");
        for (std::size_t i=0;i<n;i++)
            direction[i] = preconditioned[i]+beta*direction[i];
        residualPreconditioned = nextResidualPreconditioned;
    }
    if (!converged)
        failedStep(std::string(systemName) + " conjugate gradients did not converge.");

    // Certify the returned solution with an independently recomputed residual.
    std::vector<double> freshResidual = multiply(matrix,x);
    for (std::size_t i=0;i<n;i++)
        freshResidual[i] = rhs[i]-freshResidual[i];
    if (stableNorm(freshResidual,"fresh linear-system residual") > target)
        failedStep(std::string(systemName) + " fresh residual failed certification.");
    return x;
}

void validateNodeIds(const std::vector<QString> &ids,
                     QHash<QString,std::size_t> &index,
                     const char *description)
{
    if (ids.empty())
        invalidInput(std::string(description) + " requires at least one node.");
    index.reserve(qsizetype(ids.size()));
    for (std::size_t i=0;i<ids.size();i++)
    {
        const QString &id = ids[i];
        if (id.isEmpty() || id != id.trimmed() || index.contains(id))
            invalidInput(std::string(description) + " node IDs must be unique, nonempty and unpadded.");
        index.insert(id,i);
    }
}

double evaluateProperty(const std::function<double(double)> &property,
                        double temperature,
                        const char *name)
{
    if (!property)
        invalidInput(std::string(name) + " callback is missing.");
    double value = 0.0;
    try
    {
        value = property(temperature);
    }
    catch (...)
    {
        invalidInput(std::string(name) + " callback threw an exception.");
    }
    if (!std::isfinite(value))
        invalidInput(std::string(name) + " callback returned a non-finite value.");
    return value;
}

double gaussExpansion(const std::function<double(double)> &expansion,
                      double fromKelvin,
                      double toKelvin,
                      bool sixteenPoint)
{
    static constexpr double abscissae8[] = {
        0.18343464249564980494,0.52553240991632898582,
        0.79666647741362673959,0.96028985649753623168
    };
    static constexpr double weights8[] = {
        0.36268378337836198297,0.31370664587788728734,
        0.22238103445337447054,0.10122853629037625915
    };
    static constexpr double abscissae16[] = {
        0.095012509837637440185,0.28160355077925891323,
        0.45801677765722738634,0.61787624440264374845,
        0.75540440835500303390,0.86563120238783174388,
        0.94457502307323257608,0.98940093499164993260
    };
    static constexpr double weights16[] = {
        0.18945061045506849629,0.18260341504492358887,
        0.16915651939500253819,0.14959598881657673208,
        0.12462897125553387205,0.09515851168249278481,
        0.062253523938647892863,0.027152459411754094852
    };
    const double midpoint = 0.5*(fromKelvin+toKelvin);
    const double halfRange = 0.5*(toKelvin-fromKelvin);
    double sum = 0.0;
    const std::size_t count = sixteenPoint ? 8U : 4U;
    for (std::size_t i=0;i<count;i++)
    {
        const double abscissa = sixteenPoint ? abscissae16[i] : abscissae8[i];
        const double weight = sixteenPoint ? weights16[i] : weights8[i];
        const double offset = halfRange*abscissa;
        const double lower = evaluateProperty(expansion,midpoint-offset,"thermal expansion");
        const double upper = evaluateProperty(expansion,midpoint+offset,"thermal expansion");
        sum += weight*(lower+upper);
        if (!std::isfinite(sum))
            failedStep("integrated thermal expansion overflowed.");
    }
    const double result = halfRange*sum;
    if (!std::isfinite(result))
        failedStep("integrated thermal expansion is non-finite.");
    return result;
}

double adaptiveExpansion(const std::function<double(double)> &expansion,
                         double fromKelvin,
                         double toKelvin,
                         double absoluteTolerance,
                         double relativeTolerance,
                         unsigned int remainingDepth)
{
    const double estimate8 = gaussExpansion(expansion,fromKelvin,toKelvin,false);
    const double estimate16 = gaussExpansion(expansion,fromKelvin,toKelvin,true);
    const double errorEstimate = std::fabs(estimate16-estimate8);
    const double tolerance = absoluteTolerance+relativeTolerance*std::fabs(estimate16);
    if (!std::isfinite(errorEstimate) || !std::isfinite(tolerance))
        failedStep("adaptive thermal-expansion quadrature produced a non-finite estimate.");
    if (errorEstimate <= tolerance)
        return estimate16;
    if (remainingDepth == 0)
        failedStep("adaptive thermal-expansion quadrature did not converge.");
    const double midpoint = 0.5*(fromKelvin+toKelvin);
    if (midpoint == fromKelvin || midpoint == toKelvin)
        failedStep("adaptive thermal-expansion quadrature reached floating-point resolution.");
    const double left = adaptiveExpansion(expansion,fromKelvin,midpoint,
                                          0.5*absoluteTolerance,relativeTolerance,
                                          remainingDepth-1);
    const double right = adaptiveExpansion(expansion,midpoint,toKelvin,
                                           0.5*absoluteTolerance,relativeTolerance,
                                           remainingDepth-1);
    const double result = left+right;
    if (!std::isfinite(result))
        failedStep("adaptive thermal-expansion integral overflowed.");
    return result;
}

double integratedExpansion(const std::function<double(double)> &expansion,
                           double fromKelvin,
                           double toKelvin)
{
    // Compare 8- and 16-point Gauss-Legendre rules and subdivide when their
    // difference exceeds tolerance. Arbitrary non-smooth callbacks may fail
    // explicitly at the bounded recursion depth instead of silently passing.
    return adaptiveExpansion(expansion,fromKelvin,toKelvin,1.0e-14,1.0e-10,16);
}

ValidatedInput validate(const RStructuralThermoCoupled1DModel &model,
                        const std::vector<RStructuralThermoCoupled1DNodeState> &acceptedState,
                        const RStructuralThermoCoupled1DConditions &conditions,
                        const RStructuralThermoCoupled1DOptions &options)
{
    ValidatedInput data;
    validateNodeIds(model.nodeIds,data.nodeIndex,"Thermo-structural model");
    const std::size_t nNodes = model.nodeIds.size();
    if (acceptedState.size() != nNodes)
        invalidInput("accepted state must contain one record per physical node.");
    if (model.members.empty())
        invalidInput("at least one rod member is required.");
    if (!std::isfinite(model.referenceTemperatureKelvin) || model.referenceTemperatureKelvin <= 0.0)
        invalidInput("reference temperature must be finite and above absolute zero.");
    if (!std::isfinite(options.dt) || options.dt <= 0.0 ||
        options.maxLinearIterations == 0 ||
        (options.twoWayDissipationHeat && options.maxCouplingIterations == 0))
        invalidInput("time step and iteration limits must be positive and finite.");

    const double tolerances[] = {
        options.thermalAbsoluteResidualTolerance,options.thermalRelativeResidualTolerance,
        options.mechanicalAbsoluteResidualTolerance,options.mechanicalRelativeResidualTolerance,
        options.energyAbsoluteTolerance,options.energyRelativeTolerance,
        options.couplingPowerAbsoluteTolerance,options.couplingPowerRelativeTolerance
    };
    for (double tolerance : tolerances)
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0)
            invalidInput("all residual, energy and coupling tolerances must be finite and nonnegative.");
    }
    if (options.thermalAbsoluteResidualTolerance == 0.0 &&
        options.thermalRelativeResidualTolerance == 0.0)
        invalidInput("at least one thermal residual tolerance must be positive.");
    if (options.mechanicalAbsoluteResidualTolerance == 0.0 &&
        options.mechanicalRelativeResidualTolerance == 0.0)
        invalidInput("at least one mechanical residual tolerance must be positive.");
    if (options.energyAbsoluteTolerance == 0.0 && options.energyRelativeTolerance == 0.0)
        invalidInput("at least one energy-balance tolerance must be positive.");
    if (options.couplingPowerAbsoluteTolerance == 0.0 &&
        options.couplingPowerRelativeTolerance == 0.0)
        invalidInput("at least one coupling-power tolerance must be positive.");

    QHash<QString,std::size_t> stateIndex;
    stateIndex.reserve(qsizetype(nNodes));
    data.oldTemperature.assign(nNodes,0.0);
    data.oldDisplacement.assign(nNodes,0.0);
    for (const auto &node : acceptedState)
    {
        if (node.nodeId.isEmpty() || node.nodeId != node.nodeId.trimmed() ||
            !data.nodeIndex.contains(node.nodeId) || stateIndex.contains(node.nodeId))
            invalidInput("accepted-state node IDs must match the physical mesh exactly.");
        if (!std::isfinite(node.temperatureKelvin) || node.temperatureKelvin <= 0.0 ||
            !std::isfinite(node.axialDisplacement))
            invalidInput("accepted temperatures must be positive finite kelvin and displacements finite metres.");
        const std::size_t index = data.nodeIndex.value(node.nodeId);
        data.oldTemperature[index] = node.temperatureKelvin;
        data.oldDisplacement[index] = node.axialDisplacement;
        stateIndex.insert(node.nodeId,index);
    }
    if (stateIndex.size() != qsizetype(nNodes))
        invalidInput("accepted-state node IDs do not cover the physical mesh.");

    data.members.reserve(model.members.size());
    data.nodalCapacity.assign(nNodes,0.0);
    DisjointSets components(nNodes);
    std::vector<bool> nodeUsed(nNodes,false);
    QSet<QString> elementIds;
    for (const auto &member : model.members)
    {
        if (member.elementId.isEmpty() || member.elementId != member.elementId.trimmed() ||
            elementIds.contains(member.elementId))
            invalidInput("element IDs must be unique, nonempty and unpadded.");
        elementIds.insert(member.elementId);
        if (!data.nodeIndex.contains(member.node0Id) || !data.nodeIndex.contains(member.node1Id))
            invalidInput("a rod endpoint ID is absent from the model node set.");
        const std::size_t node0 = data.nodeIndex.value(member.node0Id);
        const std::size_t node1 = data.nodeIndex.value(member.node1Id);
        if (node0 == node1 || !std::isfinite(member.area) || member.area <= 0.0 ||
            !std::isfinite(member.length) || member.length <= 0.0 ||
            !std::isfinite(member.conductivity) || member.conductivity <= 0.0 ||
            !std::isfinite(member.density) || member.density <= 0.0 ||
            !std::isfinite(member.specificHeat) || member.specificHeat <= 0.0 ||
            !std::isfinite(member.internalViscosity) || member.internalViscosity < 0.0 ||
            !member.youngModulusAtKelvin || !member.expansionAtKelvin)
            invalidInput("rod geometry and thermal/mechanical coefficients are invalid.");
        const double capacity = member.density*member.specificHeat*member.area*member.length;
        const double conductance = member.conductivity*member.area/member.length;
        const double viscosityCoefficient = member.internalViscosity*member.area/member.length;
        if (!std::isfinite(capacity) || capacity <= 0.0 || !std::isfinite(conductance) || conductance <= 0.0 ||
            !std::isfinite(viscosityCoefficient) || viscosityCoefficient < 0.0)
            invalidInput("rod thermal capacity, conductance or viscous coefficient overflowed.");
        const double meanTemperature = 0.5*(data.oldTemperature[node0]+data.oldTemperature[node1]);
        const double youngModulus = evaluateProperty(member.youngModulusAtKelvin,meanTemperature,"Young modulus");
        if (youngModulus <= 0.0)
            invalidInput("Young modulus must be positive at the accepted temperature.");
        (void)evaluateProperty(member.expansionAtKelvin,meanTemperature,"thermal expansion");
        if (!std::isfinite(youngModulus*member.area/member.length) ||
            youngModulus*member.area/member.length <= 0.0)
            invalidInput("rod axial stiffness is nonpositive or non-finite.");

        data.members.push_back({node0,node1});
        data.nodalCapacity[node0] += 0.5*capacity;
        data.nodalCapacity[node1] += 0.5*capacity;
        if (!std::isfinite(data.nodalCapacity[node0]) || !std::isfinite(data.nodalCapacity[node1]))
            invalidInput("assembled nodal thermal capacity overflowed.");
        nodeUsed[node0] = true;
        nodeUsed[node1] = true;
        components.unite(node0,node1);
    }
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!nodeUsed[node] || !std::isfinite(data.nodalCapacity[node]) || data.nodalCapacity[node] <= 0.0)
            invalidInput("every physical node must belong to a rod with positive thermal capacity.");
    }

    data.heaterPower.assign(nNodes,0.0);
    QSet<QString> heaterNodes;
    for (const auto &heater : conditions.heaters)
    {
        if (!data.nodeIndex.contains(heater.nodeId))
            invalidInput("heater refers to an unknown physical node ID.");
        if (heaterNodes.contains(heater.nodeId))
            invalidInput("each physical node may have at most one heater command per step.");
        heaterNodes.insert(heater.nodeId);
        if (!std::isfinite(heater.powerWatts) || !std::isfinite(heater.maximumPowerWatts) ||
            heater.maximumPowerWatts < 0.0 || heater.powerWatts < 0.0 ||
            heater.powerWatts > heater.maximumPowerWatts)
            invalidInput("heater power must be finite, nonnegative and within its configured watt limit.");
        data.heaterPower[data.nodeIndex.value(heater.nodeId)] = heater.powerWatts;
    }

    data.externalForce.assign(nNodes,0.0);
    for (const auto &force : conditions.mechanicalForces)
    {
        if (!data.nodeIndex.contains(force.nodeId) || !std::isfinite(force.forceNewtons))
            invalidInput("mechanical load refers to an unknown node or is non-finite.");
        const std::size_t node = data.nodeIndex.value(force.nodeId);
        data.externalForce[node] += force.forceNewtons;
        if (!std::isfinite(data.externalForce[node]))
            invalidInput("combined nodal mechanical load overflowed.");
    }

    data.displacementPrescribed.assign(nNodes,false);
    data.prescribedDisplacement.assign(nNodes,0.0);
    QSet<QString> prescribedNodes;
    for (const auto &boundary : conditions.prescribedDisplacements)
    {
        if (!data.nodeIndex.contains(boundary.nodeId) || !std::isfinite(boundary.displacementMetres))
            invalidInput("prescribed displacement refers to an unknown node or is non-finite.");
        if (prescribedNodes.contains(boundary.nodeId))
            invalidInput("a physical node has duplicate prescribed-displacement records.");
        prescribedNodes.insert(boundary.nodeId);
        const std::size_t node = data.nodeIndex.value(boundary.nodeId);
        data.displacementPrescribed[node] = true;
        data.prescribedDisplacement[node] = boundary.displacementMetres;
    }
    QSet<std::size_t> anchoredComponents;
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (data.displacementPrescribed[node])
            anchoredComponents.insert(components.find(node));
    }
    for (std::size_t node=0;node<nNodes;node++)
    {
        if (!anchoredComponents.contains(components.find(node)))
            invalidInput("every mechanical connected component requires a prescribed-displacement gauge.");
        if (!data.displacementPrescribed[node])
        {
            data.globalToFree.push_back(data.freeDofs.size());
            data.freeDofs.push_back(node);
        }
        else
        {
            data.globalToFree.push_back(std::numeric_limits<std::size_t>::max());
        }
    }
    return data;
}

struct ThermalSolution
{
    std::vector<double> temperature;
    double residualNorm = 0.0;
};

ThermalSolution solveThermal(const RStructuralThermoCoupled1DModel &model,
                             const ValidatedInput &data,
                             const std::vector<double> &dissipationPower,
                             const RStructuralThermoCoupled1DOptions &options)
{
    const std::size_t nNodes = model.nodeIds.size();
    SparseRows matrix(nNodes);
    std::vector<double> rhs(nNodes,0.0);
    for (std::size_t node=0;node<nNodes;node++)
    {
        const double capacityRate = data.nodalCapacity[node]/options.dt;
        const double source = data.heaterPower[node]+dissipationPower[node];
        if (!std::isfinite(capacityRate) || capacityRate <= 0.0 || !std::isfinite(source))
            failedStep("thermal capacity rate or source power is invalid.");
        add(matrix,node,node,capacityRate);
        rhs[node] = capacityRate*data.oldTemperature[node]+source;
        if (!std::isfinite(rhs[node]))
            failedStep("thermal right-hand side overflowed.");
    }
    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto &member = model.members[e];
        const auto [node0,node1] = data.members[e];
        const double conductance = member.conductivity*member.area/member.length;
        add(matrix,node0,node0,conductance);
        add(matrix,node1,node1,conductance);
        add(matrix,node0,node1,-conductance);
        add(matrix,node1,node0,-conductance);
    }
    ThermalSolution solution;
    solution.temperature = solveCg(matrix,rhs,options.maxLinearIterations,
                                   options.thermalAbsoluteResidualTolerance,
                                   options.thermalRelativeResidualTolerance,"thermal");
    for (double temperature : solution.temperature)
    {
        if (!std::isfinite(temperature) || temperature <= 0.0)
            failedStep("backward-Euler solve produced a nonpositive absolute temperature.");
    }
    std::vector<double> freshResidual = multiply(matrix,solution.temperature);
    std::vector<double> lhs = freshResidual;
    for (std::size_t i=0;i<nNodes;i++)
        freshResidual[i] -= rhs[i];
    solution.residualNorm = stableNorm(freshResidual,"fresh thermal residual");
    const double physicalScale = std::max(stableNorm(rhs,"thermal right-hand side"),
                                          stableNorm(lhs,"thermal left-hand side"));
    const double tolerance = options.thermalAbsoluteResidualTolerance+
                             options.thermalRelativeResidualTolerance*physicalScale;
    if (!std::isfinite(tolerance) || solution.residualNorm > tolerance)
        failedStep("fresh thermal residual failed its W-based certificate.");
    return solution;
}

struct MechanicalSolution
{
    std::vector<double> displacement;
    std::vector<double> internalForce;
    std::vector<double> elasticInternalForce;
    std::vector<double> viscousInternalForce;
    std::vector<double> memberForce;
    std::vector<double> elasticMemberForce;
    std::vector<double> viscousMemberForce;
    double elasticStrainEnergy = 0.0;
    double residualNorm = 0.0;
};

MechanicalSolution solveMechanical(const RStructuralThermoCoupled1DModel &model,
                                   const ValidatedInput &data,
                                   const std::vector<double> &temperature,
                                   const RStructuralThermoCoupled1DOptions &options)
{
    const std::size_t nNodes = model.nodeIds.size();
    const std::size_t nFree = data.freeDofs.size();
    SparseRows matrix(nFree);
    std::vector<double> rhs(nFree,0.0);
    for (std::size_t free=0;free<nFree;free++)
        rhs[free] = data.externalForce[data.freeDofs[free]];

    std::vector<double> thermalStrains(model.members.size(),0.0);
    std::vector<double> axialStiffness(model.members.size(),0.0);
    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto &member = model.members[e];
        const auto [node0,node1] = data.members[e];
        const double meanTemperature = 0.5*(temperature[node0]+temperature[node1]);
        const double young = evaluateProperty(member.youngModulusAtKelvin,meanTemperature,"Young modulus");
        if (young <= 0.0)
            failedStep("temperature-dependent Young modulus is not positive.");
        const double stiffness = young*member.area/member.length;
        const double expansion = integratedExpansion(member.expansionAtKelvin,
                                                       model.referenceTemperatureKelvin,
                                                       meanTemperature);
        const double dampingCoefficient = member.internalViscosity*member.area/member.length;
        const double dampingRate = dampingCoefficient/options.dt;
        const double effectiveStiffness = stiffness+dampingRate;
        const double previousRelativeDisplacement = data.oldDisplacement[node1]-data.oldDisplacement[node0];
        const double equivalentForce = stiffness*member.length*expansion+
                                      dampingRate*previousRelativeDisplacement;
        if (!std::isfinite(stiffness) || stiffness <= 0.0 || !std::isfinite(expansion) ||
            !std::isfinite(dampingCoefficient) || dampingCoefficient < 0.0 ||
            !std::isfinite(dampingRate) || dampingRate < 0.0 ||
            !std::isfinite(effectiveStiffness) || effectiveStiffness <= 0.0 ||
            !std::isfinite(equivalentForce))
            failedStep("temperature-dependent member stiffness or free thermal strain is invalid.");
        thermalStrains[e] = expansion;
        axialStiffness[e] = stiffness;

        const std::size_t free0 = data.globalToFree[node0];
        const std::size_t free1 = data.globalToFree[node1];
        if (free0 != std::numeric_limits<std::size_t>::max())
        {
            add(matrix,free0,free0,effectiveStiffness);
            rhs[free0] -= equivalentForce;
            if (free1 == std::numeric_limits<std::size_t>::max())
                rhs[free0] += effectiveStiffness*data.prescribedDisplacement[node1];
            else
                add(matrix,free0,free1,-effectiveStiffness);
        }
        if (free1 != std::numeric_limits<std::size_t>::max())
        {
            add(matrix,free1,free1,effectiveStiffness);
            rhs[free1] += equivalentForce;
            if (free0 == std::numeric_limits<std::size_t>::max())
                rhs[free1] += effectiveStiffness*data.prescribedDisplacement[node0];
            else
                add(matrix,free1,free0,-effectiveStiffness);
        }
    }

    MechanicalSolution solution;
    solution.displacement = data.prescribedDisplacement;
    const std::vector<double> freeDisplacement = solveCg(
        matrix,rhs,options.maxLinearIterations,
        options.mechanicalAbsoluteResidualTolerance,
        options.mechanicalRelativeResidualTolerance,"mechanical");
    for (std::size_t free=0;free<nFree;free++)
        solution.displacement[data.freeDofs[free]] = freeDisplacement[free];

    solution.internalForce.assign(nNodes,0.0);
    solution.elasticInternalForce.assign(nNodes,0.0);
    solution.viscousInternalForce.assign(nNodes,0.0);
    solution.memberForce.resize(model.members.size(),0.0);
    solution.elasticMemberForce.resize(model.members.size(),0.0);
    solution.viscousMemberForce.resize(model.members.size(),0.0);
    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto [node0,node1] = data.members[e];
        const auto &member = model.members[e];
        const double mechanicalExtension = solution.displacement[node1]-solution.displacement[node0]-
                                           member.length*thermalStrains[e];
        const double displacementIncrement =
            (solution.displacement[node1]-solution.displacement[node0])-
            (data.oldDisplacement[node1]-data.oldDisplacement[node0]);
        const double elasticForce = axialStiffness[e]*mechanicalExtension;
        const double viscousCoefficient = member.internalViscosity*member.area/member.length;
        const double viscousForce = viscousCoefficient*displacementIncrement/options.dt;
        const double axialForce = elasticForce+viscousForce;
        const double elasticEnergy = 0.5*axialStiffness[e]*mechanicalExtension*mechanicalExtension;
        if (!std::isfinite(mechanicalExtension) || !std::isfinite(displacementIncrement) ||
            !std::isfinite(elasticForce) || !std::isfinite(viscousForce) || !std::isfinite(axialForce) ||
            !std::isfinite(elasticEnergy) || elasticEnergy < 0.0)
            failedStep("recovered axial member force is non-finite.");
        solution.memberForce[e] = axialForce;
        solution.elasticMemberForce[e] = elasticForce;
        solution.viscousMemberForce[e] = viscousForce;
        solution.elasticStrainEnergy += elasticEnergy;
        if (!std::isfinite(solution.elasticStrainEnergy))
            failedStep("elastic strain energy overflowed.");
        solution.elasticInternalForce[node0] -= elasticForce;
        solution.elasticInternalForce[node1] += elasticForce;
        solution.viscousInternalForce[node0] -= viscousForce;
        solution.viscousInternalForce[node1] += viscousForce;
        solution.internalForce[node0] -= axialForce;
        solution.internalForce[node1] += axialForce;
        if (!std::isfinite(solution.elasticInternalForce[node0]) ||
            !std::isfinite(solution.elasticInternalForce[node1]) ||
            !std::isfinite(solution.viscousInternalForce[node0]) ||
            !std::isfinite(solution.viscousInternalForce[node1]) ||
            !std::isfinite(solution.internalForce[node0]) || !std::isfinite(solution.internalForce[node1]))
            failedStep("assembled nodal internal force is non-finite.");
    }

    std::vector<double> freeResidual;
    std::vector<double> freeExternal;
    std::vector<double> freeInternal;
    freeResidual.reserve(nFree);
    freeExternal.reserve(nFree);
    freeInternal.reserve(nFree);
    for (std::size_t node : data.freeDofs)
    {
        freeResidual.push_back(solution.internalForce[node]-data.externalForce[node]);
        freeExternal.push_back(data.externalForce[node]);
        freeInternal.push_back(solution.internalForce[node]);
    }
    solution.residualNorm = stableNorm(freeResidual,"free mechanical force residual");
    const double forceScale = std::max(stableNorm(freeExternal,"free mechanical loads"),
                                       stableNorm(freeInternal,"free mechanical internal forces"));
    const double tolerance = options.mechanicalAbsoluteResidualTolerance+
                             options.mechanicalRelativeResidualTolerance*forceScale;
    if (!std::isfinite(tolerance) || solution.residualNorm > tolerance)
        failedStep("fresh mechanical residual failed its N-based certificate.");
    return solution;
}

struct ElasticState
{
    std::vector<double> internalForce;
    double strainEnergy = 0.0;
};

ElasticState recoverElasticState(const RStructuralThermoCoupled1DModel &model,
                                 const ValidatedInput &data,
                                 const std::vector<double> &temperature,
                                 const std::vector<double> &displacement)
{
    ElasticState state;
    state.internalForce.assign(model.nodeIds.size(),0.0);
    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto &member = model.members[e];
        const auto [node0,node1] = data.members[e];
        const double meanTemperature = 0.5*(temperature[node0]+temperature[node1]);
        const double young = evaluateProperty(member.youngModulusAtKelvin,meanTemperature,"Young modulus");
        if (young <= 0.0)
            failedStep("temperature-dependent Young modulus is not positive.");
        const double thermalStrain = integratedExpansion(member.expansionAtKelvin,
                                                          model.referenceTemperatureKelvin,
                                                          meanTemperature);
        const double stiffness = young*member.area/member.length;
        const double mechanicalExtension = displacement[node1]-displacement[node0]-
                                           member.length*thermalStrain;
        const double force = stiffness*mechanicalExtension;
        const double energy = 0.5*stiffness*mechanicalExtension*mechanicalExtension;
        if (!std::isfinite(stiffness) || stiffness <= 0.0 || !std::isfinite(mechanicalExtension) ||
            !std::isfinite(force) || !std::isfinite(energy) || energy < 0.0)
            failedStep("recovered accepted elastic state is non-finite.");
        state.internalForce[node0] -= force;
        state.internalForce[node1] += force;
        state.strainEnergy += energy;
        if (!std::isfinite(state.internalForce[node0]) || !std::isfinite(state.internalForce[node1]) ||
            !std::isfinite(state.strainEnergy))
            failedStep("accepted elastic force or strain energy overflowed.");
    }
    return state;
}

std::vector<double> nodalDissipationPower(const RStructuralThermoCoupled1DModel &model,
                                          const ValidatedInput &data,
                                          const MechanicalSolution &mechanical,
                                          const ValidatedInput &accepted,
                                          double dt,
                                          std::vector<double> &memberEnergy,
                                          std::vector<double> &memberPower)
{
    std::vector<double> nodalPower(model.nodeIds.size(),0.0);
    memberEnergy.assign(model.members.size(),0.0);
    memberPower.assign(model.members.size(),0.0);
    for (std::size_t e=0;e<model.members.size();e++)
    {
        const auto &member = model.members[e];
        const auto [node0,node1] = data.members[e];
        const double extensionIncrement =
            (mechanical.displacement[node1]-mechanical.displacement[node0])-
            (accepted.oldDisplacement[node1]-accepted.oldDisplacement[node0]);
        const double coefficient = member.internalViscosity*member.area/member.length;
        const double power = coefficient*extensionIncrement*extensionIncrement/(dt*dt);
        const double energy = power*dt;
        if (!std::isfinite(extensionIncrement) || !std::isfinite(coefficient) || coefficient < 0.0 ||
            !std::isfinite(power) || !std::isfinite(energy) || power < 0.0 || energy < 0.0)
            failedStep("physical Kelvin-Voigt dissipation must be finite and nonnegative.");
        // Independently recover work from the viscous force and increment.
        const double viscousForce = coefficient*extensionIncrement/dt;
        const double signedWork = viscousForce*extensionIncrement;
        const double workTolerance = 64.0*std::numeric_limits<double>::epsilon()*
                                     std::max(energy,std::numeric_limits<double>::min());
        if (!std::isfinite(viscousForce) || !std::isfinite(signedWork) || signedWork < 0.0 ||
            std::fabs(signedWork-energy) > workTolerance)
            failedStep("independent viscous-force work check failed its nonnegative N m certificate.");
        memberPower[e] = power;
        memberEnergy[e] = energy;
        nodalPower[node0] += 0.5*power;
        nodalPower[node1] += 0.5*power;
        if (!std::isfinite(nodalPower[node0]) || !std::isfinite(nodalPower[node1]))
            failedStep("assembled nodal dissipation heat overflowed.");
    }
    return nodalPower;
}

}

RStructuralThermoCoupled1DResult RStructuralThermoCoupled1D::step(
    const RStructuralThermoCoupled1DModel &model,
    const std::vector<RStructuralThermoCoupled1DNodeState> &acceptedState,
    const RStructuralThermoCoupled1DConditions &conditions,
    const RStructuralThermoCoupled1DOptions &options)
{
    const ValidatedInput data = validate(model,acceptedState,conditions,options);
    const std::size_t nNodes = model.nodeIds.size();
    std::vector<double> guessedNodalDissipation(nNodes,0.0);
    ThermalSolution thermal;
    MechanicalSolution mechanical;
    std::vector<double> memberDissipationEnergy;
    std::vector<double> memberDissipationPower;
    std::vector<double> actualNodalDissipation;
    double couplingResidual = 0.0;
    unsigned int couplingIterations = 0;
    bool couplingConverged = false;

    const unsigned int maximumIterations = options.twoWayDissipationHeat
        ? options.maxCouplingIterations : 1U;
    for (unsigned int iteration=0;iteration<maximumIterations;iteration++)
    {
        thermal = solveThermal(model,data,guessedNodalDissipation,options);
        mechanical = solveMechanical(model,data,thermal.temperature,options);
        actualNodalDissipation = nodalDissipationPower(
            model,data,mechanical,data,options.dt,memberDissipationEnergy,memberDissipationPower);
        couplingIterations = iteration+1;

        if (!options.twoWayDissipationHeat)
        {
            couplingResidual = 0.0;
            couplingConverged = true;
            break;
        }

        std::vector<double> powerDifference(nNodes,0.0);
        for (std::size_t node=0;node<nNodes;node++)
            powerDifference[node] = actualNodalDissipation[node]-guessedNodalDissipation[node];
        couplingResidual = stableNorm(powerDifference,"two-way dissipation-power residual");
        const double powerScale = std::max(stableNorm(actualNodalDissipation,"physical dissipation power"),
                                           stableNorm(guessedNodalDissipation,"iterated dissipation power"));
        const double powerTolerance = options.couplingPowerAbsoluteTolerance+
                                      options.couplingPowerRelativeTolerance*powerScale;
        if (!std::isfinite(powerTolerance))
            failedStep("two-way coupling-power tolerance overflowed.");
        if (couplingResidual <= powerTolerance)
        {
            couplingConverged = true;
            break;
        }
        guessedNodalDissipation = actualNodalDissipation;
    }
    if (!couplingConverged)
        failedStep("two-way dissipation/temperature iteration did not converge.");

    RStructuralThermoCoupled1DResult result;
    result.thermalResidualNormWatts = thermal.residualNorm;
    result.mechanicalResidualNormNewtons = mechanical.residualNorm;
    result.couplingPowerResidualWatts = couplingResidual;
    result.couplingIterations = couplingIterations;
    result.memberAxialForcesNewtons = mechanical.memberForce;
    result.memberElasticAxialForcesNewtons = mechanical.elasticMemberForce;
    result.memberViscousAxialForcesNewtons = mechanical.viscousMemberForce;
    result.memberDissipatedEnergyJoules = memberDissipationEnergy;
    result.memberDissipationPowerWatts = memberDissipationPower;
    const ElasticState oldElastic = recoverElasticState(model,data,data.oldTemperature,data.oldDisplacement);
    result.elasticStrainEnergyChangeJoules = mechanical.elasticStrainEnergy-oldElastic.strainEnergy;
    result.nodes.reserve(nNodes);
    result.nextState.reserve(nNodes);
    std::vector<double> thermalEnergyIncrements(nNodes,0.0);
    std::vector<double> displacementIncrements(nNodes,0.0);
    result.mechanicalExternalWorkJoules = 0.0;
    result.mechanicalSupportWorkJoules = 0.0;
    for (std::size_t node=0;node<nNodes;node++)
    {
        result.nextState.push_back({model.nodeIds[node],thermal.temperature[node],mechanical.displacement[node]});
        displacementIncrements[node] = mechanical.displacement[node]-data.oldDisplacement[node];
        const double internalForce = mechanical.internalForce[node];
        const double support = data.displacementPrescribed[node]
            ? internalForce-data.externalForce[node] : 0.0;
        const double externalWork = data.externalForce[node]*displacementIncrements[node];
        if (!std::isfinite(support) || !std::isfinite(displacementIncrements[node]) ||
            !std::isfinite(externalWork))
            failedStep("prescribed-displacement support reaction is non-finite.");
        result.nodes.push_back({model.nodeIds[node],internalForce,support});
        result.mechanicalExternalWorkJoules += externalWork;
        if (data.displacementPrescribed[node])
        {
            const double oldElasticReaction = oldElastic.internalForce[node]-data.externalForce[node];
            const double newElasticReaction = mechanical.elasticInternalForce[node]-data.externalForce[node];
            const double elasticSupportWork = 0.5*(oldElasticReaction+newElasticReaction)*
                                             displacementIncrements[node];
            const double viscousSupportWork = mechanical.viscousInternalForce[node]*
                                              displacementIncrements[node];
            result.mechanicalSupportWorkJoules += elasticSupportWork+viscousSupportWork;
        }
        if (!std::isfinite(result.mechanicalExternalWorkJoules) ||
            !std::isfinite(result.mechanicalSupportWorkJoules))
            failedStep("mechanical external or support work overflowed.");
        thermalEnergyIncrements[node] = data.nodalCapacity[node]*
                                       (thermal.temperature[node]-data.oldTemperature[node]);
        if (!std::isfinite(thermalEnergyIncrements[node]))
            failedStep("nodal thermal energy increment is non-finite.");
    }
    result.thermalEnergyChangeJoules = 0.0;
    for (double energy : thermalEnergyIncrements)
    {
        result.thermalEnergyChangeJoules += energy;
        if (!std::isfinite(result.thermalEnergyChangeJoules))
            failedStep("total thermal energy change overflowed.");
    }
    result.heaterEnergyJoules = 0.0;
    for (double power : data.heaterPower)
    {
        result.heaterEnergyJoules += power*options.dt;
        if (!std::isfinite(result.heaterEnergyJoules))
            failedStep("bounded heater energy overflowed.");
    }
    result.viscousDissipationEnergyJoules = 0.0;
    for (double energy : memberDissipationEnergy)
    {
        if (energy < 0.0 || !std::isfinite(energy))
            failedStep("member viscous work is negative or non-finite.");
        result.viscousDissipationEnergyJoules += energy;
        if (!std::isfinite(result.viscousDissipationEnergyJoules))
            failedStep("total physical viscous work overflowed.");
    }
    const double thermalDissipationEnergy = options.twoWayDissipationHeat
        ? result.viscousDissipationEnergyJoules : 0.0;
    result.dissipationHeatEnergyJoules = thermalDissipationEnergy;
    result.thermalEnergyBalanceResidualJoules = result.thermalEnergyChangeJoules-
                                                result.heaterEnergyJoules-
                                                thermalDissipationEnergy;
    const double energyScale = std::max(std::fabs(result.thermalEnergyChangeJoules),
                                        result.heaterEnergyJoules+thermalDissipationEnergy);
    const double energyTolerance = options.energyAbsoluteTolerance+
                                   options.energyRelativeTolerance*energyScale+
                                   options.dt*couplingResidual;
    if (!std::isfinite(result.thermalEnergyBalanceResidualJoules) ||
        !std::isfinite(energyTolerance) ||
        std::fabs(result.thermalEnergyBalanceResidualJoules) > energyTolerance)
        failedStep("thermal energy balance failed its independent J=N m certificate.");
    result.mechanicalWorkBalanceResidualJoules = result.elasticStrainEnergyChangeJoules-
                                                 result.mechanicalExternalWorkJoules-
                                                 result.mechanicalSupportWorkJoules+
                                                 result.viscousDissipationEnergyJoules;
    if (!std::isfinite(result.elasticStrainEnergyChangeJoules) ||
        !std::isfinite(result.mechanicalWorkBalanceResidualJoules))
        failedStep("mechanical work diagnostic is non-finite.");
    return result;
}
