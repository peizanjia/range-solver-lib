#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

#include <omp.h>

#include "rsolverheat.h"
#include "rsolverfluidheat.h"
#include "rconvection.h"
#include "rmatrixsolver.h"

namespace
{

double stableNorm(const RRVector &values)
{
    double norm = 0.0;
    for (uint i=0;i<values.getNRows();i++)
    {
        if (!std::isfinite(values[i]))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        norm = std::hypot(norm,values[i]);
    }
    return norm;
}

void requireFiniteVector(const RRVector &values, const char *name)
{
    for (uint i=0;i<values.getNRows();i++)
    {
        if (!std::isfinite(values[i]))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Heat solve contains a non-finite %s value at index %u.",name,i);
        }
    }
}

double referenceMassJacobian(const RElement &element,
                             const std::vector<RNode> &referenceNodes,
                             uint integrationPoint,
                             uint elementId)
{
    const RElementShapeFunction &shapeFunction =
        RElement::getShapeFunction(element.getType(),integrationPoint);
    const RRMatrix &dN = shapeFunction.getDN();
    double characteristicLength = 0.0;
    for (uint i=0;i<element.size();i++)
    {
        const RNode &first = referenceNodes[element.getNodeId(i)];
        if (!std::isfinite(first.getX()) || !std::isfinite(first.getY()) ||
            !std::isfinite(first.getZ()))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Reference thermal mass has non-finite geometry on element %u.",
                         elementId);
        }
        for (uint j=i+1;j<element.size();j++)
        {
            const RNode &second = referenceNodes[element.getNodeId(j)];
            const double distance = std::hypot(first.getX()-second.getX(),
                                               first.getY()-second.getY(),
                                               first.getZ()-second.getZ());
            if (!std::isfinite(distance))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Reference thermal mass has non-finite geometry on element %u.",
                             elementId);
            }
            characteristicLength = std::max(characteristicLength,distance);
        }
    }

    const uint dimension = R_ELEMENT_TYPE_IS_SURFACE(element.getType()) ? 2U : 3U;
    const double scale = std::pow(characteristicLength,double(dimension));
    const double minimumJacobian = 128.0*std::numeric_limits<double>::epsilon()*scale;
    if (!std::isfinite(scale) || characteristicLength <= 0.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Reference thermal mass has degenerate or invalid geometry on element %u.",
                     elementId);
    }

    // Compute the geometric measure directly in global coordinates. Surface
    // elements use the tangent cross product, avoiding an arbitrary corner
    // frame; volumes use the signed 3D isoparametric Jacobian.
    double detJ = 0.0;
    const RNode &origin = referenceNodes[element.getNodeId(0)];
    if (R_ELEMENT_TYPE_IS_SURFACE(element.getType()))
    {
        RR3Vector firstTangent(0.0,0.0,0.0);
        RR3Vector secondTangent(0.0,0.0,0.0);
        for (uint node=0;node<element.size();node++)
        {
            const RNode &position = referenceNodes[element.getNodeId(node)];
            firstTangent[0] += dN[node][0]*(position.getX()-origin.getX());
            firstTangent[1] += dN[node][0]*(position.getY()-origin.getY());
            firstTangent[2] += dN[node][0]*(position.getZ()-origin.getZ());
            secondTangent[0] += dN[node][1]*(position.getX()-origin.getX());
            secondTangent[1] += dN[node][1]*(position.getY()-origin.getY());
            secondTangent[2] += dN[node][1]*(position.getZ()-origin.getZ());
        }
        const double crossX = firstTangent[1]*secondTangent[2]-firstTangent[2]*secondTangent[1];
        const double crossY = firstTangent[2]*secondTangent[0]-firstTangent[0]*secondTangent[2];
        const double crossZ = firstTangent[0]*secondTangent[1]-firstTangent[1]*secondTangent[0];
        detJ = std::hypot(crossX,crossY,crossZ);
    }
    else
    {
        RRMatrix jacobian(dimension,dimension);
        jacobian.fill(0.0);
        for (uint naturalAxis=0;naturalAxis<dimension;naturalAxis++)
        {
            for (uint physicalAxis=0;physicalAxis<dimension;physicalAxis++)
            {
                for (uint node=0;node<element.size();node++)
                {
                    const RNode &position = referenceNodes[element.getNodeId(node)];
                    const double coordinate = physicalAxis == 0
                        ? position.getX()-origin.getX()
                        : physicalAxis == 1 ? position.getY()-origin.getY()
                        : position.getZ()-origin.getZ();
                    jacobian[naturalAxis][physicalAxis] += dN[node][naturalAxis]*coordinate;
                }
            }
        }
        detJ = jacobian.getDeterminant();
    }
    if (!std::isfinite(detJ) || detJ <= minimumJacobian)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Reference thermal mass has degenerate or invalid geometry on element %u.",
                     elementId);
    }
    return detJ;
}

}

RSolverHeat::RSolverHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , solveCertified(false)
{
    this->problemType = R_PROBLEM_HEAT;
}

RSolverHeat::~RSolverHeat()
{
}

bool RSolverHeat::hasConverged() const
{
    return this->solveCertified;
}

void RSolverHeat::setControlledPointHeat(
    const std::vector<std::pair<uint,double>> &powers)
{
    std::vector<bool> seen(this->pModel->getNElements(),false);
    for (const auto &power : powers)
    {
        if (power.first >= this->pModel->getNElements() ||
            this->pModel->getElement(power.first).getType() != R_ELEMENT_POINT ||
            !std::isfinite(power.second) || power.second < 0.0 || seen[power.first])
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Controlled heat requires unique POINT elements and finite nonnegative W inputs.");
        }
        seen[power.first] = true;
    }
    this->controlledPointHeat = powers;
}

void RSolverHeat::setReferenceLineThermalMass(bool enabled)
{
    this->referenceLineThermalMass = enabled;
}

void RSolverHeat::setReferenceSurfaceVolumeThermalMass(bool enabled)
{
    this->referenceSurfaceVolumeThermalMass = enabled;
}

void RSolverHeat::findComputableElements(RProblemType problemType)
{
    RSolverGeneric::findComputableElements(problemType);
    if (problemType == R_PROBLEM_HEAT)
        for (const auto &power : this->controlledPointHeat)
            this->computableElements[power.first] = true;
}

void RSolverHeat::validateRunConfiguration() const
{
    // RSolverGeneric calls this before every run, before recovery or assembly
    // can fail. A previous successful solve must not make a later failed run
    // appear converged to an enclosing problem-task iteration.
    this->solveCertified = false;

    const double tolerance = this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getSolverCvgValue();
    if (!std::isfinite(tolerance) || tolerance < 0.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Heat CG convergence tolerance must be finite and nonnegative.");
    }
}

void RSolverHeat::recoverSharedData()
{
    this->RSolverGeneric::recoverSharedData();

    // Forced convection walls are driven by the fluid heat solver result. It is
    // kept apart from the shared element temperature, which this solver
    // overwrites over the whole mesh once it has run.
    this->fluidNodeTemperature.clear();
    if (this->pSharedData->hasData(RSolverFluidHeat::fluidNodeTemperatureKey,this->pModel->getNNodes()))
    {
        this->fluidNodeTemperature = this->pSharedData->findData(RSolverFluidHeat::fluidNodeTemperatureKey);
    }

    this->fluidNodeVelocity.clear();
    if (this->pSharedData->hasData(RSolverFluidHeat::fluidNodeVelocityKey,this->pModel->getNNodes()))
    {
        this->fluidNodeVelocity = this->pSharedData->findData(RSolverFluidHeat::fluidNodeVelocityKey);
    }
}

bool RSolverHeat::findFluidTemperature(uint elementId, double &fluidTemperature) const
{
    if (this->fluidNodeTemperature.size() != this->pModel->getNNodes())
    {
        return false;
    }
    if (elementId >= this->fluidElements.size())
    {
        return false;
    }

    uint fluidElementID = this->fluidElements[elementId];
    if (fluidElementID == RConstants::eod)
    {
        return false;
    }

    // The bulk temperature of the fluid element behind the wall - not the wall
    // node temperature, which the conformal mesh shares with the solid.
    const RElement &rFluidElement = this->pModel->getElement(fluidElementID);
    double temperature = 0.0;
    for (uint i=0;i<rFluidElement.size();i++)
    {
        temperature += this->fluidNodeTemperature[rFluidElement.getNodeId(i)];
    }
    fluidTemperature = temperature / double(rFluidElement.size());

    return true;
}

bool RSolverHeat::findFluidVelocity(uint elementId, double &fluidVelocity) const
{
    if (this->fluidNodeVelocity.size() != this->pModel->getNNodes())
    {
        return false;
    }
    if (elementId >= this->fluidElements.size())
    {
        return false;
    }

    uint fluidElementID = this->fluidElements[elementId];
    if (fluidElementID == RConstants::eod)
    {
        return false;
    }

    const RElement &rFluidElement = this->pModel->getElement(fluidElementID);
    double velocity = 0.0;
    for (uint i=0;i<rFluidElement.size();i++)
    {
        velocity += this->fluidNodeVelocity[rFluidElement.getNodeId(i)];
    }
    fluidVelocity = velocity / double(rFluidElement.size());

    return true;
}

void RSolverHeat::findFluidElements()
{
    this->fluidElements.resize(this->pModel->getNElements());
    this->fluidElements.fill(RConstants::eod);

    // An element belongs to the fluid domain when its material carries the
    // properties the fluid heat solver needs - the dynamic viscosity being the
    // one a solid does not have.
    QList<RMaterialProperty::Type> propList;
    for (uint i=0;i<RMaterialProperty::nTypes;i++)
    {
        if (R_PROBLEM_FLUID_HEAT & RMaterialProperty::getProblemTypeMask(RMaterialProperty::Type(i)))
        {
            propList.push_back(RMaterialProperty::Type(i));
        }
    }

    // Index fluid volume elements by node so the search below stays local.
    std::vector<std::vector<uint>> nodeToFluidElements(this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup->getMaterial().hasProperties(propList))
        {
            continue;
        }
        for (uint j=0;j<pElementGroup->size();j++)
        {
            uint elementID = pElementGroup->get(j);
            const RElement &rElement = this->pModel->getElement(elementID);
            if (!R_ELEMENT_TYPE_IS_VOLUME(rElement.getType()))
            {
                continue;
            }
            for (uint k=0;k<rElement.size();k++)
            {
                nodeToFluidElements[rElement.getNodeId(k)].push_back(elementID);
            }
        }
    }

    // A wall is a surface element whose nodes are all shared with a fluid volume
    // element - the mesh is conformal, so the fluid element on the other side of
    // the wall contains the whole face.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);

        for (uint j=0;j<rSurface.size();j++)
        {
            uint elementID = rSurface.get(j);
            const RElement &rElement = this->pModel->getElement(elementID);

            for (uint fluidElementID : nodeToFluidElements[rElement.getNodeId(0)])
            {
                const RElement &rFluidElement = this->pModel->getElement(fluidElementID);
                uint nNodesFound = 0;
                for (uint k=0;k<rElement.size();k++)
                {
                    if (rFluidElement.hasNodeId(rElement.getNodeId(k)))
                    {
                        nNodesFound++;
                    }
                }
                if (nNodesFound == rElement.size())
                {
                    this->fluidElements[elementID] = fluidElementID;
                    break;
                }
            }
        }
    }
}

double RSolverHeat::findTemperatureScale() const
{
    return 1.0;
}

void RSolverHeat::initialize()
{
}

void RSolverHeat::updateScales()
{
    this->scales.setMetre(this->findMeshScale());
    this->scales.setKelvin(this->findTemperatureScale());
}

void RSolverHeat::recover()
{
    this->recoverVariable(R_VARIABLE_TEMPERATURE,
                          R_VARIABLE_APPLY_NODE,
                          this->pModel->getNNodes(),
                          0,
                          this->nodeTemperature,
                          RVariable::getInitValue(R_VARIABLE_TEMPERATURE));
    this->recoverVariable(R_VARIABLE_HEAT,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementHeat,
                          0.0);
    this->recoverVariable(R_VARIABLE_HEAT_RADIATION,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementRadiativeHeat,
                          0.0);
    this->recoverVariable(R_VARIABLE_JOULE_HEAT,
                          R_VARIABLE_APPLY_ELEMENT,
                          this->pModel->getNElements(),
                          0,
                          this->elementJouleHeat,
                          0.0);
}

void RSolverHeat::prepare()
{
    const bool timeSolverEnabled = this->pModel->getTimeSolver().getEnabled();
    const bool useReferenceSurfaceVolumeMass = timeSolverEnabled &&
        this->referenceSurfaceVolumeThermalMass;

    // RSolverGeneric::run applies the stored displacement before a heat solve.
    // Reconstruct the reference coordinates once, while the model is in the
    // same scaled coordinate system used by the current-geometry assembly.
    std::vector<RNode> referenceNodes;
    if (useReferenceSurfaceVolumeMass)
    {
        referenceNodes = this->pModel->getNodes();
        const uint displacementPosition = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
        if (displacementPosition != RConstants::eod)
        {
            const RVariable &displacement = this->pModel->getVariable(displacementPosition);
            if (displacement.getApplyType() != R_VARIABLE_APPLY_NODE ||
                displacement.getNValues() != this->pModel->getNNodes() ||
                displacement.getNVectors() != 3)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Reference thermal mass requires three-component nodal displacement data.");
            }
            for (uint node=0;node<this->pModel->getNNodes();node++)
            {
                const double ux = displacement.getValue(0,node);
                const double uy = displacement.getValue(1,node);
                const double uz = displacement.getValue(2,node);
                if (!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz))
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Reference thermal mass has non-finite displacement at node %u.",node);
                }
                referenceNodes[node].move(RR3Vector(-ux,-uy,-uz));
            }
        }
    }

    RBVector temperatureSetValues;
    RBVector heatSetValues;
    RBVector heatRateAreaSetValues;
    RBVector heatRateVolumeSetValues;

    if (this->fluidElements.size() != this->pModel->getNElements())
    {
        this->findFluidElements();
    }

    this->generateNodeBook(R_PROBLEM_HEAT);

//    this->pModel->convertNodeToElementVector(this->nodeTemperature,this->elementTemperature);

    this->generateVariableVector(R_VARIABLE_TEMPERATURE,this->elementTemperature,temperatureSetValues,true,this->firstRun,this->firstRun);
    // The Heat boundary condition prescribes the total heat input for the whole
    // entity - it is spread over the entity measure to give the source density
    // the assembly below integrates.
    this->generateHeatVector(this->elementHeat,heatSetValues);
    for (const auto &power : this->controlledPointHeat)
    {
        this->elementHeat[power.first] += power.second;
        if (!std::isfinite(this->elementHeat[power.first]))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Controlled point heat plus existing heat input is non-finite.");
    }
    // Heat rate per unit area / volume are boundary conditions only - they are
    // applied to surface and volume elements respectively, where the source
    // term is already integrated over the element measure.
    this->generateVariableVector(R_VARIABLE_HEAT_RATE_AREA,this->elementHeatRateArea,heatRateAreaSetValues,true,false,false);
    this->generateVariableVector(R_VARIABLE_HEAT_RATE_VOLUME,this->elementHeatRateVolume,heatRateVolumeSetValues,true,false,false);
    this->generateMaterialVecor(RMaterialProperty::ThermalConductivity,this->elementConduction);
    this->generateMaterialVecor(RMaterialProperty::HeatCapacity,this->elementCapacity);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);

    // Validate and cache the reference mass integration before entering any
    // OpenMP element loops. Invalid reference geometry then fails the entire
    // candidate solve before matrix assembly can partially proceed.
    std::vector<std::vector<double>> referenceMassJacobians(this->pModel->getNElements());
    if (useReferenceSurfaceVolumeMass)
    {
        auto cacheGroupJacobians = [&](const RElementGroup &group) {
            for (uint i=0;i<group.size();i++)
            {
                const uint elementId = group.get(i);
                if (!this->computableElements[elementId])
                {
                    continue;
                }
                const RElement &element = this->pModel->getElement(elementId);
                const uint nIntegrationPoints = RElement::getNIntegrationPoints(element.getType());
                std::vector<double> &jacobians = referenceMassJacobians[elementId];
                jacobians.resize(nIntegrationPoints);
                for (uint point=0;point<nIntegrationPoints;point++)
                {
                    jacobians[point] = referenceMassJacobian(element,referenceNodes,point,elementId);
                }
            }
        };
        for (uint i=0;i<this->pModel->getNSurfaces();i++)
        {
            cacheGroupJacobians(this->pModel->getSurface(i));
        }
        for (uint i=0;i<this->pModel->getNVolumes();i++)
        {
            cacheGroupJacobians(this->pModel->getVolume(i));
        }
    }

    this->pModel->convertElementToNodeVector(this->elementTemperature,temperatureSetValues,this->nodeTemperature,true);

    uint nEnabled = this->nodeBook.getNEnabled();

    this->b.resize(nEnabled);
    this->x.resize(nEnabled);

    this->A.clear();
    this->A.setNRows(nEnabled);
    this->b.fill(0.0);
    this->x.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b once at the end.
    int np = omp_get_max_threads();
    std::vector<RSparseMatrix> Ap(np);
    std::vector<RRVector> bp(np);
    for (int t=0;t<np;t++)
    {
        Ap[t].setNRows(nEnabled);
        bp[t].resize(nEnabled);
        bp[t].fill(0.0);
    }

    // Prepare point elements.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        RPoint &point = this->pModel->getPoint(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(point.size());j++)
        {
            try
            {
                uint elementID = point.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_POINT(element.getType()));
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                // POINT has one node and no shape-function integration points.
                // Its heat condition is a total W input (not a density).
                if (timeSolverEnabled)
                    Me[0][0] = this->elementDensity[elementID] *
                               this->elementCapacity[elementID] * point.getVolume();
                fe[0] = this->elementHeat[elementID] + this->elementJouleHeat[elementID];
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
            }
            catch (const RError &rError)
            {
                #pragma omp critical
                {
                    RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
                    abort = true;
                }
            }
        }
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            try
            {
                uint elementID = line.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());

                RRMatrix B(element.size(),1);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                double referenceMassRatio = 1.0;
                if (timeSolverEnabled && this->referenceLineThermalMass)
                {
                    if (element.getType() != R_ELEMENT_TRUSS1 || element.size() != 2)
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                     "Reference thermal mass requires two-node TRUSS1 members.");
                    const auto &first = this->pModel->getNode(element.getNodeId(0));
                    const auto &second = this->pModel->getNode(element.getNodeId(1));
                    double referenceDx = second.getX()-first.getX();
                    double referenceDy = second.getY()-first.getY();
                    double referenceDz = second.getZ()-first.getZ();
                    const double currentLength = std::hypot(referenceDx,referenceDy,referenceDz);
                    const uint displacementPosition = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
                    if (displacementPosition != RConstants::eod)
                    {
                        const RVariable &displacement = this->pModel->getVariable(displacementPosition);
                        referenceDx -= displacement.getValue(0,element.getNodeId(1))-
                                       displacement.getValue(0,element.getNodeId(0));
                        referenceDy -= displacement.getValue(1,element.getNodeId(1))-
                                       displacement.getValue(1,element.getNodeId(0));
                        referenceDz -= displacement.getValue(2,element.getNodeId(1))-
                                       displacement.getValue(2,element.getNodeId(0));
                    }
                    const double referenceLength = std::hypot(referenceDx,referenceDy,referenceDz);
                    if (!std::isfinite(currentLength) || !std::isfinite(referenceLength) ||
                        currentLength <= 0.0 || referenceLength <= 0.0)
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                     "Thermal rod has invalid current or reference length.");
                    referenceMassRatio = referenceLength/currentLength;
                }

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += dN[m][0]*J[0][0];
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            double kcnd = (B[m][0]*B[n][0]) * line.getCrossArea() * this->elementConduction[elementID];

                            Ke[m][n] += kcnd * detJ * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * detJ * referenceMassRatio
                                         * shapeFunc.getW()
                                         * line.getCrossArea();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
            }
            catch (const RError &rError)
            {
                #pragma omp critical
                {
                    RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
                    abort = true;
                }
            }
        }
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

        double surfaceHtc = 0.0;
        double surfaceHtt = 0.0;

        this->getSimpleConvection(surface,surfaceHtc,surfaceHtt);
        this->reportForcedConvection(surface);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            try
            {
                uint elementID = surface.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),2);

                // htc/htt must be per-iteration - the correlated conditions
                // overwrite them per element and the loop runs in parallel.
                double htc = surfaceHtc;
                double htt = surfaceHtt;
                this->getForcedConvection(surface,elementID,htc,htt);
                this->getNaturalConvection(surface,elementID,htc,htt);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    const double massDetJ = useReferenceSurfaceVolumeMass
                        ? referenceMassJacobians[elementID][k] : detJ;

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]);
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            double kcnd = (B[m][0]*B[n][0]+B[m][1]*B[n][1]) * surface.getThickness() * this->elementConduction[elementID];
                            // Convection
                            double kcnv = N[m] * N[n] * htc;

                            Ke[m][n] += (kcnd + kcnv) * detJ * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * massDetJ
                                         * shapeFunc.getW()
                                         * surface.getThickness();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementHeatRateArea[elementID] + this->elementRadiativeHeat[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }

                // Convection force
                double elementArea = 0.0;
                if (element.findArea(this->pModel->getNodes(),elementArea))
                {
                    for (unsigned m=0;m<element.size();m++)
                    {
                        fe[m] += htc * htt * elementArea / element.size();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
            }
            catch (const RError &rError)
            {
                #pragma omp critical
                {
                    RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
                    abort = true;
                }
            }
        }
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Prepare volume elements.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        RVolume &volume = this->pModel->getVolume(i);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(volume.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = volume.get(j);

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size(),element.size());
                RRMatrix Ke(element.size(),element.size());
                RRVector fe(element.size());
                RRMatrix B(element.size(),3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                // Conduction
                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);
                    const double massDetJ = useReferenceSurfaceVolumeMass
                        ? referenceMassJacobians[elementID][k] : detJ;

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (unsigned m=0;m<element.size();m++)
                    {
                        for (unsigned n=0;n<element.size();n++)
                        {
                            // Conduction
                            Ke[m][n] += (B[m][0]*B[n][0] + B[m][1]*B[n][1] + B[m][2]*B[n][2])
                                     * this->elementConduction[elementID]
                                     * detJ
                                     * shapeFunc.getW();

                            // Mass
                            if (timeSolverEnabled)
                            {
                                Me[m][n] += N[m] * N[n]
                                         * this->elementDensity[elementID]
                                         * this->elementCapacity[elementID]
                                         * massDetJ
                                         * shapeFunc.getW();
                            }
                        }
                        // Force
                        fe[m] += (this->elementHeat[elementID] + this->elementHeatRateVolume[elementID] + this->elementJouleHeat[elementID]) * N[m] * detJ * shapeFunc.getW();
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())]);
            }
            catch (const RError &rError)
            {
                #pragma omp critical
                {
                    RLogger::error("%s\n",rError.getMessage().toUtf8().constData());
                    abort = true;
                }
            }
        }
        if (abort)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to prepare matrix system.");
        }
    }

    // Merge per-thread assembly buffers.
    #pragma omp parallel for default(shared)
    for (int64_t i=0;i<int64_t(this->A.getNRows());i++)
    {
        for (int t=0;t<np;t++)
        {
            this->A.getVector(uint(i)).addVector(Ap[t].getVector(uint(i)));
            this->b[uint(i)] += bp[t][uint(i)];
        }
    }
}

void RSolverHeat::solve()
{
    this->solveCertified = false;
    try
    {
        RLogger::indent();
        RMatrixSolver matrixSolver(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
        matrixSolver.solve(this->A,this->b,this->x,R_MATRIX_PRECONDITIONER_JACOBI,1);
        RLogger::unindent();
    }
    catch (const RError &)
    {
        RLogger::unindent();
        throw;
    }

    const uint nUnknowns = this->b.getNRows();
    if (this->A.getNRows() != nUnknowns || this->x.getNRows() != nUnknowns)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Heat solve returned inconsistent matrix or vector dimensions.");
    }
    requireFiniteVector(this->b,"right-hand side");
    requireFiniteVector(this->x,"solution");
    for (uint row=0;row<this->A.getNRows();row++)
    {
        const RSparseVector<double> &entries = this->A.getVector(row);
        for (uint entry=0;entry<entries.size();entry++)
        {
            if (entries.getIndex(entry) >= nUnknowns || !std::isfinite(entries.getValue(entry)))
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Heat solve matrix contains an invalid entry at row %u.",row);
            }
        }
    }

    RRVector product;
    RSparseMatrix::mlt(this->A,this->x,product);
    if (product.getNRows() != nUnknowns)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Heat solve residual has an invalid dimension.");
    }

    RRVector residual(nUnknowns,0.0);
    for (uint i=0;i<nUnknowns;i++)
    {
        if (!std::isfinite(product[i]))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Heat solve produced a non-finite matrix product at index %u.",i);
        }
        residual[i] = product[i]-this->b[i];
    }

    const double residualNorm = stableNorm(residual);
    const double rhsNorm = stableNorm(this->b);
    const double configuredTolerance = this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getSolverCvgValue();
    // Certify the requested relative tolerance.  The separate machine-scale
    // allowance covers rounding in the independently recomputed A*x-b.
    const double roundingAllowance = 64.0*std::numeric_limits<double>::epsilon()
                                     *std::max(1.0,rhsNorm);
    const double allowedResidual = configuredTolerance*rhsNorm+roundingAllowance;
    if (!std::isfinite(residualNorm) || !std::isfinite(rhsNorm) || !std::isfinite(allowedResidual))
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Heat solve produced a non-finite residual or convergence threshold.");
    }
    if (residualNorm > allowedResidual)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Heat solve did not converge (residual %g, allowed %g).",
                     residualNorm,allowedResidual);
    }

    if (this->nodeTemperature.getNRows() != this->pModel->getNNodes())
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Heat solve has an invalid nodal temperature vector.");
    }
    RRVector candidateNodeTemperature(this->nodeTemperature);

    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        if (this->nodeBook.getValue(i,position))
        {
            if (position >= nUnknowns)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Heat solve node book contains an invalid unknown index at node %u.",i);
            }
            candidateNodeTemperature[i] = this->x[position];
        }
    }

    requireFiniteVector(candidateNodeTemperature,"nodal temperature");
    RRVector candidateElementTemperature;
    this->pModel->convertNodeToElementVector(candidateNodeTemperature,candidateElementTemperature);
    requireFiniteVector(candidateElementTemperature,"element temperature");

    this->nodeTemperature = candidateNodeTemperature;
    this->elementTemperature = candidateElementTemperature;
    this->solveCertified = true;
}

void RSolverHeat::process()
{
//    RLogger::info("Processing results\n");
    double Qx;
    double Qy;
    double Qz;

    // Initialize heat flux vector vector
    this->elementHeatFlux.resize(this->pModel->getNElements(),RR3Vector(0.0,0.0,0.0));

    // Initialize heat transfer coefficient vector. It stays zero on every
    // element which carries no convection boundary condition.
    this->elementHeatTransferCoefficient.resize(this->pModel->getNElements());
    this->elementHeatTransferCoefficient.fill(0.0);

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);

        for (uint j=0;j<line.size();j++)
        {
            uint elementID = line.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRVector B(element.size());

            Qx = Qy = Qz = 0.0;

            // Conduction
            if (line.getCrossArea() > 0.0)
            {
                double Qi = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    // Heat flux is a density - the shape function derivatives
                    // are averaged over the integration points and must not be
                    // weighted by the Jacobian determinant.
                    this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    if (line.getCrossArea() != 0.0)
                    {
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m] += dN[m][0] * J[0][0] / double(nInp);
                        }
                    }
                }

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    Qi -= B[k] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                Qx += R[0][0]*Qi;
                Qy += R[1][0]*Qi;
                Qz += R[2][0]*Qi;
            }

            this->elementHeatFlux[elementID][0] = Qx;
            this->elementHeatFlux[elementID][1] = Qy;
            this->elementHeatFlux[elementID][2] = Qz;
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);

        double surfaceHtc = 0.0;
        double surfaceHtt = 0.0;

        this->getSimpleConvection(surface,surfaceHtc,surfaceHtt);

        for (uint j=0;j<surface.size();j++)
        {
            uint elementID = surface.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRMatrix B(element.size(),2);

            double htc = surfaceHtc;
            double htt = surfaceHtt;
            this->getForcedConvection(surface,elementID,htc,htt);
            this->getNaturalConvection(surface,elementID,htc,htt);

            this->elementHeatTransferCoefficient[elementID] = htc;

            Qx = Qy = Qz = 0.0;

            // Conduction
            if (surface.getThickness() > 0.0)
            {
                double Qi = 0.0;
                double Qj = 0.0;

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    // Heat flux is a density - see the line element loop above.
                    this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                    if (surface.getThickness() != 0.0)
                    {
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1]) / double(nInp);
                            B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1]) / double(nInp);
                        }
                    }
                }

                for (uint k=0;k<element.size();k++)
                {
                    uint nodeID = element.getNodeId(k);

                    Qi -= B[k][0] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                    Qj -= B[k][1] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                }

                RRMatrix R;
                RRVector t;
                this->pModel->getElement(elementID).findTransformationMatrix(this->pModel->getNodes(),R,t);

                Qx += R[0][0]*Qi + R[0][1]*Qj;
                Qy += R[1][0]*Qi + R[1][1]*Qj;
                Qz += R[2][0]*Qi + R[2][1]*Qj;
            }

            // Convection - Newton's law of cooling gives a flux density, so the
            // element area must not enter here.
            if (htc > 0.0)
            {
                RR3Vector normal;
                if (element.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]))
                {
                    double Qhe = htc * (htt - this->elementTemperature[elementID]);
                    Qx += normal[0] * Qhe;
                    Qy += normal[1] * Qhe;
                    Qz += normal[2] * Qhe;
                }
            }


            this->elementHeatFlux[elementID][0] = Qx;
            this->elementHeatFlux[elementID][1] = Qy;
            this->elementHeatFlux[elementID][2] = Qz;
        }
    }

    // Process volume elements.
    for (uint i=0;i<this->pModel->getNVolumes();i++)
    {
        RVolume &volume = this->pModel->getVolume(i);

        for (uint j=0;j<volume.size();j++)
        {
            uint elementID = volume.get(j);

            if (!this->computableElements[elementID])
            {
                continue;
            }

            const RElement &element = this->pModel->getElement(elementID);
            uint nInp = RElement::getNIntegrationPoints(element.getType());
            RRMatrix B(element.size(),3);

            Qx = Qy = Qz = 0.0;

            B.fill(0.0);

            // Conduction
            for (uint k=0;k<nInp;k++)
            {
                const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                const RRMatrix &dN = shapeFunc.getDN();
                RRMatrix J, Rt;
                // Heat flux is a density - see the line element loop above.
                this->pModel->getElement(elementID).findJacobian(this->pModel->getNodes(),k,J,Rt);

                for (uint m=0;m<dN.getNRows();m++)
                {
                    B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]) / double(nInp);
                    B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]) / double(nInp);
                    B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]) / double(nInp);
                }
            }

            for (uint m=0;m<element.size();m++)
            {
                uint nodeID = element.getNodeId(m);

                Qx -= B[m][0] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                Qy -= B[m][1] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
                Qz -= B[m][2] * this->elementConduction[elementID] * this->nodeTemperature[nodeID];
            }

            this->elementHeatFlux[elementID][0] = Qx;
            this->elementHeatFlux[elementID][1] = Qy;
            this->elementHeatFlux[elementID][2] = Qz;
        }
    }
}

void RSolverHeat::store()
{
    RLogger::info("Storing results\n");
    RLogger::indent();

    // Temperature
    uint temperaturePos = this->pModel->findVariable(R_VARIABLE_TEMPERATURE);
    if (temperaturePos == RConstants::eod)
    {
        temperaturePos = this->pModel->addVariable(R_VARIABLE_TEMPERATURE);
        this->pModel->getVariable(temperaturePos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->nodeTemperature),
                    RStatistics::findMaximumValue(this->nodeTemperature));
    }
    RVariable &temperature =  this->pModel->getVariable(temperaturePos);

    temperature.setApplyType(R_VARIABLE_APPLY_NODE);
    temperature.resize(1,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        temperature.setValue(0,i,this->nodeTemperature[i]);
    }

    // Heat flux
    uint heatFluxPos = this->pModel->findVariable(R_VARIABLE_HEAT_FLUX);
    if (heatFluxPos == RConstants::eod)
    {
        heatFluxPos = this->pModel->addVariable(R_VARIABLE_HEAT_FLUX);
        this->pModel->getVariable(heatFluxPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumMagnitude(this->elementHeatFlux),
                    RStatistics::findMaximumMagnitude(this->elementHeatFlux));
    }
    RVariable &heatFlux =  this->pModel->getVariable(heatFluxPos);

    heatFlux.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    heatFlux.resize(3,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        heatFlux.setValue(0,i,this->elementHeatFlux[i][0]);
        heatFlux.setValue(1,i,this->elementHeatFlux[i][1]);
        heatFlux.setValue(2,i,this->elementHeatFlux[i][2]);
    }

    // Heat transfer coefficient
    uint heatTransferCoefficientPos = this->pModel->findVariable(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
    if (heatTransferCoefficientPos == RConstants::eod)
    {
        heatTransferCoefficientPos = this->pModel->addVariable(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
        this->pModel->getVariable(heatTransferCoefficientPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementHeatTransferCoefficient),
                    RStatistics::findMaximumValue(this->elementHeatTransferCoefficient));
    }
    RVariable &heatTransferCoefficient = this->pModel->getVariable(heatTransferCoefficientPos);

    heatTransferCoefficient.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    heatTransferCoefficient.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        heatTransferCoefficient.setValue(0,i,this->elementHeatTransferCoefficient[i]);
    }

    RLogger::unindent();
}

void RSolverHeat::statistics()
{
    this->printStats(R_VARIABLE_TEMPERATURE);
    this->printStats(R_VARIABLE_HEAT_FLUX);
    this->printStats(R_VARIABLE_HEAT_TRANSFER_COEFFICIENT);
    this->processMonitoringPoints();
}

void RSolverHeat::assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp)
{
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    const RElement &element = this->pModel->getElement(elementID);

    RRMatrix Ae(element.size(),element.size());
    RRVector be(element.size());

    Ae.fill(0.0);
    be.fill(0.0);

    if (this->pModel->getTimeSolver().getEnabled())
    {
        for (unsigned m=0;m<element.size();m++)
        {
            be[m] = dt * fe[m];
            for (unsigned n=0;n<element.size();n++)
            {
                Ae[m][n] = Me[m][n] + alpha * dt * Ke[m][n];
                be[m] += (Me[m][n] - (1.0 - alpha) * dt * Ke[m][n]) * this->nodeTemperature[element.getNodeId(n)];
            }
        }
    }
    else
    {
        Ae = Ke;
        be = fe;
    }

    // Apply explicit boundary conditions.
    for (uint m=0;m<element.size();m++)
    {
        uint position;
        uint nodeID = element.getNodeId(m);
        if (!this->nodeBook.getValue(nodeID,position))
        {
            for (uint n=0;n<element.size();n++)
            {
                be[n] -= Ae[n][m] * this->nodeTemperature[nodeID];
            }
        }
    }


    // Assembly final matrix system
    for (uint m=0;m<element.size();m++)
    {
        uint mp;

        if (this->nodeBook.getValue(element.getNodeId(m),mp))
        {
            bp[mp] += be[m];
            for (uint n=0;n<element.size();n++)
            {
                uint np = 0;

                if (this->nodeBook.getValue(element.getNodeId(n),np))
                {
                    Ap.addValue(mp,np,Ae[m][n]);
                }
            }
        }
    }
}

bool RSolverHeat::getSimpleConvection(const RElementGroup &elementGroup, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE))
    {
        return false;
    }

    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE);
    uint cPos = 0;

    cPos = bc.findComponentPosition(R_VARIABLE_CONVECTION_COEFFICIENT);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_CONVECTION_COEFFICIENT).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE).toUtf8().constData());
    }
    htc = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    cPos = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_FLUID_TEMPERATURE).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_SIMPLE).toUtf8().constData());
    }
    htt = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    return true;
}

void RSolverHeat::checkConvectionInput(double value,
                                       RVariableType variableType,
                                       RBoundaryConditionType boundaryConditionType,
                                       const RElementGroup &elementGroup) const
{
    if (std::fabs(value) > RConstants::eps)
    {
        return;
    }

    throw RError(RError::Type::Application,R_ERROR_REF,
                 "Value of \'%s\' configured in \'%s\' boundary condition on entity \'%s\' is zero - the convection correlation can not be evaluated.",
                 RVariable::getName(variableType).toUtf8().constData(),
                 RBoundaryCondition::getName(boundaryConditionType).toUtf8().constData(),
                 elementGroup.getName().toUtf8().constData());
}

void RSolverHeat::reportForcedConvection(const RElementGroup &elementGroup)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED))
    {
        return;
    }

    double htc = 0.0;
    double htt = 0.0;
    double fastestFluid = 0.0;
    bool temperatureFromFluid = false;
    bool velocityFromFluid = false;
    bool applied = false;

    for (uint i=0;i<elementGroup.size();i++)
    {
        if (this->findFluidTemperature(elementGroup.get(i),htt))
        {
            temperatureFromFluid = true;
            break;
        }
    }
    for (uint i=0;i<elementGroup.size();i++)
    {
        double velocity = 0.0;
        if (this->findFluidVelocity(elementGroup.get(i),velocity))
        {
            velocityFromFluid = true;
            fastestFluid = std::max(fastestFluid,std::fabs(velocity));
        }
    }
    for (uint i=0;i<elementGroup.size();i++)
    {
        if (this->getForcedConvection(elementGroup,elementGroup.get(i),htc,htt))
        {
            applied = true;
            break;
        }
    }

    if (!applied)
    {
        if (velocityFromFluid && fastestFluid < RConstants::eps)
        {
            RLogger::warning("Forced convection on entity \'%s\' is disabled - the fluid heat result has the fluid at rest, so the correlation transfers no heat.\n",
                             elementGroup.getName().toUtf8().constData());
        }
        else
        {
            RLogger::warning("Forced convection on entity \'%s\' is ignored - no fluid heat result covers this surface and the condition carries no fluid temperature to fall back on.\n",
                             elementGroup.getName().toUtf8().constData());
        }
        return;
    }
    if (!temperatureFromFluid || !velocityFromFluid)
    {
        const char *configured = (!temperatureFromFluid && !velocityFromFluid) ? "fluid temperature and velocity"
                               : (!temperatureFromFluid)                       ? "fluid temperature"
                                                                               : "velocity";
        RLogger::info("Forced convection on entity \'%s\' uses the %s configured on the condition - no fluid heat result covers this surface.\n",
                      elementGroup.getName().toUtf8().constData(),
                      configured);
    }
}

bool RSolverHeat::getForcedConvection(const RElementGroup &elementGroup, uint elementId, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED))
    {
        return false;
    }

    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_FORCED);
    uint cPos = 0;

    // The fluid temperature is normally the one the fluid heat solver computed
    // on the fluid side of the wall. The value configured on the condition is a
    // fall-back for a surface no such result covers - one bordering no meshed
    // fluid, or a run where the fluid heat solver has not produced a result yet.
    if (!this->findFluidTemperature(elementId,htt))
    {
        uint fluidTemperaturePosition = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
        if (fluidTemperaturePosition == RConstants::eod)
        {
            // Conditions stored before the component existed carry no value to
            // fall back on, so there is nothing to exchange heat with.
            return false;
        }
        htt = bc.getComponent(fluidTemperaturePosition).get(this->pModel->getTimeSolver().getCurrentTime());
    }

    // Density
    cPos = bc.findComponentPosition(R_VARIABLE_DENSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DENSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double ro = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(ro,R_VARIABLE_DENSITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Dynamic viscosity
    cPos = bc.findComponentPosition(R_VARIABLE_DYNAMIC_VISCOSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DYNAMIC_VISCOSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double mu = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(mu,R_VARIABLE_DYNAMIC_VISCOSITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Heat capacity
    cPos = bc.findComponentPosition(R_VARIABLE_HEAT_CAPACITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HEAT_CAPACITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double c = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(c,R_VARIABLE_HEAT_CAPACITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Hydraulic diameter
    cPos = bc.findComponentPosition(R_VARIABLE_HYDRAULIC_DIAMETER);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HYDRAULIC_DIAMETER).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double d = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(d,R_VARIABLE_HYDRAULIC_DIAMETER,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Thermal conductivity
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_CONDUCTIVITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_CONDUCTIVITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
    }
    double k = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(k,R_VARIABLE_THERMAL_CONDUCTIVITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);

    // Mean velocity - like the fluid temperature it is normally the one the
    // fluid solver computed behind the wall, and the configured value is the
    // fall-back for a surface no such result covers.
    double v = 0.0;
    if (this->findFluidVelocity(elementId,v))
    {
        if (std::fabs(v) < RConstants::eps)
        {
            // The fluid solver has the fluid at rest. That is a result, not a
            // misconfiguration, so the wall is left unconvected for this pass
            // rather than falling back to a value the flow contradicts.
            return false;
        }
    }
    else
    {
        cPos = bc.findComponentPosition(R_VARIABLE_VELOCITY);
        if (cPos == RConstants::eod)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Failed to find \'%s\' component in \'%s\' boundary condition.",
                         RVariable::getName(R_VARIABLE_VELOCITY).toUtf8().constData(),
                         RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_FORCED).toUtf8().constData());
        }
        v = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
        this->checkConvectionInput(v,R_VARIABLE_VELOCITY,R_BOUNDARY_CONDITION_CONVECTION_FORCED,elementGroup);
    }

    RConvection convection;

    convection.setType(R_CONVECTION_FORCED_EXTERNAL);
    convection.setMaterial("Custom material",mu,ro,k,c,0.0);
    convection.setDiameter(d);
    convection.setVelocity(v);
    convection.setFluidTemp(htt);

    htc = convection.calculateHtc();

    return true;
}

bool RSolverHeat::getNaturalConvection(const RElementGroup &elementGroup, uint elementId, double &htc, double &htt)
{
    if (!elementGroup.hasBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_NATURAL))
    {
        return false;
    }

    RBoundaryCondition bc = elementGroup.getBoundaryCondition(R_BOUNDARY_CONDITION_CONVECTION_NATURAL);
    uint cPos = 0;

    // Density
    cPos = bc.findComponentPosition(R_VARIABLE_DENSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DENSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double ro = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    // Dynamic viscosity
    cPos = bc.findComponentPosition(R_VARIABLE_DYNAMIC_VISCOSITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_DYNAMIC_VISCOSITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double mu = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(mu,R_VARIABLE_DYNAMIC_VISCOSITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Fluid temperature
    cPos = bc.findComponentPosition(R_VARIABLE_FLUID_TEMPERATURE);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_FLUID_TEMPERATURE).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    htt = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    // Heat capacity
    cPos = bc.findComponentPosition(R_VARIABLE_HEAT_CAPACITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HEAT_CAPACITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double c = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(c,R_VARIABLE_HEAT_CAPACITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Hydraulic diameter
    cPos = bc.findComponentPosition(R_VARIABLE_HYDRAULIC_DIAMETER);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_HYDRAULIC_DIAMETER).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double d = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(d,R_VARIABLE_HYDRAULIC_DIAMETER,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Thermal conductivity
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_CONDUCTIVITY);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_CONDUCTIVITY).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double k = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());
    this->checkConvectionInput(k,R_VARIABLE_THERMAL_CONDUCTIVITY,R_BOUNDARY_CONDITION_CONVECTION_NATURAL,elementGroup);

    // Thermal expansion coefficient
    cPos = bc.findComponentPosition(R_VARIABLE_THERMAL_EXPANSION_COEFFICIENT);
    if (cPos == RConstants::eod)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to find \'%s\' component in \'%s\' boundary condition.",
                     RVariable::getName(R_VARIABLE_THERMAL_EXPANSION_COEFFICIENT).toUtf8().constData(),
                     RBoundaryCondition::getName(R_BOUNDARY_CONDITION_CONVECTION_NATURAL).toUtf8().constData());
    }
    double b = bc.getComponent(cPos).get(this->pModel->getTimeSolver().getCurrentTime());

    RConvection convection;

    convection.setType(R_CONVECTION_NATURAL_EXTERNAL_HORIZONTAL_PLATES);
    convection.setMaterial("Custom material",mu,ro,k,c,b);
    convection.setDiameter(d);
    convection.setSurfTemp(this->elementTemperature[elementId]);
    convection.setFluidTemp(htt);

    htc = convection.calculateHtc();

    return true;
}
