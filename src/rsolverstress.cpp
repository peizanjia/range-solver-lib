#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <map>

#include <QCryptographicHash>
#include <omp.h>

#include "rsolverstress.h"
#include "rmatrixsolver.h"
#include "rmodal_damping.h"
#include "reigenvaluesolver.h"
#include "rstructural_dynamics_integrator.h"
#include "rstructural_control_force_map.h"
#include "rstructural_cable_model_adapter.h"
#include "rstructural_static_truss_model.h"

#include <rml_problem_task_item.h>
#include <rml_structural_piezo_setup.h>

namespace
{

constexpr double axisymmetricFullAngle = 2.0*3.14159265358979323846;

bool finiteStressVector(const RRVector &vector)
{
    for (uint i=0;i<vector.getNRows();i++)
    {
        if (!std::isfinite(vector[i]))
        {
            return false;
        }
    }
    return true;
}

RDofType translationDof(uint component)
{
    R_ERROR_ASSERT(component < 3);
    return component == 0 ? RDofType::UX
         : component == 1 ? RDofType::UY
                          : RDofType::UZ;
}

bool findTranslationComponent(RDofType dof, uint &component)
{
    switch (dof)
    {
        case RDofType::UX: component = 0; return true;
        case RDofType::UY: component = 1; return true;
        case RDofType::UZ: component = 2; return true;
        default: return false;
    }
}

bool findRotationComponent(RDofType dof, uint &component)
{
    switch (dof)
    {
        case RDofType::RX: component = 0; return true;
        case RDofType::RY: component = 1; return true;
        case RDofType::RZ: component = 2; return true;
        default: return false;
    }
}

void addSymmetric(RRMatrix &matrix, uint row, uint column, double value)
{
    matrix[row][column] += value;
    if (row != column)
    {
        matrix[column][row] += value;
    }
}

//! Build a right-handed frame transformation. T maps local frame DOFs to
//! global DOFs and uses [u,v,w,rx,ry,rz] at each of the two nodes.
double buildBeamTransformation(const RLine &line,
                               const RElement &element,
                               const std::vector<RNode> &nodes,
                               RRMatrix &rotation,
                               RRMatrix &transformation)
{
    R_ERROR_ASSERT(element.getType() == R_ELEMENT_BEAM1 && element.size() == 2);
    const RR3Vector p0(nodes[element.getNodeId(0)].toVector());
    const RR3Vector p1(nodes[element.getNodeId(1)].toVector());
    RR3Vector ex;
    RR3Vector::subtract(p1,p0,ex);
    const double length = ex.normalize();
    if (length <= RConstants::eps)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Beam element %u has zero length.",element.getNodeId(0));
    }

    RR3Vector ey(line.getLocalYDirection());
    const double projection = RRVector::dot(ey,ex);
    for (uint i=0;i<3;i++)
    {
        ey[i] -= projection*ex[i];
    }
    if (ey.normalize() <= 1.0e-12)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Beam section local Y direction is parallel to element axis.");
    }
    RR3Vector ez;
    RR3Vector::cross(ex,ey,ez);
    ez.normalize();
    RR3Vector::cross(ez,ex,ey);
    ey.normalize();

    rotation.resize(3,3,0.0);
    for (uint i=0;i<3;i++)
    {
        rotation[i][0] = ex[i];
        rotation[i][1] = ey[i];
        rotation[i][2] = ez[i];
    }

    transformation.setIdentity(12);
    for (uint node=0;node<2;node++)
    {
        transformation.setBlock(rotation,6*node,6*node);
        transformation.setBlock(rotation,6*node+3,6*node+3);
    }
    return length;
}

//! Build a two-node 3D Euler-Bernoulli frame element in global coordinates.
void buildBeam1Element(const RLine &line,
                       const RElement &element,
                       const std::vector<RNode> &nodes,
                       double elasticity,
                       double poissonRatio,
                       double density,
                       double thermalExpansion,
                       double temperatureChange,
                       const RR3Vector &distributedForce,
                       bool needsMass,
                       bool lumpedMass,
                       RRMatrix &mass,
                       RRMatrix &stiffness,
                       RRVector &force)
{
    const double area = line.getCrossArea();
    const double iy = line.getSecondMomentY();
    const double iz = line.getSecondMomentZ();
    const double torsion = line.getTorsionalConstant();
    if (area <= 0.0 || iy <= 0.0 || iz <= 0.0 || torsion <= 0.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "BEAM1 requires positive A, Iy, Iz and J on line '%s'.",
                     line.getName().toUtf8().constData());
    }
    if (poissonRatio <= -1.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Beam Poisson ratio must be greater than -1.");
    }

    RRMatrix rotation,transformation;
    const double length = buildBeamTransformation(line,element,nodes,rotation,transformation);
    const double length2 = length*length;
    const double length3 = length2*length;
    const double shearModulus = elasticity/(2.0*(1.0+poissonRatio));

    RRMatrix localK(12,12,0.0);
    const double axial = elasticity*area/length;
    addSymmetric(localK,0,0,axial);
    addSymmetric(localK,6,6,axial);
    addSymmetric(localK,0,6,-axial);

    const double torsional = shearModulus*torsion/length;
    addSymmetric(localK,3,3,torsional);
    addSymmetric(localK,9,9,torsional);
    addSymmetric(localK,3,9,-torsional);

    const double bz = elasticity*iz/length3;
    localK[1][1] += 12.0*bz; localK[5][5] += 4.0*length2*bz;
    localK[7][7] += 12.0*bz; localK[11][11] += 4.0*length2*bz;
    addSymmetric(localK,1,5,6.0*length*bz);
    addSymmetric(localK,1,7,-12.0*bz);
    addSymmetric(localK,1,11,6.0*length*bz);
    addSymmetric(localK,5,7,-6.0*length*bz);
    addSymmetric(localK,5,11,2.0*length2*bz);
    addSymmetric(localK,7,11,-6.0*length*bz);

    const double by = elasticity*iy/length3;
    localK[2][2] += 12.0*by; localK[4][4] += 4.0*length2*by;
    localK[8][8] += 12.0*by; localK[10][10] += 4.0*length2*by;
    addSymmetric(localK,2,4,-6.0*length*by);
    addSymmetric(localK,2,8,-12.0*by);
    addSymmetric(localK,2,10,-6.0*length*by);
    addSymmetric(localK,4,8,6.0*length*by);
    addSymmetric(localK,4,10,2.0*length2*by);
    addSymmetric(localK,8,10,6.0*length*by);

    RRMatrix transformationT;
    transformationT.transpose(transformation);
    RRMatrix temporary;
    RRMatrix::mlt(transformation,localK,temporary);
    RRMatrix::mlt(temporary,transformationT,stiffness);

    RRVector localForce(12,0.0);
    localForce[0] -= elasticity*area*thermalExpansion*temperatureChange;
    localForce[6] += elasticity*area*thermalExpansion*temperatureChange;

    RRMatrix rotationT;
    rotationT.transpose(rotation);
    RR3Vector localDistributedForce;
    RRMatrix::mlt(rotationT,distributedForce,localDistributedForce);
    localForce[0] += localDistributedForce[0]*length/2.0;
    localForce[6] += localDistributedForce[0]*length/2.0;
    localForce[1] += localDistributedForce[1]*length/2.0;
    localForce[7] += localDistributedForce[1]*length/2.0;
    localForce[5] += localDistributedForce[1]*length2/12.0;
    localForce[11] -= localDistributedForce[1]*length2/12.0;
    localForce[2] += localDistributedForce[2]*length/2.0;
    localForce[8] += localDistributedForce[2]*length/2.0;
    localForce[4] -= localDistributedForce[2]*length2/12.0;
    localForce[10] += localDistributedForce[2]*length2/12.0;
    RRMatrix::mlt(transformation,localForce,force);

    mass.resize(12,12,0.0);
    if (!needsMass)
    {
        return;
    }

    RRMatrix localMass(12,12,0.0);
    if (lumpedMass)
    {
        const double nodalMass = density*area*length/2.0;
        for (uint node=0;node<2;node++)
        {
            localMass[6*node+0][6*node+0] = nodalMass;
            localMass[6*node+1][6*node+1] = nodalMass;
            localMass[6*node+2][6*node+2] = nodalMass;
            localMass[6*node+3][6*node+3] = density*(iy+iz)*length/2.0;
            localMass[6*node+4][6*node+4] = density*iy*length/2.0;
            localMass[6*node+5][6*node+5] = density*iz*length/2.0;
        }
    }
    else
    {
        const double axialMass = density*area*length/6.0;
        localMass[0][0] += 2.0*axialMass; localMass[6][6] += 2.0*axialMass;
        addSymmetric(localMass,0,6,axialMass);

        const double torsionalMass = density*(iy+iz)*length/6.0;
        localMass[3][3] += 2.0*torsionalMass; localMass[9][9] += 2.0*torsionalMass;
        addSymmetric(localMass,3,9,torsionalMass);

        const double bendingMass = density*area*length/420.0;
        const uint v[4] = {1,5,7,11};
        const double mv[4][4] = {
            {156.0,22.0*length,54.0,-13.0*length},
            {22.0*length,4.0*length2,13.0*length,-3.0*length2},
            {54.0,13.0*length,156.0,-22.0*length},
            {-13.0*length,-3.0*length2,-22.0*length,4.0*length2}};
        const uint w[4] = {2,4,8,10};
        const double mw[4][4] = {
            {156.0,-22.0*length,54.0,13.0*length},
            {-22.0*length,4.0*length2,-13.0*length,-3.0*length2},
            {54.0,-13.0*length,156.0,22.0*length},
            {13.0*length,-3.0*length2,22.0*length,4.0*length2}};
        for (uint i=0;i<4;i++)
        {
            for (uint j=0;j<4;j++)
            {
                localMass[v[i]][v[j]] += bendingMass*mv[i][j];
                localMass[w[i]][w[j]] += bendingMass*mw[i][j];
            }
        }
    }
    RRMatrix::mlt(transformation,localMass,temporary);
    RRMatrix::mlt(temporary,transformationT,mass);
}

//! Build either the two-node frame or the three-node frame contract. BEAM2
//! uses its midside node as a real frame node and is assembled from two exact
//! Euler-Bernoulli spans in the stable [end 0, end 1, middle] node ordering.
void buildBeamElement(const RLine &line,
                      const RElement &element,
                      const std::vector<RNode> &nodes,
                      double elasticity,
                      double poissonRatio,
                      double density,
                      double thermalExpansion,
                      double temperatureChange,
                      const RR3Vector &distributedForce,
                      bool needsMass,
                      bool lumpedMass,
                      RRMatrix &mass,
                      RRMatrix &stiffness,
                      RRVector &force)
{
    if (element.getType() == R_ELEMENT_BEAM1)
    {
        buildBeam1Element(line,element,nodes,elasticity,poissonRatio,density,
                          thermalExpansion,temperatureChange,distributedForce,
                          needsMass,lumpedMass,mass,stiffness,force);
        return;
    }
    if (element.getType() != R_ELEMENT_BEAM2 || element.size() != 3)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Unsupported beam element contract.");
    }

    mass.resize(18,18,0.0);
    stiffness.resize(18,18,0.0);
    force.resize(18,0.0);
    const uint parentNodes[2][2] = {{0,2},{2,1}};
    for (uint span=0;span<2;span++)
    {
        RElement child(R_ELEMENT_BEAM1);
        child.setNodeId(0,element.getNodeId(parentNodes[span][0]));
        child.setNodeId(1,element.getNodeId(parentNodes[span][1]));
        RRMatrix childMass,childStiffness;
        RRVector childForce;
        buildBeam1Element(line,child,nodes,elasticity,poissonRatio,density,
                          thermalExpansion,temperatureChange,distributedForce,
                          needsMass,lumpedMass,childMass,childStiffness,childForce);
        for (uint a=0;a<2;a++)
        {
            const uint parentA = parentNodes[span][a];
            for (uint da=0;da<6;da++)
            {
                const uint row = 6*parentA+da;
                const uint childRow = 6*a+da;
                force[row] += childForce[childRow];
                for (uint b=0;b<2;b++)
                {
                    const uint parentB = parentNodes[span][b];
                    for (uint db=0;db<6;db++)
                    {
                        const uint column = 6*parentB+db;
                        const uint childColumn = 6*b+db;
                        stiffness[row][column] += childStiffness[childRow][childColumn];
                        mass[row][column] += childMass[childRow][childColumn];
                    }
                }
            }
        }
    }
}

void lumpMassMatrix(RRMatrix &matrix)
{
    R_ERROR_ASSERT(matrix.getNRows() == matrix.getNColumns());
    RRVector rowSums(matrix.getNRows(),0.0);
    for (uint i=0;i<matrix.getNRows();i++)
    {
        for (uint j=0;j<matrix.getNColumns();j++)
        {
            rowSums[i] += matrix[i][j];
        }
    }
    matrix.fill(0.0);
    for (uint i=0;i<matrix.getNRows();i++)
    {
        matrix[i][i] = rowSums[i];
    }
}

void buildSurfaceConstitutiveMatrix(RElementFormulation formulation,
                                    double elasticityModulus,
                                    double poissonRatio,
                                    RRMatrix &matrix)
{
    matrix.resize(3,3,0.0);
    matrix.fill(0.0);

    if (formulation == RElementFormulation::PlaneStress)
    {
        const double scale = elasticityModulus/(1.0-poissonRatio*poissonRatio);
        matrix[0][0] = scale;
        matrix[0][1] = scale*poissonRatio;
        matrix[1][0] = scale*poissonRatio;
        matrix[1][1] = scale;
        matrix[2][2] = elasticityModulus/(2.0*(1.0+poissonRatio));
        return;
    }

    if (formulation == RElementFormulation::PlaneStrain)
    {
        matrix[0][0] = 1.0-poissonRatio;
        matrix[0][1] = poissonRatio;
        matrix[1][0] = poissonRatio;
        matrix[1][1] = 1.0-poissonRatio;
        matrix[2][2] = (1.0-2.0*poissonRatio)/2.0;
        matrix *= elasticityModulus/((1.0+poissonRatio)*(1.0-2.0*poissonRatio));
        return;
    }

    if (formulation == RElementFormulation::Axisymmetric)
    {
        matrix.resize(4,4,0.0);
        matrix.fill(0.0);
        matrix[0][0] = 1.0-poissonRatio;
        matrix[0][1] = poissonRatio;
        matrix[0][2] = poissonRatio;
        matrix[1][0] = poissonRatio;
        matrix[1][1] = 1.0-poissonRatio;
        matrix[1][2] = poissonRatio;
        matrix[2][0] = poissonRatio;
        matrix[2][1] = poissonRatio;
        matrix[2][2] = 1.0-poissonRatio;
        matrix[3][3] = (1.0-2.0*poissonRatio)/2.0;
        matrix *= elasticityModulus/((1.0+poissonRatio)*(1.0-2.0*poissonRatio));
        return;
    }

    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                 "Structural surface formulation '%s' is not supported by the stress solver.",
                 RElementDescriptor::formulationToId(formulation).toUtf8().constData());
}

void validateSurfaceFormulation(const RSurface &surface)
{
    if (surface.getStructuralFormulation() != RElementFormulation::PlaneStress &&
        surface.getStructuralFormulation() != RElementFormulation::PlaneStrain &&
        surface.getStructuralFormulation() != RElementFormulation::Axisymmetric)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Unsupported structural formulation on surface '%s'.",
                     surface.getName().toUtf8().constData());
    }
}

//! Build the radius/axis derivatives and the local-to-global displacement map
//! for the first axisymmetric contract. The meridian is the global XY plane:
//! X is radius, Y is the revolution axis and Z is circumferential.
double findAxisymmetricKinematics(const RElement &element,
                                  const std::vector<RNode> &nodes,
                                  const RElementShapeFunction &shapeFunction,
                                  RRMatrix &derivatives,
                                  RRMatrix &transformation,
                                  double &radius)
{
    const RRVector &N = shapeFunction.getN();
    const RRMatrix &dN = shapeFunction.getDN();
    double dxDksi = 0.0;
    double dxDeta = 0.0;
    double dyDksi = 0.0;
    double dyDeta = 0.0;
    radius = 0.0;
    double coordinateScale = 1.0;

    for (uint i=0;i<element.size();i++)
    {
        const RNode &node = nodes[element.getNodeId(i)];
        coordinateScale = std::max(coordinateScale,
                                   std::max(std::fabs(node.getX()),std::fabs(node.getY())));
        if (std::fabs(node.getZ()) > 1.0e-10*coordinateScale)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Axisymmetric elements must lie in the global XY meridian plane "
                         "(X = radius, Y = axis); node %u has Z = %g.",
                         element.getNodeId(i),node.getZ());
        }
        if (node.getX() < -1.0e-10*coordinateScale)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Axisymmetric radius must be non-negative; node %u has X = %g.",
                         element.getNodeId(i),node.getX());
        }
        radius += N[i]*node.getX();
        dxDksi += dN[i][0]*node.getX();
        dxDeta += dN[i][1]*node.getX();
        dyDksi += dN[i][0]*node.getY();
        dyDeta += dN[i][1]*node.getY();
    }

    const double signedDeterminant = dxDksi*dyDeta-dxDeta*dyDksi;
    const double determinant = std::fabs(signedDeterminant);
    if (determinant <= 1.0e-14*coordinateScale*coordinateScale)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Axisymmetric element has a singular meridian Jacobian.");
    }
    if (radius <= 1.0e-12*coordinateScale)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Axisymmetric integration point lies on or beyond the revolution axis (radius = %g).",
                     radius);
    }

    const double ksiX = dyDeta/signedDeterminant;
    const double ksiY = -dxDeta/signedDeterminant;
    const double etaX = -dyDksi/signedDeterminant;
    const double etaY = dxDksi/signedDeterminant;
    derivatives.resize(element.size(),2,0.0);
    for (uint i=0;i<element.size();i++)
    {
        derivatives[i][0] = dN[i][0]*ksiX+dN[i][1]*etaX;
        derivatives[i][1] = dN[i][0]*ksiY+dN[i][1]*etaY;
    }

    transformation.resize(3*element.size(),2*element.size(),0.0);
    transformation.fill(0.0);
    for (uint i=0;i<element.size();i++)
    {
        transformation[3*i+0][2*i+0] = 1.0;
        transformation[3*i+1][2*i+1] = 1.0;
    }
    return determinant;
}

double findAxisymmetricMeasure(const RSurface &surface,
                               const std::vector<RNode> &nodes,
                               const std::vector<RElement> &elements)
{
    double measure = 0.0;
    RRMatrix derivatives;
    RRMatrix transformation;
    for (uint i=0;i<surface.size();i++)
    {
        const RElement &element = elements[surface.get(i)];
        const uint nIntegrationPoints = RElement::getNIntegrationPoints(element.getType());
        for (uint k=0;k<nIntegrationPoints;k++)
        {
            const RElementShapeFunction &shapeFunction = RElement::getShapeFunction(element.getType(),k);
            double radius = 0.0;
            const double determinant = findAxisymmetricKinematics(element,nodes,shapeFunction,
                                                                   derivatives,transformation,radius);
            measure += axisymmetricFullAngle*radius*determinant*shapeFunction.getW();
        }
    }
    return measure;
}

}

struct RSolverStress::DynamicsStepSnapshot
{
    bool hadIntegrator = false;
    RStructuralDynamicsCheckpoint integratorCheckpoint;
    RRVector previousForce;
    RRVector currentForce;
    RRVector initialForce;
    bool preparingInitialForce = false;
    RSparseMatrix massReference;
    RSparseMatrix stiffnessReference;
    RSparseMatrix motionStiffness;
    RSparseMatrix motionMass;
    RSparseMatrix motionStiffnessReference;
    RSparseMatrix motionMassReference;
    std::vector<RStructuralPrescribedMotion> prescribedMotions;
    std::vector<int> motionIndexByNodeDof;
    QString motionDefinitionReference;
    RRMatrix constraintNullspaceReference;
    RRVector constraintParticularReference;
    std::vector<double> boundaryReference;
    bool boundaryReferenceReady = false;
    double previousTime = 0.0;
    double currentTime = 0.0;
    RSparseMatrix supportStiffness;
    RSparseMatrix supportMass;
    RSparseMatrix supportStiffnessReference;
    RSparseMatrix supportMassReference;
    RRVector supportForce;
    RRVector initialSupportForce;
    std::vector<RStructuralDofRef> heldDofs;
    std::vector<int> heldIndexByNodeDof;
    std::vector<bool> pointSupportMask;
    std::vector<RStructuralSupportSample> acceptedSupportSamples;
    std::vector<RStructuralSupportSample> trialSupportSamples;
    bool hadControlMap = false;
    QString cachedControlFingerprint;
    QString cachedControlCatalogFingerprint;
    std::vector<RStructuralControlChannelSample> controlLedgerSamples;
    RStructuralControlConditioner::State acceptedConditionerState;
    RStructuralControlConditioner::State candidateConditionerState;
    std::vector<RStructuralControlCommand> activeControlCommands;
    RRVector activeControlForce;
    QByteArray acceptedControllerState;
    QByteArray candidateControllerState;
    RStructuralControlStepRecord trialControlRecord;
    std::vector<ControlSensorMap> controlSensors;
    RStructuralCableDynamicsResult structuralCableDynamicsResult;
    bool structuralCableDynamicsResultReady = false;
    bool hasCableAcceptedTime = false;
    uint cableAcceptedTimeStep = 0;
    double cableAcceptedComputedTime = 0.0;
};

RSolverStress::RSolverStress(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData, bool modalAnalysis)
    : RSolverGeneric(pModel,modelFileName,convergenceFileName,sharedData)
    , nLinearConstraintEquations(0)
{
    this->problemType = modalAnalysis ? R_PROBLEM_STRESS_MODAL : R_PROBLEM_STRESS;
    this->structuralCableModeAtConstruction =
        this->pModel->getProblemSetup().getStructuralCableSetup().enabled;
    this->structuralCableDynamicsModeAtConstruction =
        this->structuralCableModeAtConstruction &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    this->structuralPiezoModeAtConstruction =
        this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled;
    if (this->structuralCableModeAtConstruction || this->structuralPiezoModeAtConstruction)
    {
        if (modalAnalysis)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Enabled structural cable/piezo mode requires a non-modal static stress solver.");
        }
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
    }
}

bool sameControlTimeSolver(const RTimeSolver &left,const RTimeSolver &right)
{
    return left.getEnabled()==right.getEnabled() &&
        left.getTimeMarchApproximation()==right.getTimeMarchApproximation() &&
        left.getTimes()==right.getTimes() &&
        left.getInputNTimeSteps()==right.getInputNTimeSteps() &&
        left.getInputStartTime()==right.getInputStartTime() &&
        left.getInputTimeStepSize()==right.getInputTimeStepSize() &&
        left.getCurrentTimeStep()==right.getCurrentTimeStep() &&
        left.getComputedTime()==right.getComputedTime() &&
        left.getOutputFrequency()==right.getOutputFrequency();
}

bool sameControlSetup(const RStructuralControlSetup &left,
                      const RStructuralControlSetup &right)
{
    if (left.enabled!=right.enabled ||
        left.controlProtocolVersion!=right.controlProtocolVersion ||
        left.controllerId!=right.controllerId ||
        left.deterministicSeed!=right.deterministicSeed ||
        left.channelConditioning.size()!=right.channelConditioning.size())
        return false;
    for (std::size_t i=0;i<left.channelConditioning.size();i++)
    {
        const auto &a = left.channelConditioning[i];
        const auto &b = right.channelConditioning[i];
        if (a.actuatorId!=b.actuatorId || a.delaySteps!=b.delaySteps ||
            a.holdSteps!=b.holdSteps || a.noiseAmplitude!=b.noiseAmplitude)
            return false;
    }
    return true;
}

namespace
{

void appendProblemTaskLeaves(const RProblemTaskItem &item,std::vector<RProblemType> &leaves)
{
    if (item.getProblemType() != R_PROBLEM_NONE)
    {
        leaves.push_back(item.getProblemType());
        return;
    }
    for (uint child=0;child<item.getNChildren();child++)
    {
        appendProblemTaskLeaves(item.getChild(child),leaves);
    }
}

//! Prestress is accepted only from an executed Heat -> Stress sequence that
//! precedes the modal task in the current task tree. A shared-data entry left
//! by a prior solver run must not silently activate geometric stiffness.
bool hasCurrentThermalStressBeforeModal(const RProblemTaskItem &taskTree)
{
    std::vector<RProblemType> leaves;
    appendProblemTaskLeaves(taskTree,leaves);
    bool heatSeen = false;
    bool stressAfterHeatSeen = false;
    for (RProblemType type : leaves)
    {
        if (type == R_PROBLEM_STRESS_MODAL)
        {
            return heatSeen && stressAfterHeatSeen;
        }
        if (type == R_PROBLEM_HEAT && !stressAfterHeatSeen)
        {
            heatSeen = true;
        }
        else if (type == R_PROBLEM_STRESS && heatSeen)
        {
            stressAfterHeatSeen = true;
        }
    }
    return false;
}

void addTrussGeometricStiffness(const RElement &element,
                                const std::vector<RNode> &nodes,
                                double axialForce,
                                RRMatrix &stiffness)
{
    if (element.getType() != R_ELEMENT_TRUSS1 || element.size() != 2 ||
        stiffness.getNRows() != 6 || stiffness.getNColumns() != 6)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Thermal geometric stiffness supports only TRUSS1 elements.");
    }

    RR3Vector direction;
    RR3Vector::subtract(nodes[element.getNodeId(1)].toVector(),
                        nodes[element.getNodeId(0)].toVector(),direction);
    const double length = direction.normalize();
    if (!std::isfinite(length) || length <= RConstants::eps || !std::isfinite(axialForce))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "TRUSS1 geometric stiffness requires finite axial force and nonzero length.");
    }

    // The consistent initial-stress block is N/L times the transverse
    // projector, with positive N denoting tension.
    for (uint rowNode=0;rowNode<2;rowNode++)
    {
        const double rowSign = rowNode == 0 ? 1.0 : -1.0;
        for (uint columnNode=0;columnNode<2;columnNode++)
        {
            const double columnSign = columnNode == 0 ? 1.0 : -1.0;
            for (uint row=0;row<3;row++)
            {
                for (uint column=0;column<3;column++)
                {
                    const double projector = (row == column ? 1.0 : 0.0)
                                           - direction[row]*direction[column];
                    stiffness[3*rowNode+row][3*columnNode+column] +=
                        axialForce/length*rowSign*columnSign*projector;
                }
            }
        }
    }
}

void addBeam1GeometricStiffness(const RLine &line,
                                const RElement &element,
                                const std::vector<RNode> &nodes,
                                double axialForce,
                                RRMatrix &stiffness)
{
    if (element.getType() != R_ELEMENT_BEAM1 || element.size() != 2 ||
        stiffness.getNRows() != 12 || stiffness.getNColumns() != 12)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Thermal geometric stiffness supports only BEAM1 elements.");
    }

    RRMatrix rotation,transformation;
    const double length = buildBeamTransformation(line,element,nodes,rotation,transformation);
    if (!std::isfinite(axialForce) || !std::isfinite(length) || length <= RConstants::eps)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "BEAM1 geometric stiffness requires finite axial force and nonzero length.");
    }

    // Integral of N_i' N_j' for cubic Hermite bending interpolation. Tension
    // adds stiffness; compression subtracts it. The local z-rotation used by
    // the w-plane beam DOFs is opposite to dw/dx, so that block changes sign
    // on its rotational rows and columns.
    const double length2 = length*length;
    const double geometric[4][4] = {
        {36.0,3.0*length,-36.0,3.0*length},
        {3.0*length,4.0*length2,-3.0*length,-length2},
        {-36.0,-3.0*length,36.0,-3.0*length},
        {3.0*length,-length2,-3.0*length,4.0*length2}};
    const uint vDofs[4] = {1,5,7,11};
    const uint wDofs[4] = {2,4,8,10};
    const double vSigns[4] = {1.0,1.0,1.0,1.0};
    const double wSigns[4] = {1.0,-1.0,1.0,-1.0};
    RRMatrix local(12,12,0.0);
    const double factor = axialForce/(30.0*length);
    for (uint row=0;row<4;row++)
    {
        for (uint column=0;column<4;column++)
        {
            local[vDofs[row]][vDofs[column]] += factor*geometric[row][column]
                                               * vSigns[row]*vSigns[column];
            local[wDofs[row]][wDofs[column]] += factor*geometric[row][column]
                                               * wSigns[row]*wSigns[column];
        }
    }

    RRMatrix transformationT,temporary,global;
    transformationT.transpose(transformation);
    RRMatrix::mlt(transformation,local,temporary);
    RRMatrix::mlt(temporary,transformationT,global,true);
    for (uint row=0;row<12;row++)
    {
        for (uint column=0;column<12;column++)
        {
            stiffness[row][column] += global[row][column];
        }
    }
}

}

// RMaterialProperty uses piecewise-linear interpolation with endpoint clamping.
// The mean alpha over the actual temperature interval makes alpha*DeltaT equal
// the integrated free thermal strain, including when alpha varies with T.
double meanThermalExpansion(const RMaterialProperty &property,
                            double referenceKelvin,double currentKelvin)
{
    if (!std::isfinite(referenceKelvin) || referenceKelvin < 0.0 ||
        !std::isfinite(currentKelvin) || currentKelvin < 0.0 || property.empty())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Thermal expansion requires finite absolute-kelvin temperatures and a nonempty alpha table.");
    }
    if (referenceKelvin == currentKelvin)
        return property.get(currentKelvin);

    const double low = std::min(referenceKelvin,currentKelvin);
    const double high = std::max(referenceKelvin,currentKelvin);
    double previous = low;
    double integral = 0.0;
    for (uint point=0;point<property.size();point++)
    {
        const double knot = property.getKey(point);
        if (!std::isfinite(knot))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Thermal expansion table has a non-finite temperature key.");
        if (knot <= previous || knot >= high) continue;
        const double left = property.get(previous);
        const double right = property.get(knot);
        if (!std::isfinite(left) || !std::isfinite(right))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Thermal expansion table has a non-finite coefficient.");
        integral += 0.5*(left+right)*(knot-previous);
        previous = knot;
    }
    const double left = property.get(previous);
    const double right = property.get(high);
    if (!std::isfinite(left) || !std::isfinite(right))
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Thermal expansion table has a non-finite coefficient.");
    integral += 0.5*(left+right)*(high-previous);
    const double mean = integral/(high-low);
    if (!std::isfinite(mean))
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Integrated thermal expansion is non-finite.");
    return mean;
}

RSolverStress::~RSolverStress()
{

}

void RSolverStress::beginStepTransaction()
{
    if (!this->pendingControlObserverFailure.isEmpty())
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Control step observer failed after a prior committed result: %s",
                     this->pendingControlObserverFailure.toUtf8().constData());
    std::unique_ptr<DynamicsStepSnapshot> snapshot(new DynamicsStepSnapshot);
    snapshot->hadIntegrator = bool(this->dynamicsIntegrator);
    if (snapshot->hadIntegrator)
    {
        snapshot->integratorCheckpoint = this->dynamicsIntegrator->checkpoint();
    }
    snapshot->previousForce = this->previousDynamicForce;
    snapshot->currentForce = this->currentDynamicForce;
    snapshot->initialForce = this->initialDynamicForce;
    snapshot->preparingInitialForce = this->preparingInitialForce;
    snapshot->massReference = this->dynamicMassReference;
    snapshot->stiffnessReference = this->dynamicStiffnessReference;
    snapshot->motionStiffness = this->motionStiffness;
    snapshot->motionMass = this->motionMass;
    snapshot->motionStiffnessReference = this->dynamicMotionStiffnessReference;
    snapshot->motionMassReference = this->dynamicMotionMassReference;
    snapshot->prescribedMotions = this->prescribedMotions;
    snapshot->motionIndexByNodeDof = this->motionIndexByNodeDof;
    snapshot->motionDefinitionReference = this->dynamicMotionDefinitionReference;
    snapshot->constraintNullspaceReference = this->dynamicConstraintNullspaceReference;
    snapshot->constraintParticularReference = this->dynamicConstraintParticularReference;
    snapshot->boundaryReference = this->dynamicBoundaryReference;
    snapshot->boundaryReferenceReady = this->dynamicBoundaryReferenceReady;
    snapshot->previousTime = this->previousDynamicTime;
    snapshot->currentTime = this->currentDynamicTime;
    snapshot->supportStiffness = this->supportStiffness;
    snapshot->supportMass = this->supportMass;
    snapshot->supportStiffnessReference = this->supportStiffnessReference;
    snapshot->supportMassReference = this->supportMassReference;
    snapshot->supportForce = this->supportForce;
    snapshot->initialSupportForce = this->initialSupportForce;
    snapshot->heldDofs = this->heldDofs;
    snapshot->heldIndexByNodeDof = this->heldIndexByNodeDof;
    snapshot->pointSupportMask = this->pointSupportMask;
    snapshot->acceptedSupportSamples = this->acceptedSupportSamples;
    snapshot->trialSupportSamples = this->trialSupportSamples;
    snapshot->hadControlMap = bool(this->controlForceMap);
    snapshot->cachedControlFingerprint = this->cachedControlFingerprint;
    snapshot->cachedControlCatalogFingerprint = this->cachedControlCatalogFingerprint;
    if (this->controlLedger)
        snapshot->controlLedgerSamples = this->controlLedger->accepted();
    snapshot->acceptedConditionerState = this->acceptedConditionerState;
    snapshot->candidateConditionerState = this->candidateConditionerState;
    snapshot->activeControlCommands = this->activeControlCommands;
    snapshot->activeControlForce = this->activeControlForce;
    snapshot->acceptedControllerState = this->acceptedControllerState;
    snapshot->candidateControllerState = this->candidateControllerState;
    snapshot->trialControlRecord = this->trialControlRecord;
    snapshot->controlSensors = this->controlSensors;
    snapshot->structuralCableDynamicsResult = this->structuralCableDynamicsResult;
    snapshot->structuralCableDynamicsResultReady = this->structuralCableDynamicsResultReady;
    if (this->isStructuralCableDynamicsMode())
    {
        const RTimeSolver &time = this->pModel->getTimeSolver();
        snapshot->hasCableAcceptedTime = true;
        snapshot->cableAcceptedTimeStep = time.getCurrentTimeStep() > 0
            ? time.getCurrentTimeStep()-1 : 0;
        snapshot->cableAcceptedComputedTime = time.getPreviousTime();
    }
    this->dynamicsStepSnapshot = std::move(snapshot);
}

void RSolverStress::rollbackStepTransaction()
{
    if (!this->dynamicsStepSnapshot)
    {
        return;
    }
    const DynamicsStepSnapshot &snapshot = *this->dynamicsStepSnapshot;
    if (snapshot.hadIntegrator)
    {
        this->dynamicsIntegrator->restore(snapshot.integratorCheckpoint);
    }
    else
    {
        this->dynamicsIntegrator.reset();
    }
    this->previousDynamicForce = snapshot.previousForce;
    this->currentDynamicForce = snapshot.currentForce;
    this->initialDynamicForce = snapshot.initialForce;
    this->preparingInitialForce = snapshot.preparingInitialForce;
    this->dynamicMassReference = snapshot.massReference;
    this->dynamicStiffnessReference = snapshot.stiffnessReference;
    this->motionStiffness = snapshot.motionStiffness;
    this->motionMass = snapshot.motionMass;
    this->dynamicMotionStiffnessReference = snapshot.motionStiffnessReference;
    this->dynamicMotionMassReference = snapshot.motionMassReference;
    this->prescribedMotions = snapshot.prescribedMotions;
    this->motionIndexByNodeDof = snapshot.motionIndexByNodeDof;
    this->dynamicMotionDefinitionReference = snapshot.motionDefinitionReference;
    this->dynamicConstraintNullspaceReference = snapshot.constraintNullspaceReference;
    this->dynamicConstraintParticularReference = snapshot.constraintParticularReference;
    this->dynamicBoundaryReference = snapshot.boundaryReference;
    this->dynamicBoundaryReferenceReady = snapshot.boundaryReferenceReady;
    this->previousDynamicTime = snapshot.previousTime;
    this->currentDynamicTime = snapshot.currentTime;
    this->supportStiffness = snapshot.supportStiffness;
    this->supportMass = snapshot.supportMass;
    this->supportStiffnessReference = snapshot.supportStiffnessReference;
    this->supportMassReference = snapshot.supportMassReference;
    this->supportForce = snapshot.supportForce;
    this->initialSupportForce = snapshot.initialSupportForce;
    this->heldDofs = snapshot.heldDofs;
    this->heldIndexByNodeDof = snapshot.heldIndexByNodeDof;
    this->pointSupportMask = snapshot.pointSupportMask;
    this->acceptedSupportSamples = snapshot.acceptedSupportSamples;
    this->trialSupportSamples = snapshot.trialSupportSamples;
    if (!snapshot.hadControlMap)
    {
        this->controlForceMap.reset();
        this->controlLedger.reset();
        this->controlConditioner.reset();
    }
    else if (this->controlLedger)
    {
        this->controlLedger->restore(snapshot.controlLedgerSamples);
    }
    this->activeControlCommands = snapshot.activeControlCommands;
    this->acceptedConditionerState = snapshot.acceptedConditionerState;
    this->candidateConditionerState = snapshot.candidateConditionerState;
    this->cachedControlFingerprint = snapshot.cachedControlFingerprint;
    this->cachedControlCatalogFingerprint = snapshot.cachedControlCatalogFingerprint;
    this->activeControlForce = snapshot.activeControlForce;
    this->acceptedControllerState = snapshot.acceptedControllerState;
    this->candidateControllerState = snapshot.candidateControllerState;
    this->trialControlRecord = snapshot.trialControlRecord;
    this->controlSensors = snapshot.controlSensors;
    this->structuralCableDynamicsResult = snapshot.structuralCableDynamicsResult;
    this->structuralCableDynamicsResultReady = snapshot.structuralCableDynamicsResultReady;
    if (snapshot.hasCableAcceptedTime && this->isStructuralCableDynamicsMode())
    {
        RTimeSolver &time = this->pModel->getTimeSolver();
        time.setCurrentTimeStep(snapshot.cableAcceptedTimeStep);
        time.setComputedTime(snapshot.cableAcceptedComputedTime);
    }
    this->dynamicsStepSnapshot.reset();
}

void RSolverStress::commitStepTransaction()
{
    if (this->controlLedger)
    {
        this->controlLedger->accept();
        this->acceptedControllerState = std::move(this->candidateControllerState);
        this->acceptedConditionerState = std::move(this->candidateConditionerState);
    }
    this->acceptedSupportSamples = std::move(this->trialSupportSamples);
    this->trialSupportSamples.clear();
    this->dynamicsStepSnapshot.reset();
    if (this->supportResponseObserver)
    {
        try
        {
            // The result file and in-memory step are already committed.
            // Observers are notifications, never a reason to roll back a
            // published checkpoint after an external side effect.
            this->supportResponseObserver(this->acceptedSupportSamples);
        }
        catch (const std::exception &error)
        {
            RLogger::warning("Support-response observer failed after commit: %s\n",error.what());
        }
        catch (...)
        {
            RLogger::warning("Support-response observer failed after commit.\n");
        }
    }
    if (this->controlStepObserver && this->controlLedger)
    {
        try { this->controlStepObserver(this->trialControlRecord); }
        catch (const std::exception &error)
        {
            this->pendingControlObserverFailure = QString::fromUtf8(error.what());
            RLogger::warning("Control observer failed after committed result: %s\n",error.what());
        }
        catch (...)
        {
            this->pendingControlObserverFailure = "Unknown control observer failure";
            RLogger::warning("Control observer failed after committed result.\n");
        }
    }
}

void RSolverStress::setControlCallback(
    const QString &id, const QString &revision,
    std::function<RStructuralControlDecision(const RStructuralControlObservation &)> callback)
{
    if (this->isSpecializedStaticTrussMode() && callback)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control callbacks are unsupported in static cable/piezo mode.");
    }
    this->controllerId = id;
    this->controllerRevision = revision;
    this->controlCallback = std::move(callback);
}

void RSolverStress::setControlStepObserver(
    std::function<void(const RStructuralControlStepRecord &)> observer)
{
    if (this->isSpecializedStaticTrussMode() && observer)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control observers are unsupported in static cable/piezo mode.");
    }
    this->controlStepObserver = std::move(observer);
}

std::vector<RStructuralControlChannelSample> RSolverStress::getAcceptedControlSamples() const
{
    return this->controlLedger ? this->controlLedger->accepted()
                               : std::vector<RStructuralControlChannelSample>();
}

QString RSolverStress::getPendingControlObserverFailure() const
{
    return this->pendingControlObserverFailure;
}

void RSolverStress::setSupportResponseObserver(
    std::function<void(const std::vector<RStructuralSupportSample> &)> observer)
{
    if (this->isSpecializedStaticTrussMode() && observer)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural support-response observers are unsupported in static cable/piezo mode.");
    }
    this->supportResponseObserver = std::move(observer);
}

std::vector<RStructuralSupportSample> RSolverStress::getAcceptedSupportSamples() const
{
    return this->acceptedSupportSamples;
}

bool RSolverStress::hasConverged() const
{
    if (this->isStructuralCableMode())
    {
        if (this->isStructuralCableDynamicsMode())
        {
            return this->structuralCableDynamicsResultReady;
        }
        return this->structuralCableResultReady;
    }
    if (this->isStructuralPiezoMode())
    {
        return this->structuralPiezoResultReady;
    }
    return true;
}

void RSolverStress::validateRunConfiguration() const
{
    const bool cableEnabled = this->pModel->getProblemSetup().getStructuralCableSetup().enabled;
    const bool piezoEnabled = this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled;
    const bool cableDynamicsEnabled = cableEnabled &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    if (cableEnabled != this->structuralCableModeAtConstruction ||
        piezoEnabled != this->structuralPiezoModeAtConstruction ||
        (cableEnabled &&
         cableDynamicsEnabled != this->structuralCableDynamicsModeAtConstruction))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural cable/piezo mode cannot be enabled, disabled or switched between static and transient after constructing RSolverStress.");
    }
    if (!cableEnabled && !piezoEnabled)
    {
        return;
    }
    if (this->problemType != R_PROBLEM_STRESS)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Enabled structural cable/piezo mode requires a non-modal stress solver.");
    }
    RStructuralCableModelAdapter::validate(*this->pModel);
    RStructuralPiezoModelAdapter::validate(*this->pModel);
    const bool cableControlEnabled = cableDynamicsEnabled &&
        this->pModel->getProblemSetup().getStructuralControlSetup().enabled;
    if (this->supportResponseObserver ||
        ((this->controlCallback || this->controlStepObserver) && !cableControlEnabled))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural support response and control are supported only for transient cable NodalForce v1.");
    }
    if (cableDynamicsEnabled)
    {
        const RStructuralControlSetup &controlSetup =
            this->pModel->getProblemSetup().getStructuralControlSetup();
        if (controlSetup.enabled &&
            (!this->controlCallback || this->controllerId != controlSetup.controllerId ||
             this->controllerRevision.isEmpty()))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Enabled cable control requires a matching NodalForce callback and nonempty revision fingerprint.");
        }
        if (!controlSetup.enabled && this->controlCallback)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Registered cable control callback requires enabled NodalForce control setup.");
        const RStructuralDynamicsRestartData &checkpoint =
            this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
        const bool controlEnabled = controlSetup.enabled;
        if (checkpoint.available && controlEnabled != checkpoint.controlAvailable)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable dynamics control setup does not match the persisted accepted control ledger.");
        }
    }

    const RStructuralPiezoSetup &piezoSetup =
        this->pModel->getProblemSetup().getStructuralPiezoSetup();
    const bool thermalNodalPiezo = piezoEnabled &&
        piezoSetup.electricalMode == RStructuralPiezoElectricalMode::NodalPotential &&
        this->pModel->getProblemTaskTree().getProblemTypeMask() ==
            (R_PROBLEM_HEAT | R_PROBLEM_STRESS);

    const RRVector &sharedTemperature = this->pSharedData->findData("element-temperature");
    if (sharedTemperature.getNRows() != 0 &&
        sharedTemperature.getNRows() != this->pModel->getNElements())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Cable/piezo mode requires one element temperature per model element.");
    }
    for (uint element=0;element<sharedTemperature.getNRows();element++)
    {
        const double temperature = sharedTemperature[element];
        if (!std::isfinite(temperature) || temperature < 0.0 ||
            (!thermalNodalPiezo && std::fabs(temperature-293.15) > 1.0e-10))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         thermalNodalPiezo
                             ? "Heat->Piezo requires finite nonnegative absolute-kelvin element temperatures."
                             : "Cable/piezo mode is isothermal at 293.15 K; shared element temperature is unsupported.");
        }
    }
}

uint RSolverStress::getNComputedModes() const
{
    return uint(this->d.size());
}

void RSolverStress::initialize()
{
}

void RSolverStress::updateScales()
{
    if (this->isStructuralCableMode() || this->isStructuralPiezoMode())
    {
        this->scales = RScales();
        return;
    }
    this->scales.setMetre(this->findMeshScale());
}

bool RSolverStress::isStructuralCableMode() const
{
    return this->problemType == R_PROBLEM_STRESS &&
           this->pModel->getProblemSetup().getStructuralCableSetup().enabled;
}

bool RSolverStress::isStructuralCableDynamicsMode() const
{
    return this->isStructuralCableMode() &&
           this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
}

bool RSolverStress::isStructuralPiezoMode() const
{
    return this->problemType == R_PROBLEM_STRESS &&
           this->pModel->getProblemSetup().getStructuralPiezoSetup().enabled;
}

bool RSolverStress::isSpecializedStaticTrussMode() const
{
    return (this->isStructuralCableMode() && !this->isStructuralCableDynamicsMode()) ||
           this->isStructuralPiezoMode();
}

void RSolverStress::findComputableElements(RProblemType problemType)
{
    if (this->isSpecializedStaticTrussMode() || this->isStructuralCableDynamicsMode())
    {
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
        this->computableElements.resize(this->pModel->getNElements(),false);
        this->computableElements.fill(false);
        for (uint elementId=0;elementId<this->pModel->getNElements();elementId++)
        {
            if (this->pModel->getElement(elementId).getType() != R_ELEMENT_POINT)
            {
                this->computableElements[elementId] = true;
            }
        }
        return;
    }
    RSolverGeneric::findComputableElements(problemType);
}

void RSolverStress::recover()
{
    if (this->isStructuralCableDynamicsMode())
    {
        RStructuralCableModelAdapter::validate(*this->pModel);
        this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),0,this->nodeDisplacement.x,0.0);
        this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),1,this->nodeDisplacement.y,0.0);
        this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),2,this->nodeDisplacement.z,0.0);
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),0,this->nodeAcceleration.x,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),1,this->nodeAcceleration.y,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,
                              this->pModel->getNNodes(),2,this->nodeAcceleration.z,0.0);
        this->structuralCableDynamicsResultReady = false;
        this->structuralCableDynamicsResult = RStructuralCableDynamicsResult();
        return;
    }
    if (this->isSpecializedStaticTrussMode())
    {
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
        this->structuralCableResultReady = false;
        this->structuralCableResult = RStructuralCableEquilibriumResult();
        this->structuralPiezoResultReady = false;
        this->structuralPiezoResult = RStructuralPiezoModelResult();
        if (this->isStructuralPiezoMode())
            RStructuralPiezoModelAdapter::clearResults(*this->pModel);
        else
            RStructuralCableModelAdapter::clearResults(*this->pModel);
        return;
    }
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeDisplacement.x,0.0);
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeDisplacement.y,0.0);
    this->recoverVariable(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeDisplacement.z,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeForce.x,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeForce.y,0.0);
    this->recoverVariable(R_VARIABLE_FORCE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeForce.z,0.0);
    this->recoverVariable(R_VARIABLE_ROTATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeRotation.x,0.0);
    this->recoverVariable(R_VARIABLE_ROTATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeRotation.y,0.0);
    this->recoverVariable(R_VARIABLE_ROTATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeRotation.z,0.0);
    if (this->problemType == R_PROBLEM_STRESS && this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled())
    {
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeVelocity.x,0.0);
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeVelocity.y,0.0);
        this->recoverVariable(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeVelocity.z,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeAcceleration.x,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeAcceleration.y,0.0);
        this->recoverVariable(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeAcceleration.z,0.0);
        this->nodeAngularVelocity.x.resize(this->pModel->getNNodes(),0.0);
        this->nodeAngularVelocity.y.resize(this->pModel->getNNodes(),0.0);
        this->nodeAngularVelocity.z.resize(this->pModel->getNNodes(),0.0);
        this->nodeAngularAcceleration.x.resize(this->pModel->getNNodes(),0.0);
        this->nodeAngularAcceleration.y.resize(this->pModel->getNNodes(),0.0);
        this->nodeAngularAcceleration.z.resize(this->pModel->getNNodes(),0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeAngularVelocity.x,0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeAngularVelocity.y,0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_VELOCITY,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeAngularVelocity.z,0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeAngularAcceleration.x,0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeAngularAcceleration.y,0.0);
        this->recoverVariable(R_VARIABLE_ANGULAR_ACCELERATION,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeAngularAcceleration.z,0.0);
        if (this->firstRun)
        {
            const double initialTime = this->pModel->getTimeSolver().getPreviousTime();
            const RVariableType displacementTypes[3] = {R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
            const RVariableType velocityTypes[3] = {R_VARIABLE_VELOCITY_X,R_VARIABLE_VELOCITY_Y,R_VARIABLE_VELOCITY_Z};
            const RVariableType rotationTypes[3] = {R_VARIABLE_ROTATION_X,R_VARIABLE_ROTATION_Y,R_VARIABLE_ROTATION_Z};
            const RVariableType angularVelocityTypes[3] = {R_VARIABLE_ANGULAR_VELOCITY_X,R_VARIABLE_ANGULAR_VELOCITY_Y,R_VARIABLE_ANGULAR_VELOCITY_Z};
            RRVector *displacementComponents[3] = {&this->nodeDisplacement.x,&this->nodeDisplacement.y,&this->nodeDisplacement.z};
            RRVector *velocityComponents[3] = {&this->nodeVelocity.x,&this->nodeVelocity.y,&this->nodeVelocity.z};
            RRVector *rotationComponents[3] = {&this->nodeRotation.x,&this->nodeRotation.y,&this->nodeRotation.z};
            RRVector *angularVelocityComponents[3] = {&this->nodeAngularVelocity.x,&this->nodeAngularVelocity.y,&this->nodeAngularVelocity.z};
            for (uint groupId=0;groupId<this->pModel->getNElementGroups();groupId++)
            {
                const RElementGroup *group = this->pModel->getElementGroupPtr(groupId);
                if (!group)
                {
                    continue;
                }
                for (uint conditionId=0;conditionId<group->getNInitialConditions();conditionId++)
                {
                    const RInitialCondition &condition = group->getInitialCondition(conditionId);
                    const RVariableType *types = nullptr;
                    RRVector **components = nullptr;
                    if (condition.getType() == R_INITIAL_CONDITION_DISPLACEMENT)
                    {
                        types = displacementTypes;
                        components = displacementComponents;
                    }
                    else if (condition.getType() == R_INITIAL_CONDITION_VELOCITY)
                    {
                        types = velocityTypes;
                        components = velocityComponents;
                    }
                    else if (condition.getType() == R_INITIAL_CONDITION_ROTATION)
                    {
                        types = rotationTypes;
                        components = rotationComponents;
                    }
                    else if (condition.getType() == R_INITIAL_CONDITION_ANGULAR_VELOCITY)
                    {
                        types = angularVelocityTypes;
                        components = angularVelocityComponents;
                    }
                    if (!types)
                    {
                        continue;
                    }
                    for (uint component=0;component<3;component++)
                    {
                        const uint position = condition.findComponentPosition(types[component]);
                        if (position == RConstants::eod || !condition.getComponent(position).getEnabled())
                        {
                            continue;
                        }
                        const double value = condition.getComponent(position).get(initialTime);
                        for (uint member=0;member<group->size();member++)
                        {
                            const RElement &element = this->pModel->getElement(group->get(member));
                            for (uint localNode=0;localNode<element.size();localNode++)
                            {
                                (*components[component])[element.getNodeId(localNode)] = value;
                            }
                        }
                    }
                }
            }
        }
    }
    this->recoverVariable(R_VARIABLE_MOMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodeMoment.x,0.0);
    this->recoverVariable(R_VARIABLE_MOMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),1,this->nodeMoment.y,0.0);
    this->recoverVariable(R_VARIABLE_MOMENT,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),2,this->nodeMoment.z,0.0);
    this->recoverVariable(R_VARIABLE_PRESSURE,R_VARIABLE_APPLY_NODE,this->pModel->getNNodes(),0,this->nodePressure,0.0);

//    this->syncShared("node-displacement-x",this->nodeDisplacement.x);
//    this->syncShared("node-displacement-y",this->nodeDisplacement.y);
//    this->syncShared("node-displacement-z",this->nodeDisplacement.z);

}

void RSolverStress::prepare()
{
    if (this->isStructuralCableDynamicsMode())
    {
        RStructuralCableModelAdapter::validate(*this->pModel);
        return;
    }
    if (this->isSpecializedStaticTrussMode())
    {
        RStructuralCableModelAdapter::validate(*this->pModel);
        RStructuralPiezoModelAdapter::validate(*this->pModel);
        return;
    }

    RRVector modalAxialForce;
    bool useThermalGeometricStiffness = false;
    if (this->problemType == R_PROBLEM_STRESS_MODAL &&
        hasCurrentThermalStressBeforeModal(this->pModel->getProblemTaskTree()))
    {
        const uint elementCount = this->pModel->getNElements();
        const QString axialForceKey = QStringLiteral("structural-thermal-axial-force-v1");
        const QString temperatureKey = QStringLiteral("structural-thermal-axial-force-temperature-v1");
        if (!this->pSharedData->hasData(axialForceKey,elementCount) ||
            !this->pSharedData->hasData(temperatureKey,elementCount))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Heat->Stress->Modal requires an accepted axial-force field from the preceding stress solve.");
        }

        const RRVector &acceptedAxialForce = this->pSharedData->findData(axialForceKey);
        const RRVector &acceptedTemperature = this->pSharedData->findData(temperatureKey);
        const double forceScale = this->scales.findScaleFactor(R_VARIABLE_FORCE);
        if (!std::isfinite(forceScale) || forceScale <= 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Heat->Stress->Modal axial-force scale is invalid.");
        }

        modalAxialForce.resize(elementCount,0.0);
        for (uint elementId=0;elementId<elementCount;elementId++)
        {
            const double currentTemperature = this->elementTemperature[elementId];
            const double stressTemperature = acceptedTemperature[elementId];
            if (!std::isfinite(currentTemperature) || !std::isfinite(stressTemperature) ||
                std::fabs(currentTemperature-stressTemperature) >
                    1.0e-10*std::max({1.0,std::fabs(currentTemperature),std::fabs(stressTemperature)}))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Modal temperature differs from the temperature field accepted by the preceding stress solve.");
            }
            const RElementType type = this->pModel->getElement(elementId).getType();
            if (this->computableElements[elementId] && type != R_ELEMENT_POINT &&
                type != R_ELEMENT_TRUSS1 && type != R_ELEMENT_BEAM1)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Thermal-prestress modal geometric stiffness supports only TRUSS1 and BEAM1 structural elements; '%s' is unsupported.",
                             RElement::getName(type).toUtf8().constData());
            }
            if (!std::isfinite(acceptedAxialForce[elementId]))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Accepted thermal axial-force field is non-finite at element %u.",elementId);
            }
            if (type == R_ELEMENT_POINT && std::fabs(acceptedAxialForce[elementId]) > 0.0)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Accepted thermal axial-force field assigns nonzero force to point element %u.",elementId);
            }
            if (!this->computableElements[elementId])
            {
                continue;
            }

            const double scaledForce = acceptedAxialForce[elementId]*forceScale;
            if (!std::isfinite(scaledForce))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Accepted thermal axial force overflows the modal force scale at element %u.",elementId);
            }
            modalAxialForce[elementId] = scaledForce;
        }
        useThermalGeometricStiffness = true;
    }

    const bool structuralDynamics = this->problemType == R_PROBLEM_STRESS &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    if (structuralDynamics && !this->dynamicsIntegrator && !this->preparingInitialForce)
    {
        // Assemble f(t0) separately. Reusing f(t1) changes the initial
        // acceleration for a ramp, impulse or user-supplied time history.
        RTimeSolver &time = this->pModel->getTimeSolver();
        const uint step = time.getCurrentTimeStep();
        const double endTime = time.getCurrentTime();
        const double startTime = time.getPreviousTime();
        this->preparingInitialForce = true;
        time.getTimes()[step] = startTime;
        try
        {
            this->prepare();
            this->initialDynamicForce = this->b;
            this->dynamicMassReference = this->M;
            this->dynamicStiffnessReference = this->A;
            this->dynamicMotionMassReference = this->motionMass;
            this->dynamicMotionStiffnessReference = this->motionStiffness;
            this->dynamicConstraintNullspaceReference = this->constraintNullspace;
            this->dynamicConstraintParticularReference = this->constraintParticular;
            this->initialSupportForce = this->supportForce;
            this->supportStiffnessReference = this->supportStiffness;
            this->supportMassReference = this->supportMass;
        }
        catch (...)
        {
            time.getTimes()[step] = endTime;
            this->preparingInitialForce = false;
            throw;
        }
        time.getTimes()[step] = endTime;
        this->preparingInitialForce = false;
    }
    const bool needsMass = structuralDynamics || this->pModel->getTimeSolver().getEnabled() || this->problemType == R_PROBLEM_STRESS_MODAL;
    const bool lumpedMass = this->pModel->getProblemSetup().getMassMatrixType() == RMassMatrixType::Lumped;

    //! Element force vector.
    RSolverCartesianVector<RRVector> elementForce;
    RSolverCartesianVector<RBVector> forceSetValues;

    //! Element moment vector. Moments are currently accepted on point
    //! entities and assembled directly into frame rotational equations.
    RSolverCartesianVector<RRVector> elementMoment;
    RSolverCartesianVector<RBVector> momentSetValues;


    //! Element gravity vector.
    RSolverCartesianVector<RRVector> elementGravity;
    RSolverCartesianVector<RBVector> gGravitySetValues;

    //! Element pressure.
    RRVector elementPressure;
    RBVector pressureSetValues;

    //! Element traction per unit area.
    RSolverCartesianVector<RRVector> elementForceUnitArea;
    RSolverCartesianVector<RBVector> forceUnitAreaSetValues;

    //! Element weight.
    RRVector elementWeight;
    RBVector weightSetValues;

    RBVector temperatureSetValues;

    // The node frames and the directions they hold come from the constraints
    // themselves, so they have to be known before the equation book is built.
    this->generateLocalConstraints();
    if (structuralDynamics)
    {
        this->validateDynamicBoundaryState();
    }
    this->generateEquationBook();
    this->nLinearConstraintEquations = 0;

    const RStructuralInteractionSetup &interactionSetup = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    const RStructuralValidationReport validation = interactionSetup.validate(*this->pModel);
    for (const RStructuralDiagnostic &diagnostic : validation.getDiagnostics())
    {
        if (diagnostic.severity == RStructuralDiagnosticSeverity::Warning)
        {
            RLogger::warning("Structural setup [%s] %s\n",
                             diagnostic.code.toUtf8().constData(),
                             diagnostic.message.toUtf8().constData());
        }
        else if (diagnostic.severity == RStructuralDiagnosticSeverity::Error)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural setup '%s' is invalid [%s]: %s",
                         diagnostic.objectId.toUtf8().constData(),
                         diagnostic.code.toUtf8().constData(),
                         diagnostic.message.toUtf8().constData());
        }
    }
    if (structuralDynamics)
    {
        this->preparePrescribedMotions();
        this->prepareSupportResponse();
    }
    this->generateVariableVector(R_VARIABLE_FORCE_X,elementForce.x,forceSetValues.x,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_FORCE_Y,elementForce.y,forceSetValues.y,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_FORCE_Z,elementForce.z,forceSetValues.z,true,this->firstRun,this->firstRun);
    this->generateVariableVector(R_VARIABLE_MOMENT_X,elementMoment.x,momentSetValues.x,true,false,false);
    this->generateVariableVector(R_VARIABLE_MOMENT_Y,elementMoment.y,momentSetValues.y,true,false,false);
    this->generateVariableVector(R_VARIABLE_MOMENT_Z,elementMoment.z,momentSetValues.z,true,false,false);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_X,elementGravity.x,gGravitySetValues.x,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Y,elementGravity.y,gGravitySetValues.y,true,true,true);
    this->generateVariableVector(R_VARIABLE_G_ACCELERATION_Z,elementGravity.z,gGravitySetValues.z,true,true,true);
    this->generateVariableVector(R_VARIABLE_PRESSURE,elementPressure,pressureSetValues,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_X,elementForceUnitArea.x,forceUnitAreaSetValues.x,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_Y,elementForceUnitArea.y,forceUnitAreaSetValues.y,true,true,true);
    this->generateVariableVector(R_VARIABLE_FORCE_UNIT_AREA_Z,elementForceUnitArea.z,forceUnitAreaSetValues.z,true,true,true);
    this->generateVariableVector(R_VARIABLE_WEIGHT,elementWeight,weightSetValues,true,true,true);
    this->generateMaterialVecor(RMaterialProperty::ModulusOfElasiticity,this->elementElasticityModulus);
    this->generateMaterialVecor(RMaterialProperty::PoissonRatio,this->elementPoissonRatio);
    this->generateMaterialVecor(RMaterialProperty::Density,this->elementDensity);
    this->generateMaterialVecor(RMaterialProperty::ThermalExpansionCoefficient,this->elementThermalExpansion);
    this->generateVariableVector(R_VARIABLE_TEMPERATURE,this->elementEnvironmentTemperature,temperatureSetValues,false,false,true);
    if (this->elementTemperature.getNRows() == this->pModel->getNElements() &&
        this->elementEnvironmentTemperature.getNRows() == this->pModel->getNElements())
    {
        for (uint groupId=0;groupId<this->pModel->getNElementGroups();groupId++)
        {
            const RElementGroup *group = this->pModel->getElementGroupPtr(groupId);
            if (!group) continue;
            const RMaterial &material = group->getMaterial();
            const uint propertyPosition = material.findPosition(
                RMaterialProperty::ThermalExpansionCoefficient);
            if (propertyPosition >= material.size()) continue;
            const RMaterialProperty &property = material.get(propertyPosition);
            for (uint position=0;position<group->size();position++)
            {
                const uint elementId = group->get(position);
                if (elementId >= this->pModel->getNElements())
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Thermal expansion group references invalid element %u.",elementId);
                this->elementThermalExpansion[elementId] = meanThermalExpansion(
                    property,this->elementEnvironmentTemperature[elementId],
                    this->elementTemperature[elementId]);
            }
        }
    }

    this->b.resize(this->equationBook.getNEnabled());
    this->x.resize(this->equationBook.getNEnabled());

    this->M.clear();
    this->M.setNRows(this->b.size());
    this->A.clear();
    this->A.setNRows(this->b.size());
    this->motionStiffness.clear();
    this->motionMass.clear();
    const bool hasPrescribedMotion = structuralDynamics && !this->prescribedMotions.empty();
    if (hasPrescribedMotion)
    {
        this->motionStiffness.setNRows(this->b.size());
        this->motionMass.setNRows(this->b.size());
    }
    this->b.fill(0.0);
    this->x.fill(0.0);

    // Per-thread assembly buffers - elements are assembled without
    // synchronization and merged into A/b (and M for modal) at the end.
    int np = omp_get_max_threads();
    std::vector<RSparseMatrix> Ap(np);
    std::vector<RSparseMatrix> Mp(np);
    std::vector<RSparseMatrix> Kfcp(np);
    std::vector<RSparseMatrix> Mfcp(np);
    std::vector<RRVector> bp(np);
    this->supportStiffnessBuffers.resize(np);
    this->supportMassBuffers.resize(np);
    this->supportForceBuffers.resize(np);
    for (int t=0;t<np;t++)
    {
        Ap[t].setNRows(this->b.size());
        Mp[t].setNRows(this->b.size());
        if (hasPrescribedMotion)
        {
            Kfcp[t].setNRows(this->b.size());
            Mfcp[t].setNRows(this->b.size());
        }
        bp[t].resize(this->b.size());
        bp[t].fill(0.0);
        this->supportStiffnessBuffers[t].clear();
        this->supportMassBuffers[t].clear();
        this->supportStiffnessBuffers[t].setNRows(uint(this->heldDofs.size()));
        this->supportMassBuffers[t].setNRows(uint(this->heldDofs.size()));
        this->supportForceBuffers[t].resize(uint(this->heldDofs.size()),0.0);
        this->supportForceBuffers[t].fill(0.0);
    }

    this->pModel->convertElementToNodeVector(elementForce.x,forceSetValues.x,this->nodeForce.x,true);
    this->pModel->convertElementToNodeVector(elementForce.y,forceSetValues.y,this->nodeForce.y,true);
    this->pModel->convertElementToNodeVector(elementForce.z,forceSetValues.z,this->nodeForce.z,true);
    this->pModel->convertElementToNodeVector(elementPressure,pressureSetValues,this->nodePressure,true);

    // Convert node pressure to element pressure
    for (uint i=0;i<elementPressure.size();i++)
    {
        if (!pressureSetValues[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            uint nne = rElement.size();
            if (nne > 0)
            {
                elementPressure[i] = 0.0;
                for (uint j=0;j<nne;j++)
                {
                    elementPressure[i] += this->nodePressure[rElement.getNodeId(j)];
                }
                elementPressure[i] /= double(nne);
            }
        }
    }

    std::vector<RNode> nodesBkp;
    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        this->nodeInitialDisplacement.x = this->nodeDisplacement.x;
        this->nodeInitialDisplacement.y = this->nodeDisplacement.y;
        this->nodeInitialDisplacement.z = this->nodeDisplacement.z;
        RLogger::info("Moving prestressed nodes\n");
        nodesBkp = this->pModel->getNodes();
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            this->pModel->getNode(i).move(RR3Vector(this->nodeDisplacement.x[i],this->nodeDisplacement.y[i],this->nodeDisplacement.z[i]));
        }
    }

    // Prepare point elements.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        RPoint &point = this->pModel->getPoint(i);
        double pointVolume = point.getVolume();
        // Force and Weight are totals over the entity, exactly as they are for
        // a line or a surface, so they are spread over its point elements.
        double pointCount = double(std::max(point.size(),uint(1)));

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(point.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = point.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_POINT(this->pModel->getElement(elementID).getType()));
                RRMatrix Me(3,3);
                RRMatrix Ke(3,3);
                RRVector fe(3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                // Force
                fe[0] += elementForce.x[elementID] / pointCount;
                fe[1] += elementForce.y[elementID] / pointCount;
                fe[2] += elementForce.z[elementID] / pointCount;
                // Weight
                fe[0] += elementWeight[elementID] * elementGravity.x[elementID] / pointCount;
                fe[1] += elementWeight[elementID] * elementGravity.y[elementID] / pointCount;
                fe[2] += elementWeight[elementID] * elementGravity.z[elementID] / pointCount;
                // Own weight
                if (pointVolume > 0.0)
                {
                    fe[0] += elementGravity.x[elementID] * this->elementDensity[elementID] * pointVolume;
                    fe[1] += elementGravity.y[elementID] * this->elementDensity[elementID] * pointVolume;
                    fe[2] += elementGravity.z[elementID] * this->elementDensity[elementID] * pointVolume;
                }

                // Mass
                if (needsMass)
                {
                    Me.setIdentity(3);
                    Me *= this->elementDensity[elementID] * pointVolume;
                }

                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())],Mp[uint(omp_get_thread_num())],Kfcp[uint(omp_get_thread_num())],Mfcp[uint(omp_get_thread_num())]);

                if (this->problemType != R_PROBLEM_STRESS_MODAL)
                {
                    const uint nodeId = this->pModel->getElement(elementID).getNodeId(0);
                    RR3Vector moment(elementMoment.x[elementID]/pointCount,
                                     elementMoment.y[elementID]/pointCount,
                                     elementMoment.z[elementID]/pointCount);
                    if (this->rotationLocalRotations[nodeId].isActive())
                    {
                        const RRMatrix &rotation = this->rotationLocalRotations[nodeId].getR();
                        const RR3Vector globalMoment(moment);
                        for (uint c=0;c<3;c++)
                        {
                            moment[c] = rotation[0][c]*globalMoment[0]
                                      + rotation[1][c]*globalMoment[1]
                                      + rotation[2][c]*globalMoment[2];
                        }
                    }
                    for (uint c=0;c<3;c++)
                    {
                        const RDofType dof = c == 0 ? RDofType::RX
                                            : c == 1 ? RDofType::RY
                                                     : RDofType::RZ;
                        uint position = 0;
                        if (this->equationBook.getValue(nodeId,dof,position))
                        {
                            bp[uint(omp_get_thread_num())][position] += moment[c];
                        }
                    }
                }
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
        double lineCrossArea = line.getCrossArea();

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = line.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));

                double lineLength = 0.0;
                element.findLength(this->pModel->getNodes(),lineLength);
                if (element.getType() == R_ELEMENT_BEAM1 || element.getType() == R_ELEMENT_BEAM2)
                {
                    if (lineLength <= RConstants::eps)
                    {
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Beam element has zero end-to-end length.");
                    }
                    RRMatrix Me,Ke;
                    RRVector fe;
                    RR3Vector distributedForce(
                        elementForce.x[elementID]/lineLength
                            + elementWeight[elementID]*elementGravity.x[elementID]/lineLength
                            + elementGravity.x[elementID]*this->elementDensity[elementID]*lineCrossArea,
                        elementForce.y[elementID]/lineLength
                            + elementWeight[elementID]*elementGravity.y[elementID]/lineLength
                            + elementGravity.y[elementID]*this->elementDensity[elementID]*lineCrossArea,
                        elementForce.z[elementID]/lineLength
                            + elementWeight[elementID]*elementGravity.z[elementID]/lineLength
                            + elementGravity.z[elementID]*this->elementDensity[elementID]*lineCrossArea);
                    buildBeamElement(line,element,this->pModel->getNodes(),
                                     this->elementElasticityModulus[elementID],
                                     this->elementPoissonRatio[elementID],
                                     this->elementDensity[elementID],
                                     this->elementThermalExpansion[elementID],
                                     this->elementTemperature[elementID]-this->elementEnvironmentTemperature[elementID],
                                     distributedForce,needsMass,lumpedMass,Me,Ke,fe);
                    if (useThermalGeometricStiffness)
                    {
                        addBeam1GeometricStiffness(line,element,this->pModel->getNodes(),
                                                   modalAxialForce[elementID],Ke);
                    }
                    this->assemblyMatrix(elementID,Me,Ke,fe,
                                         Ap[uint(omp_get_thread_num())],
                                         bp[uint(omp_get_thread_num())],
                                         Mp[uint(omp_get_thread_num())],
                                         Kfcp[uint(omp_get_thread_num())],
                                         Mfcp[uint(omp_get_thread_num())]);
                    continue;
                }
                if (element.getType() != R_ELEMENT_TRUSS1 && element.getType() != R_ELEMENT_TRUSS2)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Line element type '%s' is not implemented by the structural solver.",
                                 RElement::getName(element.getType()).toUtf8().constData());
                }
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3,0.0);
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);

                RRMatrix Be(3*element.size(),1);
                RRMatrix BeT(1,3*element.size());
                RRMatrix Met(element.size(),element.size(),0.0);
                RRMatrix MeRt(element.size()*3,element.size(),0.0);

                double E = this->elementElasticityModulus[elementID];
                double De = E * lineCrossArea;

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    if (needsMass && lineCrossArea > 0.0)
                    {
                        Met.fill(0.0);
                        for (uint m=0;m<element.size();m++)
                        {
                            for (uint n=0;n<element.size();n++)
                            {
                                Met[m][n] = N[m] * N[n]
                                          * this->elementDensity[elementID]
                                          * detJ
                                          * shapeFunc.getW()
                                          * lineCrossArea;
                            }
                        }
                        if (lumpedMass)
                        {
                            lumpMassMatrix(Met);
                        }
                        RRMatrix RtT;
                        RtT.transpose(Rt);
                        RRMatrix::mlt(Rt,Met,MeRt);
                        RRMatrix::mlt(MeRt,RtT,Me,true);
                    }
                    if (lineCrossArea > 0.0)
                    {
                        // Strain-displacement vector of the truss - it maps the
                        // global nodal displacements onto the axial strain.
                        Be.fill(0.0);
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            Be[3*m+0][0] = Rt[3*m+0][m]*dN[m][0]*J[0][0];
                            Be[3*m+1][0] = Rt[3*m+1][m]*dN[m][0]*J[0][0];
                            Be[3*m+2][0] = Rt[3*m+2][m]*dN[m][0]*J[0][0];
                        }
                        BeT.transpose(Be);

                        // Ke += E*A * B * B^T * detJ * W, accumulated over the
                        // integration points.
                        RRMatrix BeScaled(Be);
                        BeScaled *= De * detJ * shapeFunc.getW();
                        RRMatrix::mlt(BeScaled,BeT,Ke,true);

                        // Thermal expansion force: f += E*A * alpha * dT * B * detJ * W
                        double thermalFactor = De
                                             * this->elementThermalExpansion[elementID]
                                             * dT
                                             * detJ
                                             * shapeFunc.getW();
                        for (uint m=0;m<3*element.size();m++)
                        {
                            fe[m] += thermalFactor * Be[m][0];
                        }
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        double integValue = N[m] * detJ * shapeFunc.getW();

                        // Force
                        fe[3*m+0] += (elementForce.x[elementID] / lineLength) * integValue;
                        fe[3*m+1] += (elementForce.y[elementID] / lineLength) * integValue;
                        fe[3*m+2] += (elementForce.z[elementID] / lineLength) * integValue;
                        // Weight
                        fe[3*m+0] += (elementWeight[elementID] * elementGravity.x[elementID] / lineLength) * integValue;
                        fe[3*m+1] += (elementWeight[elementID] * elementGravity.y[elementID] / lineLength) * integValue;
                        fe[3*m+2] += (elementWeight[elementID] * elementGravity.z[elementID] / lineLength) * integValue;
                        // Own weight
                        if (lineCrossArea > 0.0)
                        {
                            fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                            fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                            fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * lineCrossArea * integValue;
                        }

                    }
                }
                if (useThermalGeometricStiffness)
                {
                    addTrussGeometricStiffness(element,this->pModel->getNodes(),
                                               modalAxialForce[elementID],Ke);
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())],Mp[uint(omp_get_thread_num())],Kfcp[uint(omp_get_thread_num())],Mfcp[uint(omp_get_thread_num())]);
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
        double surfaceArea = surface.findArea(this->pModel->getNodes(),this->pModel->getElements());
        double surfaceThickness = surface.getThickness();
        validateSurfaceFormulation(surface);
        const bool axisymmetric = surface.getStructuralFormulation() == RElementFormulation::Axisymmetric;
        const bool structuralSurface = axisymmetric || surfaceThickness > 0.0;
        const double surfaceMeasure = axisymmetric
                                    ? findAxisymmetricMeasure(surface,this->pModel->getNodes(),this->pModel->getElements())
                                    : surfaceArea;

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = surface.get(uint(j));

                if (!this->computableElements[elementID] && !this->includableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                if (axisymmetric &&
                    (std::fabs(elementPressure[elementID]) > 0.0 ||
                     std::fabs(elementForceUnitArea.x[elementID]) > 0.0 ||
                     std::fabs(elementForceUnitArea.y[elementID]) > 0.0 ||
                     std::fabs(elementForceUnitArea.z[elementID]) > 0.0))
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Axisymmetric pressure and surface-traction loads must be applied "
                                 "on meridian boundary lines; area loads on surface '%s' are ambiguous.",
                                 surface.getName().toUtf8().constData());
                }
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3);
                RRMatrix Ke(element.size()*3,element.size()*3);
                RRVector fe(element.size()*3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                const uint nStrainComponents = axisymmetric ? 4u : 3u;
                RRMatrix B(element.size(),2,0.0);
                RRMatrix Be(element.size()*2,nStrainComponents,0.0);
                RRMatrix BeT(nStrainComponents,element.size()*2,0.0);
                RRMatrix BeD(element.size()*2,nStrainComponents,0.0);
                RRMatrix Met(element.size()*2,element.size()*2);
                RRMatrix MeRt(element.size()*3,element.size()*2);
                RRMatrix Ket(element.size()*2,element.size()*2);
                RRMatrix KeRt(element.size()*3,element.size()*2);
                RRVector fet(element.size()*2);

                RRMatrix De;

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                buildSurfaceConstitutiveMatrix(surface.getStructuralFormulation(),E,v,De);

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                RR3Vector normal(0.0,0.0,0.0);
                if (!axisymmetric)
                {
                    element.findNormal(this->pModel->getNodes(),normal[0],normal[1],normal[2]);
                }

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double radius = 0.0;
                    double detJ = 0.0;
                    if (axisymmetric)
                    {
                        detJ = findAxisymmetricKinematics(element,this->pModel->getNodes(),shapeFunc,
                                                          B,Rt,radius);
                    }
                    else
                    {
                        detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                        B.fill(0.0);
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m][0] = dN[m][0]*J[0][0] + dN[m][1]*J[0][1];
                            B[m][1] = dN[m][0]*J[1][0] + dN[m][1]*J[1][1];
                        }
                    }
                    RtT.transpose(Rt);
                    const double integrationScale = (axisymmetric
                                                    ? axisymmetricFullAngle*radius
                                                    : surfaceThickness)
                                                  * detJ*shapeFunc.getW();

                    fet.fill(0.0);
                    if (needsMass)
                    {
                        Met.fill(0.0);
                    }

                    if (structuralSurface)
                    {
                        Be.fill(0.0);
                        for (uint m=0;m<element.size();m++)
                        {
                            if (axisymmetric)
                            {
                                Be[2*m+0][0] = B[m][0];
                                Be[2*m+0][2] = N[m]/radius;
                                Be[2*m+0][3] = B[m][1];
                                Be[2*m+1][1] = B[m][1];
                                Be[2*m+1][3] = B[m][0];
                            }
                            else
                            {
                                Be[2*m+0][0] = B[m][0];
                                Be[2*m+1][1] = B[m][1];
                                Be[2*m+0][2] = B[m][1];
                                Be[2*m+1][2] = B[m][0];
                            }
                        }
                        BeT.transpose(Be);

                        RRMatrix::mlt(Be,De,BeD);
                        RRMatrix::mlt(BeD,BeT,Ket);
                        Ket *= integrationScale;
                        RRMatrix::mlt(Rt,Ket,KeRt);
                        RRMatrix::mlt(KeRt,RtT,Ke,true);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        if (structuralSurface)
                        {
                            // Mass
                            if (needsMass)
                            {
                                for (uint n=0;n<element.size();n++)
                                {
                                    double value = N[m] * N[n]
                                                 * this->elementDensity[elementID]
                                                 * integrationScale;
                                    Met[2*m+0][2*n+0] += value;
                                    Met[2*m+1][2*n+1] += value;
                                }
                            }
                        }

                        double integValue = N[m] * (axisymmetric
                                                  ? axisymmetricFullAngle*radius*detJ*shapeFunc.getW()
                                                  : detJ*shapeFunc.getW());

                        // Pressure vector
                        if (!axisymmetric)
                        {
                            fe[3*m+0] += elementPressure[elementID] * normal[0] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                            fe[3*m+1] += elementPressure[elementID] * normal[1] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                            fe[3*m+2] += elementPressure[elementID] * normal[2] * integValue * (this->inwardElements[elementID] ? 1.0 : -1.0);
                        }
                        // Force per unit area
                        fe[3*m+0] += elementForceUnitArea.x[elementID] * integValue;
                        fe[3*m+1] += elementForceUnitArea.y[elementID] * integValue;
                        fe[3*m+2] += elementForceUnitArea.z[elementID] * integValue;
                        // Force
                        fe[3*m+0] += (elementForce.x[elementID] / surfaceMeasure) * integValue;
                        fe[3*m+1] += (elementForce.y[elementID] / surfaceMeasure) * integValue;
                        fe[3*m+2] += (elementForce.z[elementID] / surfaceMeasure) * integValue;
                        // Weight
                        fe[3*m+0] += (elementWeight[elementID] * elementGravity.x[elementID] / surfaceMeasure) * integValue;
                        fe[3*m+1] += (elementWeight[elementID] * elementGravity.y[elementID] / surfaceMeasure) * integValue;
                        fe[3*m+2] += (elementWeight[elementID] * elementGravity.z[elementID] / surfaceMeasure) * integValue;
                        // Own weight
                        if (structuralSurface)
                        {
                            const double thicknessFactor = axisymmetric ? 1.0 : surfaceThickness;
                            fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * thicknessFactor * integValue;
                            fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * thicknessFactor * integValue;
                            fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * thicknessFactor * integValue;
                        }

                        // Thermal expansion
                        if (structuralSurface)
                        {
                            const double thermalScale = this->elementThermalExpansion[elementID]
                                                      * dT * integrationScale;
                            const uint nThermalNormals = axisymmetric ? 3u : 2u;
                            for (uint n=0;n<nThermalNormals;n++)
                            {
                                fet[2*m+0] += thermalScale * BeD[2*m+0][n];
                                fet[2*m+1] += thermalScale * BeD[2*m+1][n];
                            }
                        }
                    }

                    if (structuralSurface)
                    {
                        RRMatrix::mlt(Rt,fet,fe,true);
                    }

                    // Mass
                    if (structuralSurface && needsMass)
                    {
                        if (lumpedMass)
                        {
                            lumpMassMatrix(Met);
                        }
                        RRMatrix::mlt(Rt,Met,MeRt);
                        RRMatrix::mlt(MeRt,RtT,Me,true);
                    }
                    if (!this->computableElements[elementID])
                    {
                        Me.fill(0.0);
                        Ke.fill(0.0);
                    }
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())],Mp[uint(omp_get_thread_num())],Kfcp[uint(omp_get_thread_num())],Mfcp[uint(omp_get_thread_num())]);
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
                uint elementID = volume.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Me(element.size()*3,element.size()*3);
                RRMatrix Ke(element.size()*3,element.size()*3);
                RRVector fe(element.size()*3);

                Me.fill(0.0);
                Ke.fill(0.0);
                fe.fill(0.0);

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*3,6);
                RRMatrix BeT(6,element.size()*3);
                RRMatrix BeD(element.size()*3,6);
                RRMatrix Ket(element.size()*3,element.size()*3);

                RRMatrix De(6,6);
                De.fill(0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1.0-v; De[0][1] = v;     De[0][2] = v;
                De[1][0] = v;     De[1][1] = 1.0-v; De[1][2] = v;
                De[2][0] = v;     De[2][1] = v;     De[2][2] = 1.0-v;
                De[3][3] = De[4][4] = De[5][5] = (1.0-2.0*v)/2.0;
                De *= E/((1.0+v)*(1.0-2.0*v));

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        Be[3*m+0][0] = B[m][0];   Be[3*m+1][0] = 0.0;       Be[3*m+2][0] = 0.0;
                        Be[3*m+0][1] = 0.0;       Be[3*m+1][1] = B[m][1];   Be[3*m+2][1] = 0.0;
                        Be[3*m+0][2] = 0.0;       Be[3*m+1][2] = 0.0;       Be[3*m+2][2] = B[m][2];
                        Be[3*m+0][3] = 0.0;       Be[3*m+1][3] = B[m][2];   Be[3*m+2][3] = B[m][1];
                        Be[3*m+0][4] = B[m][2];   Be[3*m+1][4] = 0.0;       Be[3*m+2][4] = B[m][0];
                        Be[3*m+0][5] = B[m][1];   Be[3*m+1][5] = B[m][0];   Be[3*m+2][5] = 0.0;
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    for (uint m=0;m<3*element.size();m++)
                    {
                        for (uint n=0;n<3*element.size();n++)
                        {
                            // Stiffness matrix
                            Ke[m][n] += Ket[m][n] * detJ * shapeFunc.getW();
                        }
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<element.size();n++)
                        {
                            // Mass
                            if (needsMass)
                            {
                                double value = N[m] * N[n]
                                             * this->elementDensity[elementID]
                                             * detJ
                                             * shapeFunc.getW();
                                Me[3*m+0][3*n+0] += value;
                                Me[3*m+1][3*n+1] += value;
                                Me[3*m+2][3*n+2] += value;
                            }
                        }

                        // Own weight
                        fe[3*m+0] += elementGravity.x[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();
                        fe[3*m+1] += elementGravity.y[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();
                        fe[3*m+2] += elementGravity.z[elementID] * this->elementDensity[elementID] * N[m] * detJ * shapeFunc.getW();

                        // Thermal expansion
                        for (uint n=0;n<3;n++)
                        {
                            fe[3*m+0] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+0][n] * detJ * shapeFunc.getW();
                            fe[3*m+1] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+1][n] * detJ * shapeFunc.getW();
                            fe[3*m+2] += this->elementThermalExpansion[elementID] * dT * BeD[3*m+2][n] * detJ * shapeFunc.getW();
                        }
                    }
                }
                if (needsMass && lumpedMass)
                {
                    lumpMassMatrix(Me);
                }
                this->assemblyMatrix(elementID,Me,Ke,fe,Ap[uint(omp_get_thread_num())],bp[uint(omp_get_thread_num())],Mp[uint(omp_get_thread_num())],Kfcp[uint(omp_get_thread_num())],Mfcp[uint(omp_get_thread_num())]);
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
            this->M.getVector(uint(i)).addVector(Mp[t].getVector(uint(i)));
            if (hasPrescribedMotion)
            {
                this->motionStiffness.getVector(uint(i)).addVector(Kfcp[t].getVector(uint(i)));
                this->motionMass.getVector(uint(i)).addVector(Mfcp[t].getVector(uint(i)));
            }
            this->b[uint(i)] += bp[t][uint(i)];
        }
    }

    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        RLogger::info("Restoring prestressed nodes\n");
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            this->pModel->setNode(i,nodesBkp[i]);
        }
    }
    if (this->supportResponseObserver)
    {
        this->supportStiffness.clear();
        this->supportMass.clear();
        this->supportStiffness.setNRows(uint(this->heldDofs.size()));
        this->supportMass.setNRows(uint(this->heldDofs.size()));
        this->supportForce.resize(uint(this->heldDofs.size()),0.0);
        this->supportForce.fill(0.0);
        for (uint row=0;row<this->heldDofs.size();row++)
        {
            for (int t=0;t<np;t++)
            {
                this->supportStiffness.getVector(row).addVector(this->supportStiffnessBuffers[t].getVector(row));
                this->supportMass.getVector(row).addVector(this->supportMassBuffers[t].getVector(row));
                this->supportForce[row] += this->supportForceBuffers[t][row];
            }
        }
    }
    this->applyLinearConstraints();
}

void RSolverStress::solveStructuralCableDynamics()
{
    if (!this->isStructuralCableDynamicsMode())
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Cable Newmark step requested outside transient TRUSS1 cable mode.");
    RStructuralCableModelAdapter::validate(*this->pModel);

    const std::size_t nDofs = std::size_t(3)*this->pModel->getNNodes();
    std::vector<double> previousDisplacements(nDofs,0.0);
    std::vector<double> previousVelocities(nDofs,0.0);
    std::vector<double> previousAccelerations(nDofs,0.0);
    for (uint node=0;node<this->pModel->getNNodes();node++)
    {
        previousDisplacements[3*node] = this->nodeDisplacement.x[node];
        previousDisplacements[3*node+1] = this->nodeDisplacement.y[node];
        previousDisplacements[3*node+2] = this->nodeDisplacement.z[node];
        previousVelocities[3*node] = this->nodeVelocity.x[node];
        previousVelocities[3*node+1] = this->nodeVelocity.y[node];
        previousVelocities[3*node+2] = this->nodeVelocity.z[node];
        previousAccelerations[3*node] = this->nodeAcceleration.x[node];
        previousAccelerations[3*node+1] = this->nodeAcceleration.y[node];
        previousAccelerations[3*node+2] = this->nodeAcceleration.z[node];
    }

    const RStructuralControlSetup &controlSetup =
        this->pModel->getProblemSetup().getStructuralControlSetup();
    if (controlSetup.enabled)
    {
        const double controlScales[] = {
            this->scales.getMetre(),
            this->scales.getKilogram(),
            this->scales.getSecond(),
            this->scales.findScaleFactor(R_VARIABLE_DISPLACEMENT_X),
            this->scales.findScaleFactor(R_VARIABLE_VELOCITY_X),
            this->scales.findScaleFactor(R_VARIABLE_ACCELERATION_X),
            this->scales.findScaleFactor(R_VARIABLE_FORCE_X)
        };
        if (!std::all_of(controlScales,controlScales+sizeof(controlScales)/sizeof(controlScales[0]),
                         [](double value){return std::isfinite(value) && value == 1.0;}))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Transient cable NodalForce control currently requires unity metre/kilogram/second FE scales.");
        }
        if (!this->controlCallback)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Enabled cable structural control has no registered NodalForce provider.");
        if (!this->controlForceMap)
            this->prepareControlMapping();
        if (!this->cachedControlCatalogFingerprint.isEmpty() &&
            this->controlCatalogFingerprint() != this->cachedControlCatalogFingerprint)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control callback changed the frozen sensor/actuator catalog.");

        const double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
        const double timeScale = this->scales.getSecond();
        if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(timeScale) || timeScale <= 0.0)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable NodalForce control requires a finite positive physical time step.");
        const double acceptedModelTime = this->pModel->getTimeSolver().getPreviousTime();
        const double nextModelTime = this->pModel->getTimeSolver().getCurrentTime();
        const double physicalDt = dt/timeScale;
        if (!std::isfinite(acceptedModelTime) || !std::isfinite(nextModelTime) ||
            std::fabs((nextModelTime-acceptedModelTime)-dt) >
                1.0e-10*std::max({1.0,std::fabs(nextModelTime),
                                   std::fabs(acceptedModelTime),dt}))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control observation time does not match the accepted Newmark interval.");
        }

        RStructuralDynamicsState accepted;
        accepted.time = acceptedModelTime;
        accepted.q = RRVector(previousDisplacements);
        accepted.v = RRVector(previousVelocities);
        accepted.a = RRVector(previousAccelerations);
        if (!std::all_of(previousDisplacements.begin(),previousDisplacements.end(),
                         [](double value){return std::isfinite(value);}) ||
            !std::all_of(previousVelocities.begin(),previousVelocities.end(),
                         [](double value){return std::isfinite(value);}) ||
            !std::all_of(previousAccelerations.begin(),previousAccelerations.end(),
                         [](double value){return std::isfinite(value);}))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control cannot observe non-finite accepted nodal kinematics.");
        }
        RStructuralControlObservation observation;
        observation.acceptedTime = accepted.time/timeScale;
        observation.nextTime = nextModelTime/timeScale;
        observation.deterministicSeed = controlSetup.deterministicSeed;
        observation.sensors = this->sampleControlSensors(accepted);
        for (const RStructuralControlSensorSample &sensor : observation.sensors)
        {
            if (!std::isfinite(sensor.value))
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable control sensor '%s' observed a non-finite accepted value.",
                             sensor.sensorId.toUtf8().constData());
        }
        observation.priorControllerState = this->acceptedControllerState;

        const RTimeSolver frozenTimeSolver = this->pModel->getTimeSolver();
        const RStructuralControlSetup frozenControlSetup = controlSetup;
        const RStructuralControlDecision decision = this->controlCallback(observation);
        if (!sameControlTimeSolver(frozenTimeSolver,this->pModel->getTimeSolver()) ||
            !sameControlSetup(frozenControlSetup,
                this->pModel->getProblemSetup().getStructuralControlSetup()))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control callback changed the frozen time solver or control setup.");
        if (this->controlCatalogFingerprint() != this->cachedControlCatalogFingerprint)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural control callback changed the frozen sensor/actuator catalog.");
        if (this->controlMappingFingerprint() != this->cachedControlFingerprint)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control callback changed the frozen physical system or NodalForce mapping.");
        if (decision.requestedCommands.size() != this->activeControlCommands.size() ||
            decision.nextControllerState.size() > 1024*1024)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control provider returned the wrong NodalForce channel count or oversized controller state.");

        const RStructuralControlConditioner::Result conditioned =
            this->controlConditioner->advance(this->acceptedConditionerState,
                                              decision.requestedCommands);
        this->candidateConditionerState = conditioned.candidate;
        std::vector<RStructuralControlCommand> applied;
        applied.reserve(conditioned.appliedCommands.size());
        const auto &actuators = this->pModel->getProblemSetup()
            .getStructuralInteractionSetup().getActuators();
        bool changed = false;
        for (size_t i=0;i<conditioned.appliedCommands.size();i++)
        {
            const RStructuralControlCommand &request = conditioned.appliedCommands[i];
            const RStructuralActuator &actuator = actuators[i];
            if (request.actuatorId != this->activeControlCommands[i].actuatorId ||
                !std::isfinite(request.value) || std::isnan(actuator.minimum) ||
                std::isnan(actuator.maximum) || std::isnan(actuator.rateLimit) ||
                actuator.minimum > actuator.maximum || actuator.rateLimit < 0.0)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable NodalForce command or actuator bounds are invalid for '%s'.",
                             actuator.id.toUtf8().constData());
            }
            const double previous = this->controlLedger->accepted()[i].command;
            const double maxDelta = actuator.rateLimit*physicalDt;
            const double lower = std::max(actuator.minimum,previous-maxDelta);
            const double upper = std::min(actuator.maximum,previous+maxDelta);
            if (lower > upper)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable NodalForce actuator '%s' has incompatible magnitude and rate limits.",
                             actuator.id.toUtf8().constData());
            const double value = std::clamp(request.value,lower,upper);
            changed = changed || value != this->activeControlCommands[i].value;
            applied.push_back({request.actuatorId,value});
        }
        if (changed)
        {
            RRVector force(uint(nDofs),0.0);
            this->controlForceMap->addForce(applied,force);
            this->activeControlForce = std::move(force);
            this->activeControlCommands = applied;
        }

        // Stage controller and conditioner state before the nonlinear solve.
        // RSolverGeneric's accepted-step transaction restores both on failure.
        this->candidateControllerState = decision.nextControllerState;
        this->trialControlRecord = {};
        this->trialControlRecord.acceptedTime = observation.acceptedTime;
        this->trialControlRecord.nextTime = observation.nextTime;
        this->trialControlRecord.sensors = observation.sensors;
        this->trialControlRecord.requestedCommands = decision.requestedCommands;
        this->trialControlRecord.conditionedCommands = conditioned.appliedCommands;
        this->trialControlRecord.controllerStateSha256 =
            QCryptographicHash::hash(decision.nextControllerState,
                                     QCryptographicHash::Sha256).toHex();

        std::vector<double> additionalNodalForces(nDofs,0.0);
        for (std::size_t dof=0;dof<nDofs;dof++)
            additionalNodalForces[dof] = this->activeControlForce[uint(dof)];
        this->structuralCableDynamicsResultReady = false;
        this->structuralCableDynamicsResult =
            RStructuralCableModelAdapter::solveDynamics(
                *this->pModel,previousDisplacements,previousVelocities,
                previousAccelerations,additionalNodalForces,changed);

        this->controlLedger->trial(*this->controlForceMap,this->activeControlCommands,
            RRVector(previousVelocities),RRVector(this->structuralCableDynamicsResult.velocities),
            physicalDt);
        this->trialControlRecord.appliedChannels = this->controlLedger->candidate();
    }
    else
    {
        const RStructuralDynamicsRestartData &restartData =
            this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
        if (restartData.available && restartData.controlAvailable)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable restart contains an accepted control ledger but no NodalForce control setup is enabled.");
        this->structuralCableDynamicsResultReady = false;
        this->structuralCableDynamicsResult =
            RStructuralCableModelAdapter::solveDynamics(
                *this->pModel,previousDisplacements,previousVelocities,
                previousAccelerations);
    }
    this->structuralCableDynamicsResultReady = true;
}

void RSolverStress::solve()
{
    try
    {
        if (this->isStructuralCableDynamicsMode())
        {
            this->solveStructuralCableDynamics();
        }
        else if (this->isStructuralCableMode())
        {
            RStructuralCableModelAdapter::validate(*this->pModel);
            this->structuralCableResultReady = false;
            this->structuralCableResult = RStructuralCableModelAdapter::solve(*this->pModel);
            this->structuralCableResultReady = true;
        }
        else if (this->isStructuralPiezoMode())
        {
            RStructuralPiezoModelAdapter::validate(*this->pModel);
            this->structuralPiezoResultReady = false;
            this->structuralPiezoResult = RStructuralPiezoModelAdapter::solve(*this->pModel);
            this->structuralPiezoResultReady = true;
        }
        else if (this->problemType == R_PROBLEM_STRESS_MODAL)
        {
            this->solveEigenValue();
        }
        else if (this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled())
        {
            this->solveStructuralDynamics();
        }
        else
        {
            this->solveStressStrain();
        }
    }
    catch (const RError &error)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve problem. %s", error.getMessage().toUtf8().constData());
    }
    catch (const std::exception &error)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve problem. %s",error.what());
    }
}

void RSolverStress::solveStressStrain()
{
    RLogger::info("Solving stress-strain problem.\n");

    try
    {
        RLogger::indent();
        const RMatrixSolverType solverType = this->nLinearConstraintEquations > 0
                                           ? RMatrixSolverConf::GMRES
                                           : RMatrixSolverConf::CG;
        const RMatrixSolverConf &matrixSolverConf = this->pModel->getMatrixSolverConf(solverType);
        RMatrixSolver matrixSolver(matrixSolverConf);
        RRVector candidate;
        matrixSolver.solve(this->A,this->b,candidate,
                           this->nLinearConstraintEquations > 0 ? R_MATRIX_PRECONDITIONER_NONE : R_MATRIX_PRECONDITIONER_JACOBI,
                           this->nLinearConstraintEquations > 0 ? 1 : 3);

        if (candidate.getNRows() != this->b.getNRows() || !finiteStressVector(candidate))
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Structural stress solve returned invalid or non-finite displacements.");
        }

        RRVector residual;
        RSparseMatrix::mlt(this->A,candidate,residual);
        if (residual.getNRows() != this->b.getNRows())
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Structural stress residual has an invalid dimension.");
        }
        for (uint i=0;i<residual.getNRows();i++)
        {
            residual[i] -= this->b[i];
        }

        const double residualNorm = RRVector::euclideanNorm(residual);
        const double rightHandSideNorm = RRVector::euclideanNorm(this->b);
        const double configuredTolerance = matrixSolverConf.getSolverCvgValue();
        if (!std::isfinite(residualNorm) || !std::isfinite(rightHandSideNorm) ||
            !std::isfinite(configuredTolerance) || configuredTolerance < 0.0)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Structural stress solve produced a non-finite residual or convergence tolerance.");
        }

        const double relativeTolerance = std::max(1.0e-8,configuredTolerance);
        const double allowedResidual = relativeTolerance*rightHandSideNorm+1.0e-12;
        if (!std::isfinite(allowedResidual) || residualNorm > allowedResidual)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Structural stress solve did not converge (residual %g, allowed %g).",
                         residualNorm,allowedResidual);
        }

        this->x = candidate;
        RLogger::unindent();
    }
    catch (const RError &error)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve stress-strain problem. %s", error.getMessage().toUtf8().constData());
    }

    this->setDisplacement(this->x);
}

void RSolverStress::validateDynamicBoundaryState()
{
    std::vector<double> current;
    current.reserve(8*this->pModel->getNNodes());
    for (uint node=0;node<this->pModel->getNNodes();node++)
    {
        const uint translational = this->nodeConstrainedDirections[node];
        const uint rotational = this->nodeHasRotationCondition[node]
            ? this->nodeConstrainedRotations[node]
            : (translational == 3 ? 7u : 0u);
        current.push_back(double(translational));
        current.push_back(double(rotational));
        for (uint component=0;component<3;component++)
        {
            current.push_back(component < translational
                              ? this->nodePrescribedDisplacement[node][component] : 0.0);
            current.push_back((rotational & (1u << component)) != 0u
                              ? this->nodePrescribedRotation[node][component] : 0.0);
        }
    }
    if (!this->dynamicBoundaryReferenceReady)
    {
        this->dynamicBoundaryReference = std::move(current);
        this->dynamicBoundaryReferenceReady = true;
        return;
    }
    if (current.size() != this->dynamicBoundaryReference.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural prescribed DOF layout changed during dynamics.");
    }
    for (uint index=0;index<current.size();index++)
    {
        const double tolerance = 1.0e-12*std::max({1.0,std::fabs(current[index]),
                                                    std::fabs(this->dynamicBoundaryReference[index])});
        if (std::fabs(current[index]-this->dynamicBoundaryReference[index]) > tolerance)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Moving prescribed displacement/rotation is not yet supported in second-order dynamics; base motion needs mass and damping coupling.");
        }
    }
}

RRVector RSolverStress::gatherNodalState(const RSolverCartesianVector<RRVector> &translation,
                                        const RSolverCartesianVector<RRVector> &rotation) const
{
    RRVector state(this->equationBook.getNEnabled(),0.0);
    for (uint node=0;node<this->pModel->getNNodes();node++)
    {
        RR3Vector globalTranslation(translation.x[node],translation.y[node],translation.z[node]);
        RR3Vector globalRotation(rotation.x[node],rotation.y[node],rotation.z[node]);
        for (uint component=0;component<3;component++)
        {
            uint equation = 0;
            if (this->equationBook.getValue(node,translationDof(component),equation))
            {
                double value = globalTranslation[component];
                if (this->localRotations[node].isActive())
                {
                    value = 0.0;
                    for (uint global=0;global<3;global++)
                    {
                        value += this->localRotations[node].getR()[global][component]*globalTranslation[global];
                    }
                }
                state[equation] = value;
            }
            const RDofType rotationDof = component == 0 ? RDofType::RX
                                          : component == 1 ? RDofType::RY : RDofType::RZ;
            if (this->equationBook.getValue(node,rotationDof,equation))
            {
                double value = globalRotation[component];
                if (this->rotationLocalRotations[node].isActive())
                {
                    value = 0.0;
                    for (uint global=0;global<3;global++)
                    {
                        value += this->rotationLocalRotations[node].getR()[global][component]*globalRotation[global];
                    }
                }
                state[equation] = value;
            }
        }
    }
    return state;
}

void RSolverStress::scatterNodalState(const RRVector &physical,
                                     RSolverCartesianVector<RRVector> &translation,
                                     RSolverCartesianVector<RRVector> &rotation) const
{
    const uint nNodes = this->pModel->getNNodes();
    translation.x.resize(nNodes,0.0);
    translation.y.resize(nNodes,0.0);
    translation.z.resize(nNodes,0.0);
    rotation.x.resize(nNodes,0.0);
    rotation.y.resize(nNodes,0.0);
    rotation.z.resize(nNodes,0.0);
    for (uint node=0;node<nNodes;node++)
    {
        RR3Vector localTranslation(0.0,0.0,0.0);
        RR3Vector localRotation(0.0,0.0,0.0);
        for (uint component=0;component<3;component++)
        {
            uint equation = 0;
            if (this->equationBook.getValue(node,translationDof(component),equation))
            {
                localTranslation[component] = physical[equation];
            }
            const RDofType rotationDof = component == 0 ? RDofType::RX
                                          : component == 1 ? RDofType::RY : RDofType::RZ;
            if (this->equationBook.getValue(node,rotationDof,equation))
            {
                localRotation[component] = physical[equation];
            }
        }
        if (this->localRotations[node].isActive())
        {
            this->localRotations[node].rotateResultsVector(localTranslation);
        }
        if (this->rotationLocalRotations[node].isActive())
        {
            this->rotationLocalRotations[node].rotateResultsVector(localRotation);
        }
        translation.x[node] = localTranslation[0];
        translation.y[node] = localTranslation[1];
        translation.z[node] = localTranslation[2];
        rotation.x[node] = localRotation[0];
        rotation.y[node] = localRotation[1];
        rotation.z[node] = localRotation[2];
    }
}

std::vector<RStructuralDofRef> RSolverStress::orderedPhysicalDofs() const
{
    std::vector<RStructuralDofRef> ordered(this->nPhysicalEquations);
    uint visited = 0;
    for (uint node=0;node<this->equationBook.getNNodes();node++)
    {
        for (RDofType dof : this->equationBook.getEnabledDofs(node))
        {
            uint equation = 0;
            if (this->equationBook.getValue(node,dof,equation) && equation < ordered.size())
            {
                ordered[equation] = {node,dof};
                visited++;
            }
        }
    }
    if (visited != ordered.size())
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Structural equation book is inconsistent with physical DOF count.");
    }
    return ordered;
}

void RSolverStress::prepareControlMapping()
{
    if (this->controlForceMap) return;
    if (this->isStructuralCableDynamicsMode())
    {
        this->prepareCableControlMapping();
        return;
    }
    const auto &catalog = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    const auto &actuators = catalog.getActuators();
    const auto &sensors = catalog.getSensors();
    if (actuators.empty())
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Enabled structural control needs at least one nodal actuator.");
    const uint n = this->A.getNRows();
    auto equationDirection = [this](uint node, const RR3Vector &global,
                                     const QString &id) -> RR3Vector
    {
        if (node >= this->pModel->getNNodes())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' references a missing node.",id.toUtf8().constData());
        if (!this->localRotations[node].isActive()) return global;
        const RRMatrix &rotation = this->localRotations[node].getR();
        RR3Vector local(0.0,0.0,0.0);
        for (uint component=0;component<3;component++)
            for (uint globalAxis=0;globalAxis<3;globalAxis++)
                local[component] += rotation[globalAxis][component]*global[globalAxis];
        return local;
    };
    auto weightsForNode = [this,n,&equationDirection](uint node, const RR3Vector &globalDirection,
                                    const QString &id) -> std::vector<std::pair<uint,double>>
    {
        if (node >= this->pModel->getNNodes())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' references a missing node.",id.toUtf8().constData());
        const RR3Vector direction = equationDirection(node,globalDirection,id);
        std::map<uint,double> nonzero;
        for (uint component=0;component<3;component++)
        {
            const double value = direction[component];
            if (std::fabs(value) <= 1.0e-14) continue;
            uint equation = 0;
            if (!this->equationBook.getValue(node,translationDof(component),equation) ||
                equation >= this->nPhysicalEquations)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Control target '%s' includes a held or absent translation DOF.",
                             id.toUtf8().constData());
            if (this->constraintNullspace.getNRows() == 0)
                nonzero[equation] += value;
            else
                for (uint reduced=0;reduced<n;reduced++)
                    if (std::fabs(this->constraintNullspace[equation][reduced]) > 1.0e-14)
                        nonzero[reduced] += value*this->constraintNullspace[equation][reduced];
        }
        std::vector<std::pair<uint,double>> result;
        for (const auto &entry : nonzero)
            if (std::fabs(entry.second) > 1.0e-14)
                result.push_back(entry);
        if (result.empty())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' has no independent reduced DOF.",
                         id.toUtf8().constData());
        return result;
    };
    auto unitDirection = [](const RR3Vector &input, const QString &id) -> RR3Vector
    {
        const double norm = input.length();
        if (!std::isfinite(norm) || norm <= 1.0e-12)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' needs a finite nonzero global direction.",
                         id.toUtf8().constData());
        return RR3Vector(input[0]/norm,input[1]/norm,input[2]/norm);
    };
    std::vector<RStructuralControlChannelMap> channels;
    channels.reserve(actuators.size());
    for (const auto &actuator : actuators)
    {
        if (actuator.type != RStructuralActuatorType::NodalForce ||
            actuator.target.scope != RStructuralTargetScope::Node ||
            actuator.frame != RStructuralFrame::Global)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control actuator '%s' needs point-node global NodalForce v1.",
                         actuator.id.toUtf8().constData());
        RR3Vector direction = actuator.direction;
        switch (actuator.variable)
        {
            case R_VARIABLE_FORCE_X: direction = RR3Vector(1.0,0.0,0.0); break;
            case R_VARIABLE_FORCE_Y: direction = RR3Vector(0.0,1.0,0.0); break;
            case R_VARIABLE_FORCE_Z: direction = RR3Vector(0.0,0.0,1.0); break;
            case R_VARIABLE_FORCE: direction = unitDirection(direction,actuator.id); break;
            default:
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Control actuator '%s' has unsupported command variable.",
                             actuator.id.toUtf8().constData());
        }
        channels.push_back({actuator.id,weightsForNode(actuator.target.nodeId,direction,actuator.id)});
    }
    this->controlSensors.clear();
    this->controlSensors.reserve(sensors.size());
    for (const auto &sensor : sensors)
    {
        if (sensor.target.scope != RStructuralTargetScope::Node ||
            sensor.frame != RStructuralFrame::Global ||
            (sensor.quantity != RStructuralSensorQuantity::Displacement &&
             sensor.quantity != RStructuralSensorQuantity::Velocity &&
             sensor.quantity != RStructuralSensorQuantity::Acceleration))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control sensor '%s' needs point-node global displacement, velocity or acceleration v1.",
                         sensor.id.toUtf8().constData());
        const RR3Vector direction = unitDirection(sensor.direction,sensor.id);
        const RR3Vector localDirection = equationDirection(sensor.target.nodeId,direction,sensor.id);
        double particular = 0.0;
        if (this->constraintNullspace.getNRows() != 0 &&
            sensor.quantity == RStructuralSensorQuantity::Displacement)
        {
            for (uint component=0;component<3;component++)
            {
                if (std::fabs(localDirection[component]) <= 1.0e-14) continue;
                uint equation = 0;
                if (!this->equationBook.getValue(sensor.target.nodeId,
                                                 translationDof(component),equation))
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Control sensor '%s' includes a held DOF.",
                                 sensor.id.toUtf8().constData());
                particular += localDirection[component]*this->constraintParticular[equation];
            }
        }
        this->controlSensors.push_back({sensor.id,sensor.quantity,sensor.unit,
            weightsForNode(sensor.target.nodeId,direction,sensor.id),particular});
    }
    this->controlForceMap.reset(new RStructuralControlForceMap(
        n,this->scales.findScaleFactor(R_VARIABLE_FORCE_X),
        this->scales.findScaleFactor(R_VARIABLE_VELOCITY_X),std::move(channels)));
    this->controlLedger.reset(new RStructuralControlLedger(*this->controlForceMap));
    std::vector<RStructuralControlConditioningPolicy> policies;
    policies.reserve(actuators.size());
    std::map<QString,size_t> positions;
    for (const auto &channel : this->controlForceMap->getChannels())
    {
        positions.emplace(channel.actuatorId,policies.size());
        policies.push_back({channel.actuatorId,0,1,0.0});
    }
    const auto &controlSetup = this->pModel->getProblemSetup().getStructuralControlSetup();
    for (const auto &configured : controlSetup.channelConditioning)
    {
        const auto position = positions.find(configured.actuatorId);
        if (position == positions.end())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control conditioner references unknown actuator '%s'.",
                         configured.actuatorId.toUtf8().constData());
        policies[position->second] = {configured.actuatorId,configured.delaySteps,
                                      configured.holdSteps,configured.noiseAmplitude};
    }
    this->controlConditioner.reset(new RStructuralControlConditioner(
        controlSetup.deterministicSeed,std::move(policies)));
    this->acceptedConditionerState = this->controlConditioner->initialState();
    this->candidateConditionerState = this->acceptedConditionerState;
    this->activeControlCommands.clear();
    for (const auto &channel : this->controlForceMap->getChannels())
        this->activeControlCommands.push_back({channel.actuatorId,0.0});
    this->activeControlForce.resize(n,0.0);
    this->cachedControlCatalogFingerprint = this->controlCatalogFingerprint();
    this->cachedControlFingerprint = this->controlMappingFingerprint();
}

void RSolverStress::prepareCableControlMapping()
{
    if (!this->isStructuralCableDynamicsMode())
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Full-nodal control mapping is available only for transient TRUSS1 cable mode.");
    const auto &catalog = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    const auto &actuators = catalog.getActuators();
    const auto &sensors = catalog.getSensors();
    if (actuators.empty())
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Enabled cable control needs at least one point-node global NodalForce actuator.");

    const uint n = 3*this->pModel->getNNodes();
    const std::vector<bool> prescribed =
        RStructuralCableModelAdapter::prescribedDofMask(*this->pModel);
    if (prescribed.size() != n)
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Cable control map received a mismatched physical DOF mask.");

    auto unitDirection = [](const RR3Vector &input,const QString &id) -> RR3Vector
    {
        const double norm = input.length();
        if (!std::isfinite(norm) || norm <= 1.0e-12)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' needs a finite nonzero global direction.",
                         id.toUtf8().constData());
        return RR3Vector(input[0]/norm,input[1]/norm,input[2]/norm);
    };
    auto directWeights = [&prescribed,n](uint node,const RR3Vector &direction,
                                         const QString &id,bool requireFreeDofs)
        -> std::vector<std::pair<uint,double>>
    {
        if (node >= n/3)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' references a missing cable node.",
                         id.toUtf8().constData());
        std::vector<std::pair<uint,double>> weights;
        for (uint axis=0;axis<3;axis++)
        {
            const double value = direction[axis];
            if (std::fabs(value) <= 1.0e-14) continue;
            const uint dof = 3*node+axis;
            if (requireFreeDofs && prescribed[dof])
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Control target '%s' includes a fixed cable translation DOF.",
                             id.toUtf8().constData());
            weights.emplace_back(dof,value);
        }
        if (weights.empty())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control target '%s' has no free global cable translation DOF.",
                         id.toUtf8().constData());
        return weights;
    };

    std::vector<RStructuralControlChannelMap> channels;
    channels.reserve(actuators.size());
    for (const RStructuralActuator &actuator : actuators)
    {
        if (actuator.type != RStructuralActuatorType::NodalForce ||
            actuator.target.scope != RStructuralTargetScope::Node ||
            actuator.frame != RStructuralFrame::Global)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control actuator '%s' requires point-node global NodalForce v1.",
                         actuator.id.toUtf8().constData());
        }
        RR3Vector direction = actuator.direction;
        switch (actuator.variable)
        {
            case R_VARIABLE_FORCE_X: direction = RR3Vector(1.0,0.0,0.0); break;
            case R_VARIABLE_FORCE_Y: direction = RR3Vector(0.0,1.0,0.0); break;
            case R_VARIABLE_FORCE_Z: direction = RR3Vector(0.0,0.0,1.0); break;
            case R_VARIABLE_FORCE: direction = unitDirection(direction,actuator.id); break;
            default:
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable control actuator '%s' command must be a physical NodalForce component.",
                             actuator.id.toUtf8().constData());
        }
        channels.push_back({actuator.id,
            directWeights(actuator.target.nodeId,direction,actuator.id,true)});
    }

    this->controlSensors.clear();
    this->controlSensors.reserve(sensors.size());
    for (const RStructuralSensor &sensor : sensors)
    {
        if (sensor.target.scope != RStructuralTargetScope::Node ||
            sensor.aggregation != RStructuralSensorAggregation::Point ||
            sensor.frame != RStructuralFrame::Global ||
            (sensor.quantity != RStructuralSensorQuantity::Displacement &&
             sensor.quantity != RStructuralSensorQuantity::Velocity &&
             sensor.quantity != RStructuralSensorQuantity::Acceleration))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control sensor '%s' requires point-node global displacement, velocity or acceleration v1.",
                         sensor.id.toUtf8().constData());
        }
        const RR3Vector direction = unitDirection(sensor.direction,sensor.id);
        this->controlSensors.push_back({sensor.id,sensor.quantity,sensor.unit,
            directWeights(sensor.target.nodeId,direction,sensor.id,false),0.0});
    }

    this->controlForceMap.reset(new RStructuralControlForceMap(
        n,this->scales.findScaleFactor(R_VARIABLE_FORCE_X),
        this->scales.findScaleFactor(R_VARIABLE_VELOCITY_X),std::move(channels)));
    this->controlLedger.reset(new RStructuralControlLedger(*this->controlForceMap));
    std::vector<RStructuralControlConditioningPolicy> policies;
    policies.reserve(actuators.size());
    std::map<QString,size_t> positions;
    for (const auto &channel : this->controlForceMap->getChannels())
    {
        positions.emplace(channel.actuatorId,policies.size());
        policies.push_back({channel.actuatorId,0,1,0.0});
    }
    const RStructuralControlSetup &controlSetup =
        this->pModel->getProblemSetup().getStructuralControlSetup();
    for (const auto &configured : controlSetup.channelConditioning)
    {
        const auto position = positions.find(configured.actuatorId);
        if (position == positions.end())
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control conditioner references unknown actuator '%s'.",
                         configured.actuatorId.toUtf8().constData());
        policies[position->second] = {configured.actuatorId,configured.delaySteps,
                                      configured.holdSteps,configured.noiseAmplitude};
    }
    this->controlConditioner.reset(new RStructuralControlConditioner(
        controlSetup.deterministicSeed,std::move(policies)));
    this->acceptedConditionerState = this->controlConditioner->initialState();
    this->candidateConditionerState = this->acceptedConditionerState;
    this->acceptedControllerState.clear();
    this->candidateControllerState.clear();
    this->activeControlCommands.clear();
    for (const auto &channel : this->controlForceMap->getChannels())
        this->activeControlCommands.push_back({channel.actuatorId,0.0});
    this->activeControlForce.resize(n,0.0);
    this->cachedControlCatalogFingerprint = this->controlCatalogFingerprint();
    this->cachedControlFingerprint = this->controlMappingFingerprint();

    if (this->pModel->getProblemSetup().getRestart())
    {
        const RStructuralDynamicsRestartData &restartData =
            this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
        if (!restartData.controlAvailable ||
            restartData.controlProtocolVersion != controlSetup.controlProtocolVersion ||
            restartData.controlDeterministicSeed != controlSetup.deterministicSeed ||
            restartData.controlFingerprint != this->cachedControlFingerprint ||
            restartData.controlChannelIds.size() != this->activeControlCommands.size())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control restart protocol, seed, actuator map or ledger does not match.");
        }
        std::vector<RStructuralControlChannelSample> samples;
        samples.reserve(this->activeControlCommands.size());
        for (size_t i=0;i<this->activeControlCommands.size();i++)
        {
            if (restartData.controlChannelIds[i] != this->activeControlCommands[i].actuatorId)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Cable control restart actuator order changed.");
            samples.push_back({restartData.controlChannelIds[i],
                               restartData.controlLastCommands[i],
                               restartData.controlPowers[i],
                               restartData.cumulativeControlWork[i]});
            this->activeControlCommands[i].value = restartData.controlLastCommands[i];
        }
        this->controlLedger->restore(samples);
        const RStructuralControlStateEnvelope envelope =
            RStructuralControlStateEnvelope::deserialize(restartData.controlControllerState);
        this->acceptedControllerState = envelope.controllerState;
        this->candidateControllerState = this->acceptedControllerState;
        this->acceptedConditionerState = this->controlConditioner->deserialize(
            envelope.conditionerState);
        this->candidateConditionerState = this->acceptedConditionerState;
        if (this->acceptedConditionerState.completedSteps != restartData.completedSteps)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Cable control restart conditioner step count differs from its accepted checkpoint.");
        this->activeControlForce.fill(0.0);
        this->controlForceMap->addForce(this->activeControlCommands,
                                        this->activeControlForce);
    }
}

QString RSolverStress::controlCatalogFingerprint() const
{
    QCryptographicHash hash(QCryptographicHash::Sha256);
    hash.addData("range-fea/control-catalog/v1;");
    auto string = [&hash](const QString &value)
    {
        const QByteArray bytes = value.toUtf8();
        hash.addData(QByteArray::number(bytes.size())); hash.addData(":");
        hash.addData(bytes);
    };
    auto number = [&hash](quint64 value)
    {
        hash.addData(QByteArray::number(value)); hash.addData(";");
    };
    auto real = [&hash](double value)
    {
        hash.addData(QByteArray::number(value,'g',17)); hash.addData(";");
    };
    const auto &setup = this->pModel->getProblemSetup().getStructuralControlSetup();
    number(setup.enabled ? 1 : 0);
    number(setup.controlProtocolVersion);
    string(setup.controllerId);
    number(setup.deterministicSeed);
    number(setup.channelConditioning.size());
    for (const auto &policy : setup.channelConditioning)
    {
        string(policy.actuatorId);
        number(policy.delaySteps);
        number(policy.holdSteps);
        real(policy.noiseAmplitude);
    }
    const auto &catalog = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    number(catalog.getActuators().size());
    for (const auto &actuator : catalog.getActuators())
    {
        string(actuator.id);
        number(uint(actuator.type));
        number(uint(actuator.variable));
        number(uint(actuator.frame));
        number(uint(actuator.target.scope));
        number(actuator.target.nodeId);
        string(actuator.target.setName);
        number(uint(actuator.target.entityType));
        string(actuator.unit);
        for (uint axis=0;axis<3;axis++) real(actuator.direction[axis]);
        real(actuator.minimum); real(actuator.maximum); real(actuator.rateLimit);
    }
    number(catalog.getSensors().size());
    for (const auto &sensor : catalog.getSensors())
    {
        string(sensor.id);
        number(uint(sensor.quantity));
        number(uint(sensor.aggregation));
        number(uint(sensor.frame));
        number(uint(sensor.target.scope));
        number(sensor.target.nodeId);
        string(sensor.target.setName);
        number(uint(sensor.target.entityType));
        string(sensor.unit);
        for (uint axis=0;axis<3;axis++) real(sensor.direction[axis]);
    }
    return QString::fromLatin1(hash.result().toHex());
}

QString RSolverStress::controlMappingFingerprint() const
{
    QCryptographicHash hash(QCryptographicHash::Sha256);
    hash.addData("range-fea/control-map/v1;");
    const auto &setup = this->pModel->getProblemSetup().getStructuralControlSetup();
    hash.addData(setup.controllerId.toUtf8()); hash.addData(";");
    hash.addData(this->controllerRevision.toUtf8()); hash.addData(";");
    hash.addData(this->controlConditioner->configurationSha256().toHex()); hash.addData(";");
    hash.addData(QByteArray::number(setup.deterministicSeed)); hash.addData(";");
    hash.addData((this->isStructuralCableDynamicsMode()
        ? RStructuralCableModelAdapter::computeDynamicsSystemFingerprint(*this->pModel)
        : this->dynamicSystemFingerprint()).toUtf8()); hash.addData(";");
    if (this->isStructuralCableDynamicsMode())
    {
        // The exact catalog digest also binds fields which are intentionally
        // not representable in the NodalForce map (for example named targets).
        hash.addData(this->controlCatalogFingerprint().toUtf8()); hash.addData(";");
    }
    const auto &actuators = this->pModel->getProblemSetup().getStructuralInteractionSetup().getActuators();
    for (size_t i=0;i<actuators.size();i++)
    {
        const auto &a = actuators[i];
        hash.addData(a.id.toUtf8()); hash.addData(";");
        hash.addData(QByteArray::number(int(a.type))); hash.addData(";");
        hash.addData(QByteArray::number(int(a.variable))); hash.addData(";");
        hash.addData(QByteArray::number(int(a.frame))); hash.addData(";");
        hash.addData(QByteArray::number(a.target.nodeId)); hash.addData(";");
        hash.addData(a.unit.toUtf8()); hash.addData(";");
        for (uint axis=0;axis<3;axis++)
        {
            hash.addData(QByteArray::number(a.direction[axis],'g',17)); hash.addData(";");
        }
        for (double limit : {a.minimum,a.maximum,a.rateLimit})
        {
            hash.addData(QByteArray::number(limit,'g',17)); hash.addData(";");
        }
        for (const auto &weight : this->controlForceMap->getChannels()[i].reducedWeights)
        {
            hash.addData(QByteArray::number(weight.first)); hash.addData(":");
            hash.addData(QByteArray::number(weight.second,'g',17)); hash.addData(";");
        }
    }
    const auto &sensorCatalog = this->pModel->getProblemSetup().getStructuralInteractionSetup().getSensors();
    for (size_t i=0;i<this->controlSensors.size();i++)
    {
        const auto &sensor = this->controlSensors[i];
        const auto &definition = sensorCatalog[i];
        hash.addData(sensor.id.toUtf8()); hash.addData(":");
        hash.addData(QByteArray::number(int(sensor.quantity))); hash.addData(":");
        hash.addData(QByteArray::number(int(definition.frame))); hash.addData(":");
        hash.addData(QByteArray::number(definition.target.nodeId)); hash.addData(":");
        hash.addData(sensor.unit.toUtf8()); hash.addData(":");
        for (uint axis=0;axis<3;axis++)
        {
            hash.addData(QByteArray::number(definition.direction[axis],'g',17)); hash.addData(";");
        }
        hash.addData(QByteArray::number(sensor.displacementOffset,'g',17)); hash.addData(";");
        for (const auto &weight : sensor.weights)
        {
            hash.addData(QByteArray::number(weight.first)); hash.addData(":");
            hash.addData(QByteArray::number(weight.second,'g',17)); hash.addData(";");
        }
    }
    return QString::fromLatin1(hash.result().toHex());
}

std::vector<RStructuralControlSensorSample> RSolverStress::sampleControlSensors(
    const RStructuralDynamicsState &state) const
{
    std::vector<RStructuralControlSensorSample> samples;
    samples.reserve(this->controlSensors.size());
    for (const auto &sensor : this->controlSensors)
    {
        const RRVector *source = nullptr;
        RVariableType scaleType = R_VARIABLE_DISPLACEMENT_X;
        if (sensor.quantity == RStructuralSensorQuantity::Displacement)
            source = &state.q;
        else if (sensor.quantity == RStructuralSensorQuantity::Velocity)
        {
            source = &state.v; scaleType = R_VARIABLE_VELOCITY_X;
        }
        else
        {
            source = &state.a; scaleType = R_VARIABLE_ACCELERATION_X;
        }
        double value = 0.0;
        for (const auto &weight : sensor.weights)
            value += weight.second*(*source)[weight.first];
        if (sensor.quantity == RStructuralSensorQuantity::Displacement)
            value += sensor.displacementOffset;
        value /= this->scales.findScaleFactor(scaleType);
        samples.push_back({sensor.id,value,sensor.unit});
    }
    return samples;
}

QString RSolverStress::motionDefinitionFingerprint() const
{
    QCryptographicHash hash(QCryptographicHash::Sha256);
    const auto &motions = this->pModel->getProblemSetup().getStructuralInteractionSetup()
        .getPrescribedMotions();
    hash.addData("range-fea/dynamics/motion/v1;");
    auto appendInteger = [&hash](quint64 value)
    {
        hash.addData(QByteArray::number(value));
        hash.addData(";",1);
    };
    auto appendReal = [&hash](double value)
    {
        hash.addData(QByteArray::number(value == 0.0 ? 0.0 : value,'g',17));
        hash.addData(";",1);
    };
    appendInteger(motions.size());
    for (const RStructuralPrescribedMotion &motion : motions)
    {
        const QByteArray id = motion.id.toUtf8();
        appendInteger(id.size());
        hash.addData(id);
        appendInteger(motion.target.nodeId);
        appendInteger(static_cast<uint>(motion.target.dof));
        appendInteger(motion.knots.size());
        for (const RStructuralMotionKnot &knot : motion.knots)
        {
            appendReal(knot.time);
            appendReal(knot.displacement);
            appendReal(knot.velocity);
            appendReal(knot.acceleration);
        }
    }
    return QString::fromLatin1(hash.result().toHex());
}

void RSolverStress::preparePrescribedMotions()
{
    const auto &interaction = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    const QString fingerprint = this->motionDefinitionFingerprint();
    if (!this->dynamicMotionDefinitionReference.isEmpty() &&
        this->dynamicMotionDefinitionReference != fingerprint)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Prescribed support motion changed during fixed-matrix dynamics.");
    }
    this->dynamicMotionDefinitionReference = fingerprint;
    this->prescribedMotions = interaction.getPrescribedMotions();
    this->motionIndexByNodeDof.assign(3*this->pModel->getNNodes(),-1);
    if (this->prescribedMotions.empty())
    {
        return;
    }
    if (!interaction.getLinearConstraints().empty() || !interaction.getCouplings().empty())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Prescribed support motion with MPC/coupling is not supported by the first FE implementation.");
    }
    const RTimeSolver &time = this->pModel->getTimeSolver();
    // The t0 force prepass temporarily replaces times[0] with t0; asking
    // getPreviousTime() there would fabricate t0-dt outside the knot domain.
    const double firstPhysicalTime = (this->preparingInitialForce
        ? time.getCurrentTime() : time.getPreviousTime())/this->scales.getSecond();
    const double lastPhysicalTime = time.getTimes().back()/this->scales.getSecond();
    for (uint index=0;index<this->prescribedMotions.size();index++)
    {
        const RStructuralPrescribedMotion &motion = this->prescribedMotions[index];
        uint component = 0;
        if (!findTranslationComponent(motion.target.dof,component))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Prescribed support motion '%s' must target a global translational DOF.",
                         motion.id.toUtf8().constData());
        }
        const uint node = motion.target.nodeId;
        if (this->localRotations[node].isActive() ||
            this->rotationLocalRotations[node].isActive())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Prescribed support motion '%s' cannot use a local node frame.",
                         motion.id.toUtf8().constData());
        }
        uint equation = 0;
        if (this->equationBook.getValue(node,motion.target.dof,equation) ||
            component >= this->nodeConstrainedDirections[node])
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Prescribed support motion '%s' requires an existing held displacement DOF.",
                         motion.id.toUtf8().constData());
        }
        const double baseline = this->nodePrescribedDisplacement[node][component];
        if (!std::isfinite(baseline) || baseline != 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Prescribed support motion '%s' requires a finite, exactly zero baseline displacement BC.",
                         motion.id.toUtf8().constData());
        }
        // The ordinary BC supplies a held-DOF mask only. Reject a delayed
        // nonzero BC knot now, before accepting any time step from two sources.
        const RVariableType displacementTypes[3] =
            {R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
        bool hasGlobalMask = false;
        for (uint groupId=0;groupId<this->pModel->getNElementGroups();groupId++)
        {
            const RElementGroup *group = this->pModel->getElementGroupPtr(groupId);
            if (!group)
            {
                continue;
            }
            bool containsNode = false;
            for (uint elementPosition=0;elementPosition<group->size() && !containsNode;elementPosition++)
            {
                const RElement &element = this->pModel->getElement(group->get(elementPosition));
                for (uint localNode=0;localNode<element.size();localNode++)
                {
                    containsNode = containsNode || element.getNodeId(localNode) == node;
                }
            }
            if (!containsNode)
            {
                continue;
            }
            for (uint bcIndex=0;bcIndex<group->getNBoundaryConditions();bcIndex++)
            {
                const RBoundaryCondition &bc = group->getBoundaryCondition(bcIndex);
                if (bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT ||
                    !RSolverStress::isComponentEnabled(bc,displacementTypes[component]))
                {
                    continue;
                }
                const RConditionComponent &mask =
                    bc.getComponent(bc.findComponentPosition(displacementTypes[component]));
                if (bc.getHasLocalDirection())
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Prescribed support motion '%s' cannot use a local displacement BC.",
                                 motion.id.toUtf8().constData());
                }
                hasGlobalMask = true;
                for (uint valueIndex=0;valueIndex<mask.size();valueIndex++)
                {
                    const double value = mask.getValue(valueIndex);
                    if (!std::isfinite(value) || value != 0.0)
                    {
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                     "Prescribed support motion '%s' requires every baseline BC knot to be finite and exactly zero.",
                                     motion.id.toUtf8().constData());
                    }
                }
            }
        }
        if (!hasGlobalMask)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Prescribed support motion '%s' cannot use a local-frame BC; an explicit global displacement BC mask is required.",
                         motion.id.toUtf8().constData());
        }
        // Use the model evaluator itself for the endpoint contract. Its
        // 16-ulp rounding tolerance is narrower than an absolute 1e-12 s;
        // a looser preflight could accept early steps then fail at the end.
        motion.evaluate(firstPhysicalTime);
        motion.evaluate(lastPhysicalTime);
        this->motionIndexByNodeDof[3*node+component] = int(index);
    }
}

void RSolverStress::prepareSupportResponse()
{
    this->heldDofs.clear();
    this->pointSupportMask.clear();
    this->heldIndexByNodeDof.assign(6*this->pModel->getNNodes(),-1);
    if (!this->supportResponseObserver) return;
    const auto &interaction = this->pModel->getProblemSetup().getStructuralInteractionSetup();
    if (!interaction.getLinearConstraints().empty() || !interaction.getCouplings().empty())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Support reaction sampling with MPC/coupling is not supported.");
    }
    if (!this->pModel->getProblemSetup().getStructuralDynamicsSetup()
             .getModalDampingRatios().empty())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Support reaction sampling with modal damping has no defined constrained damping rows.");
    }
    for (uint elementId=0;elementId<this->pModel->getNElements();elementId++)
    {
        if (!this->computableElements[elementId]) continue;
        const auto &dofs = RElementDescriptor::get(
            this->pModel->getElement(elementId).getType()).getDofsPerNode();
        if (dofs.contains(RDofType::RX) || dofs.contains(RDofType::RY) ||
            dofs.contains(RDofType::RZ))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Support reaction sampling of frame rotational DOFs/moments is not supported.");
        }
    }
    auto appendHeld = [this](uint node, RDofType dof)
    {
        const uint key = 6*node+uint(dof);
        this->heldIndexByNodeDof[key] = int(this->heldDofs.size());
        this->heldDofs.push_back({node,dof});
        this->pointSupportMask.push_back(false);
    };
    for (uint node=0;node<this->pModel->getNNodes();node++)
    {
        for (uint component=0;component<this->nodeConstrainedDirections[node];component++)
        {
            appendHeld(node,translationDof(component));
        }
        const uint rotationalMask = this->nodeHasRotationCondition[node]
            ? this->nodeConstrainedRotations[node]
            : (this->nodeConstrainedDirections[node] == 3 ? 7u : 0u);
        for (uint component=0;component<3;component++)
        {
            if ((rotationalMask & (1u << component)) != 0u)
            {
                appendHeld(node,component == 0 ? RDofType::RX :
                                component == 1 ? RDofType::RY : RDofType::RZ);
            }
        }
    }
    const RVariableType components[3] =
        {R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for (uint point=0;point<this->pModel->getNPoints();point++)
    {
        const RPoint &group = this->pModel->getPoint(point);
        for (uint bcIndex=0;bcIndex<group.getNBoundaryConditions();bcIndex++)
        {
            const RBoundaryCondition &bc = group.getBoundaryCondition(bcIndex);
            if (bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT || bc.getHasLocalDirection())
                continue;
            for (uint position=0;position<group.size();position++)
            {
                const RElement &element = this->pModel->getElement(group.get(position));
                const uint node = element.getNodeId(0);
                if (this->localRotations[node].isActive())
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Support reaction sampling requires a global point-node constraint frame.");
                }
                for (uint component=0;component<3;component++)
                {
                    if (!RSolverStress::isComponentEnabled(bc,components[component])) continue;
                    const int index = this->heldIndexByNodeDof[6*node+component];
                    if (index >= 0) this->pointSupportMask[uint(index)] = true;
                }
            }
        }
    }
}

RStructuralMotionState RSolverStress::scaledMotionState(
    const RStructuralPrescribedMotion &motion, double time) const
{
    RStructuralMotionState state = motion.evaluate(time/this->scales.getSecond());
    state.displacement *= this->scales.findScaleFactor(R_VARIABLE_DISPLACEMENT_X);
    state.velocity *= this->scales.findScaleFactor(R_VARIABLE_VELOCITY_X);
    state.acceleration *= this->scales.findScaleFactor(R_VARIABLE_ACCELERATION_X);
    return state;
}

void RSolverStress::motionEffectiveForce(double time, RRVector &force) const
{
    force.resize(this->A.getNRows());
    force.fill(0.0);
    if (this->prescribedMotions.empty())
    {
        return;
    }
    std::vector<RStructuralMotionState> states;
    states.reserve(this->prescribedMotions.size());
    for (const RStructuralPrescribedMotion &motion : this->prescribedMotions)
    {
        states.push_back(this->scaledMotionState(motion,time));
    }
    const RStructuralDynamicsSetup &setup = this->pModel->getProblemSetup()
        .getStructuralDynamicsSetup();
    // M_ff qdd + C_ff qd + K_ff q = f_f - M_fc ucdd - C_fc ucd - K_fc uc.
    // These are sparse free/support blocks; no dense identity transform is built.
    for (uint row=0;row<force.size();row++)
    {
        for (uint p=0;p<this->motionStiffness.getNColumns(row);p++)
        {
            const uint column = this->motionStiffness.getVector(row).getIndex(p);
            const double value = this->motionStiffness.getValue(row,p);
            force[row] -= value*(states[column].displacement+
                                 setup.getRayleighBeta()*states[column].velocity);
        }
        for (uint p=0;p<this->motionMass.getNColumns(row);p++)
        {
            const uint column = this->motionMass.getVector(row).getIndex(p);
            const double value = this->motionMass.getValue(row,p);
            force[row] -= value*(states[column].acceleration+
                                 setup.getRayleighAlpha()*states[column].velocity);
        }
    }
}

void RSolverStress::setPrescribedMotionResults(double time)
{
    for (const RStructuralPrescribedMotion &motion : this->prescribedMotions)
    {
        uint component = 0;
        findTranslationComponent(motion.target.dof,component);
        const RStructuralMotionState state = this->scaledMotionState(motion,time);
        const uint node = motion.target.nodeId;
        RRVector *displacement[3] = {&this->nodeDisplacement.x,&this->nodeDisplacement.y,&this->nodeDisplacement.z};
        RRVector *velocity[3] = {&this->nodeVelocity.x,&this->nodeVelocity.y,&this->nodeVelocity.z};
        RRVector *acceleration[3] = {&this->nodeAcceleration.x,&this->nodeAcceleration.y,&this->nodeAcceleration.z};
        (*displacement[component])[node] = state.displacement;
        (*velocity[component])[node] = state.velocity;
        (*acceleration[component])[node] = state.acceleration;
    }
}

std::vector<RStructuralSupportSample> RSolverStress::computeSupportSamples(
    const RStructuralDynamicsState &state, const RRVector &heldForce,
    bool initialOperator) const
{
    const RSparseMatrix &stiffness = initialOperator
        ? this->supportStiffnessReference : this->supportStiffness;
    const RSparseMatrix &mass = initialOperator
        ? this->supportMassReference : this->supportMass;
    const uint nFree = this->equationBook.getNEnabled();
    if (state.q.size() != nFree || state.v.size() != nFree || state.a.size() != nFree ||
        heldForce.size() != this->heldDofs.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Support response state/equation dimensions do not agree.");
    }
    RRVector q(nFree+uint(this->heldDofs.size()),0.0);
    RRVector v(q.size(),0.0);
    RRVector a(q.size(),0.0);
    for (uint index=0;index<nFree;index++)
    {
        q[index] = state.q[index];
        v[index] = state.v[index];
        a[index] = state.a[index];
    }
    for (uint index=0;index<this->heldDofs.size();index++)
    {
        const RStructuralDofRef &target = this->heldDofs[index];
        uint component = 0;
        if (findTranslationComponent(target.dof,component))
        {
            q[nFree+index] = this->nodePrescribedDisplacement[target.nodeId][component];
            const int motion = this->motionIndexByNodeDof[3*target.nodeId+component];
            if (motion >= 0)
            {
                const RStructuralMotionState evaluated =
                    this->scaledMotionState(this->prescribedMotions[uint(motion)],state.time);
                q[nFree+index] = evaluated.displacement;
                v[nFree+index] = evaluated.velocity;
                a[nFree+index] = evaluated.acceleration;
            }
        }
        else if (findRotationComponent(target.dof,component))
        {
            q[nFree+index] = this->nodePrescribedRotation[target.nodeId][component];
        }
    }
    const RStructuralDynamicsSetup &setup = this->pModel->getProblemSetup()
        .getStructuralDynamicsSetup();
    const double forceScale = this->scales.findScaleFactor(R_VARIABLE_FORCE_X);
    const double velocityScale = this->scales.findScaleFactor(R_VARIABLE_VELOCITY_X);
    std::vector<RStructuralSupportSample> samples;
    for (uint row=0;row<this->heldDofs.size();row++)
    {
        if (!this->pointSupportMask[row]) continue;
        double reaction = -heldForce[row];
        for (uint position=0;position<stiffness.getNColumns(row);position++)
        {
            const uint column = stiffness.getVector(row).getIndex(position);
            const double coefficient = stiffness.getValue(row,position);
            reaction += coefficient*(q[column]+setup.getRayleighBeta()*v[column]);
        }
        for (uint position=0;position<mass.getNColumns(row);position++)
        {
            const uint column = mass.getVector(row).getIndex(position);
            const double coefficient = mass.getValue(row,position);
            reaction += coefficient*(a[column]+setup.getRayleighAlpha()*v[column]);
        }
        const RStructuralDofRef &target = this->heldDofs[row];
        const uint component = uint(target.dof);
        const int motion = this->motionIndexByNodeDof[3*target.nodeId+component];
        RStructuralSupportSample sample;
        if (motion >= 0) sample.motionId = this->prescribedMotions[uint(motion)].id;
        sample.target = target;
        sample.time = state.time/this->scales.getSecond();
        sample.reaction = reaction/forceScale;
        sample.velocity = v[nFree+row]/velocityScale;
        sample.power = sample.reaction*sample.velocity;
        if (!std::isfinite(sample.reaction) || !std::isfinite(sample.power))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Support response contains a non-finite reaction or power.");
        }
        samples.push_back(std::move(sample));
    }
    return samples;
}

QString RSolverStress::dynamicSystemFingerprint() const
{
    QCryptographicHash hash(QCryptographicHash::Sha256);
    auto appendInt = [&hash](quint64 value)
    {
        hash.addData(QByteArray::number(value));
        hash.addData(";",1);
    };
    auto appendReal = [&hash](double value)
    {
        // Ordinary model ASCII files use six significant digits for material
        // values. A five-digit operator signature survives their round trip.
        hash.addData(QByteArray::number(value == 0.0 ? 0.0 : value,'g',5));
        hash.addData(";",1);
    };
    auto appendMotionReal = [&hash](double value)
    {
        // Motion knots use 17-digit ASCII serialization. Keep their restart
        // signature exact while preserving the legacy operator tolerance.
        hash.addData(QByteArray::number(value == 0.0 ? 0.0 : value,'g',17));
        hash.addData(";",1);
    };
    auto appendSparse = [&](const RSparseMatrix &matrix)
    {
        appendInt(matrix.getNRows());
        for (uint row=0;row<matrix.getNRows();row++)
        {
            std::vector<std::pair<uint,double>> entries;
            for (uint position=0;position<matrix.getNColumns(row);position++)
            {
                entries.emplace_back(matrix.getVector(row).getIndex(position),
                                     matrix.getValue(row,position));
            }
            std::sort(entries.begin(),entries.end());
            appendInt(entries.size());
            for (const auto &entry : entries)
            {
                appendInt(entry.first);
                appendReal(entry.second);
            }
        }
    };
    hash.addData("range-fea/dynamics/system/v1;");
    const RStructuralDynamicsSetup &setup = this->pModel->getProblemSetup().getStructuralDynamicsSetup();
    appendInt(static_cast<uint>(setup.getMethod()));
    appendReal(setup.getRhoInfinity());
    appendReal(setup.getRayleighAlpha());
    appendReal(setup.getRayleighBeta());
    if (!setup.getModalDampingRatios().empty())
    {
        // Modal ratios are written with 17 significant digits in new ASCII
        // files. Preserve their exact restart identity without changing the
        // legacy Rayleigh-only fingerprint.
        hash.addData("modal-damping/v1;");
        appendInt(setup.getModalDampingRatios().size());
        for (double ratio : setup.getModalDampingRatios()) appendMotionReal(ratio);
    }
    appendSparse(this->A);
    appendSparse(this->M);
    appendInt(this->constraintNullspace.getNRows());
    appendInt(this->constraintNullspace.getNColumns());
    for (uint row=0;row<this->constraintNullspace.getNRows();row++)
    {
        for (uint column=0;column<this->constraintNullspace.getNColumns();column++)
        {
            appendReal(this->constraintNullspace[row][column]);
        }
    }
    appendInt(this->constraintParticular.size());
    for (double value : this->constraintParticular) appendReal(value);
    appendInt(this->nodeConstrainedDirections.size());
    for (uint node=0;node<this->nodeConstrainedDirections.size();node++)
    {
        appendInt(this->nodeConstrainedDirections[node]);
        appendInt(this->nodeConstrainedRotations[node]);
        appendInt(this->nodeHasRotationCondition[node] ? 1 : 0);
        for (uint component=0;component<3;component++)
        {
            appendReal(this->nodePrescribedDisplacement[node][component]);
            appendReal(this->nodePrescribedRotation[node][component]);
        }
    }
    const auto &motions = this->pModel->getProblemSetup().getStructuralInteractionSetup()
        .getPrescribedMotions();
    if (!motions.empty())
    {
        hash.addData("prescribed-motion/v1;");
        appendInt(motions.size());
        for (const RStructuralPrescribedMotion &motion : motions)
        {
            const QByteArray id = motion.id.toUtf8();
            appendInt(id.size());
            hash.addData(id);
            appendInt(motion.target.nodeId);
            appendInt(static_cast<uint>(motion.target.dof));
            appendInt(motion.knots.size());
            for (const RStructuralMotionKnot &knot : motion.knots)
            {
                appendMotionReal(knot.time);
                appendMotionReal(knot.displacement);
                appendMotionReal(knot.velocity);
                appendMotionReal(knot.acceleration);
            }
        }
        appendSparse(this->motionStiffness);
        appendSparse(this->motionMass);
    }
    return QString::fromLatin1(hash.result().toHex());
}

QString RSolverStress::supportResponseFingerprint() const
{
    QCryptographicHash hash(QCryptographicHash::Sha256);
    auto appendInt = [&hash](quint64 value)
    {
        hash.addData(QByteArray::number(value));
        hash.addData(";",1);
    };
    auto appendMatrix = [&](const RSparseMatrix &matrix)
    {
        appendInt(matrix.getNRows());
        for (uint row=0;row<matrix.getNRows();row++)
        {
            std::vector<std::pair<uint,double>> entries;
            for (uint position=0;position<matrix.getNColumns(row);position++)
                entries.emplace_back(matrix.getVector(row).getIndex(position),
                                     matrix.getValue(row,position));
            std::sort(entries.begin(),entries.end());
            appendInt(entries.size());
            for (const auto &entry : entries)
            {
                appendInt(entry.first);
                hash.addData(QByteArray::number(entry.second == 0.0 ? 0.0 : entry.second,'g',5));
                hash.addData(";",1);
            }
        }
    };
    hash.addData("range-fea/dynamics/support-response/v1/5sig;");
    appendInt(this->heldDofs.size());
    for (uint index=0;index<this->heldDofs.size();index++)
    {
        appendInt(this->heldDofs[index].nodeId);
        appendInt(uint(this->heldDofs[index].dof));
        appendInt(this->pointSupportMask[index] ? 1 : 0);
    }
    appendMatrix(this->supportStiffness);
    appendMatrix(this->supportMass);
    return QString::fromLatin1(hash.result().toHex());
}

void RSolverStress::solveStructuralDynamics()
{
    const RStructuralDynamicsSetup &setup = this->pModel->getProblemSetup().getStructuralDynamicsSetup();
    const RStructuralControlSetup &controlSetup = this->pModel->getProblemSetup().getStructuralControlSetup();
    if (!this->cachedControlCatalogFingerprint.isEmpty() &&
        this->controlCatalogFingerprint() != this->cachedControlCatalogFingerprint)
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural control catalog changed after mapping was frozen.");
    if (controlSetup.enabled && !this->controlCallback)
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Enabled structural control has no registered provider.");
    const bool useModalDamping = !setup.getModalDampingRatios().empty();
    if (useModalDamping &&
        (setup.getRayleighAlpha() != 0.0 || setup.getRayleighBeta() != 0.0))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Modal and Rayleigh damping cannot be combined in the first modal-damping implementation.");
    }
    if (useModalDamping && !this->prescribedMotions.empty())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Modal damping with moving support requires modal free/support damping blocks and is not yet supported.");
    }
    if (setup.getMethod() == RStructuralTimeIntegrationMethod::Newmark &&
        (std::fabs(setup.getNewmarkGamma()-0.5) > 1.0e-12 ||
         std::fabs(setup.getNewmarkBeta()-0.25) > 1.0e-12))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "The structural dynamics baseline supports Newmark average acceleration (gamma=0.5, beta=0.25) only.");
    }
    RStructuralDynamicsIntegrator::Options options;
    switch (setup.getMethod())
    {
        case RStructuralTimeIntegrationMethod::Newmark:
            options.method = RStructuralDynamicsIntegrator::Method::NewmarkAverageAcceleration;
            break;
        case RStructuralTimeIntegrationMethod::GeneralizedAlpha:
            options.method = RStructuralDynamicsIntegrator::Method::GeneralizedAlpha;
            break;
        case RStructuralTimeIntegrationMethod::CentralDifference:
            options.method = RStructuralDynamicsIntegrator::Method::CentralDifference;
            break;
    }
    options.rhoInfinity = setup.getRhoInfinity();
    auto sameMatrix = [](const RSparseMatrix &reference, const RSparseMatrix &current)
    {
        if (reference.getNRows() != current.getNRows())
        {
            return false;
        }
        auto rowsMatch = [](const RSparseMatrix &left, const RSparseMatrix &right)
        {
            for (uint row=0;row<left.getNRows();row++)
            {
                for (uint position=0;position<left.getNColumns(row);position++)
                {
                    const uint column = left.getVector(row).getIndex(position);
                    const double value = left.getValue(row,position);
                    const double other = right.findValue(row,column);
                    if (std::fabs(value-other) >
                        1.0e-10*std::max({1.0,std::fabs(value),std::fabs(other)}))
                    {
                        return false;
                    }
                }
            }
            return true;
        };
        return rowsMatch(reference,current) && rowsMatch(current,reference);
    };
    if (!sameMatrix(this->dynamicMassReference,this->M) ||
        !sameMatrix(this->dynamicStiffnessReference,this->A) ||
        !sameMatrix(this->dynamicMotionMassReference,this->motionMass) ||
        !sameMatrix(this->dynamicMotionStiffnessReference,this->motionStiffness) ||
        (this->supportResponseObserver &&
         (!sameMatrix(this->supportStiffnessReference,this->supportStiffness) ||
          !sameMatrix(this->supportMassReference,this->supportMass))))
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural M/K changed during fixed-matrix dynamics; time-varying material or geometry needs matrix refresh.");
    }
    if (this->dynamicConstraintNullspaceReference.getNRows() != this->constraintNullspace.getNRows() ||
        this->dynamicConstraintNullspaceReference.getNColumns() != this->constraintNullspace.getNColumns() ||
        this->dynamicConstraintParticularReference.size() != this->constraintParticular.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural MPC layout changed during dynamics.");
    }
    for (uint row=0;row<this->constraintNullspace.getNRows();row++)
    {
        for (uint column=0;column<this->constraintNullspace.getNColumns();column++)
        {
            if (std::fabs(this->constraintNullspace[row][column]-
                          this->dynamicConstraintNullspaceReference[row][column]) > 1.0e-12)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Structural MPC coefficients changed during dynamics.");
            }
        }
    }
    for (uint row=0;row<this->constraintParticular.size();row++)
    {
        if (std::fabs(this->constraintParticular[row]-this->dynamicConstraintParticularReference[row]) > 1.0e-12)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Moving affine MPC is not yet supported in second-order dynamics.");
        }
    }
    this->currentDynamicForce = this->b;
    this->currentDynamicTime = this->pModel->getTimeSolver().getCurrentTime();
    if (controlSetup.enabled && !this->controlForceMap)
        this->prepareControlMapping();
    if (!this->dynamicsIntegrator)
    {
        const RStructuralDynamicsRestartData &restartData =
            this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
        const QString dofFingerprint = RStructuralDynamicsRestartData::computeDofFingerprint(
            this->orderedPhysicalDofs());
        if (this->pModel->getProblemSetup().getRestart() &&
            (!restartData.available || restartData.physicalDofCount != this->nPhysicalEquations ||
             restartData.dofFingerprint != dofFingerprint ||
             restartData.systemFingerprint != this->dynamicSystemFingerprint() ||
             restartData.q.size() != this->A.getNRows()))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural restart checkpoint does not match the equation order, reduced operator or held supports.");
        }
        if (this->supportResponseObserver && this->pModel->getProblemSetup().getRestart() &&
            !restartData.supportResponseAvailable)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Support-work restart requires a persisted accepted response ledger; this checkpoint has none.");
        }
        if (this->supportResponseObserver && this->pModel->getProblemSetup().getRestart() &&
            restartData.supportFingerprint != this->supportResponseFingerprint())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Support-work restart constrained operator or support identity changed.");
        }
        if (controlSetup.enabled && this->pModel->getProblemSetup().getRestart())
        {
            if (!restartData.controlAvailable ||
                restartData.controlProtocolVersion != controlSetup.controlProtocolVersion ||
                restartData.controlDeterministicSeed != controlSetup.deterministicSeed ||
                restartData.controlFingerprint != this->cachedControlFingerprint ||
                restartData.controlChannelIds.size() != this->activeControlCommands.size())
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Structural control restart protocol, seed, actuator map or ledger does not match.");
            std::vector<RStructuralControlChannelSample> samples;
            for (size_t i=0;i<this->activeControlCommands.size();i++)
            {
                if (restartData.controlChannelIds[i] != this->activeControlCommands[i].actuatorId)
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Structural control restart actuator order changed.");
                samples.push_back({restartData.controlChannelIds[i],
                                   restartData.controlLastCommands[i],
                                   restartData.controlPowers[i],
                                   restartData.cumulativeControlWork[i]});
                this->activeControlCommands[i].value = restartData.controlLastCommands[i];
            }
            this->controlLedger->restore(samples);
            const RStructuralControlStateEnvelope envelope =
                RStructuralControlStateEnvelope::deserialize(restartData.controlControllerState);
            this->acceptedControllerState = envelope.controllerState;
            this->acceptedConditionerState = this->controlConditioner->deserialize(
                envelope.conditionerState);
            if (this->acceptedConditionerState.completedSteps != restartData.completedSteps)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Structural control restart conditioner step count differs from solver checkpoint.");
            this->activeControlForce.fill(0.0);
            this->controlForceMap->addForce(this->activeControlCommands,this->activeControlForce);
        }
        if (this->initialDynamicForce.size() != this->currentDynamicForce.size())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural equation count changed between the initial and first time point.");
        }
        this->previousDynamicTime = this->pModel->getTimeSolver().getPreviousTime();
        this->previousDynamicForce = this->pModel->getProblemSetup().getRestart()
            ? RRVector(restartData.externalForce) : this->initialDynamicForce;
        if (this->pModel->getProblemSetup().getRestart() && !this->prescribedMotions.empty())
        {
            RRVector supportForce;
            this->motionEffectiveForce(this->previousDynamicTime,supportForce);
            for (uint row=0;row<this->previousDynamicForce.size();row++)
            {
                // Checkpoint force is effective; interpolation stores bare f_f.
                this->previousDynamicForce[row] -= supportForce[row];
            }
        }
        if (controlSetup.enabled && this->pModel->getProblemSetup().getRestart())
        {
            for (uint row=0;row<this->previousDynamicForce.size();row++)
                this->previousDynamicForce[row] -= this->activeControlForce[row];
        }
        RSparseMatrix damping;
        if (useModalDamping)
        {
            const RModalDampingResult modal = RModalDamping::synthesize(
                this->M,this->A,setup.getModalDampingRatios(),
                this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));
            damping = modal.damping;
            RLogger::info("Modal damping: %u validated low modes of %u reduced equations; unspecified higher modes have zero modal damping.\n",
                          uint(modal.eigenvalues.size()),this->M.getNRows());
        }
        else
        {
            damping.setNRows(this->M.getNRows());
            for (uint row=0;row<this->M.getNRows();row++)
            {
                for (uint position=0;position<this->M.getNColumns(row);position++)
                {
                    damping.addValue(row,this->M.getVector(row).getIndex(position),
                                     setup.getRayleighAlpha()*this->M.getValue(row,position));
                }
            }
            for (uint row=0;row<this->A.getNRows();row++)
            {
                for (uint position=0;position<this->A.getNColumns(row);position++)
                {
                    damping.addValue(row,this->A.getVector(row).getIndex(position),
                                     setup.getRayleighBeta()*this->A.getValue(row,position));
                }
            }
        }
        auto force = [this, controlEnabled=controlSetup.enabled](
            double time, const RStructuralDynamicsState &, RRVector &values)
        {
            values = this->currentDynamicForce;
            const double span = this->currentDynamicTime-this->previousDynamicTime;
            if (span > 0.0 && this->previousDynamicForce.size() == values.size())
            {
                const double weight = std::max(0.0,std::min(1.0,(time-this->previousDynamicTime)/span));
                for (uint i=0;i<values.size();i++)
                {
                    values[i] = (1.0-weight)*this->previousDynamicForce[i]+weight*values[i];
                }
            }
            if (!this->prescribedMotions.empty())
            {
                RRVector supportForce;
                this->motionEffectiveForce(time,supportForce);
                for (uint row=0;row<values.size();row++)
                {
                    values[row] += supportForce[row];
                }
            }
            if (controlEnabled)
            {
                for (uint row=0;row<values.size();row++)
                    values[row] += this->activeControlForce[row];
            }
        };
        this->dynamicsIntegrator.reset(new RStructuralDynamicsIntegrator(
            this->M,damping,this->A,this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG),options,force));
        const RRVector initialPhysical = this->gatherNodalState(this->nodeDisplacement,this->nodeRotation);
        const RRVector initialVelocityPhysical = this->gatherNodalState(this->nodeVelocity,this->nodeAngularVelocity);
        RLogger::info("Initial structural dynamics state: ||q||=%g ||v||=%g\n",
                      RRVector::euclideanNorm(initialPhysical),RRVector::euclideanNorm(initialVelocityPhysical));
        if (this->pModel->getProblemSetup().getRestart())
        {
            RStructuralDynamicsCheckpoint checkpoint;
            checkpoint.state.time = restartData.time;
            checkpoint.state.q = RRVector(restartData.q);
            checkpoint.state.v = RRVector(restartData.v);
            checkpoint.state.a = RRVector(restartData.a);
            checkpoint.state.externalForce = RRVector(restartData.externalForce);
            checkpoint.state.externalWork = restartData.cumulativeExternalWork;
            checkpoint.state.dampingDissipation = restartData.cumulativeDampingDissipation;
            checkpoint.fixedStep = restartData.fixedTimeStep;
            checkpoint.initialMechanicalEnergy = restartData.initialMechanicalEnergy;
            this->dynamicsIntegrator->restore(checkpoint);
        }
        else
        {
            this->dynamicsIntegrator->initialize(
                this->projectPhysicalState(initialPhysical,true),
                this->projectPhysicalState(initialVelocityPhysical,false),
                this->previousDynamicTime);
        }
        if (this->supportResponseObserver)
        {
            this->acceptedSupportSamples = this->computeSupportSamples(
                this->dynamicsIntegrator->acceptedState(),this->initialSupportForce,true);
            if (this->pModel->getProblemSetup().getRestart())
            {
                if (restartData.supportReactions.size() != this->acceptedSupportSamples.size())
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Support-work restart checkpoint has a different support layout.");
                }
                for (uint index=0;index<this->acceptedSupportSamples.size();index++)
                {
                    RStructuralSupportSample &sample = this->acceptedSupportSamples[index];
                    const double expected = restartData.supportReactions[index];
                    const double tolerance = 1.0e-5*std::max({1.0,std::fabs(expected),
                                                               std::fabs(sample.reaction)});
                    if (std::fabs(sample.reaction-expected) > tolerance ||
                        std::fabs(sample.power-restartData.supportPowers[index]) >
                            1.0e-5*std::max({1.0,std::fabs(sample.power),
                                             std::fabs(restartData.supportPowers[index])}))
                    {
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                     "Support-work restart reaction/power does not match the accepted model state.");
                    }
                    sample.cumulativeWork = restartData.cumulativeSupportWork[index];
                }
            }
        }
    }
    else if (this->currentDynamicForce.size() != this->dynamicsIntegrator->acceptedState().q.size())
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural dynamics equation count changed during time integration.");
    }

    const double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();
    if (!std::isfinite(dt) || dt <= 0.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Structural dynamics requires positive dt (enabled=%d, time=%g, previous=%g, input dt=%g).",
                     int(this->pModel->getTimeSolver().getEnabled()),
                     this->pModel->getTimeSolver().getCurrentTime(),
                     this->pModel->getTimeSolver().getPreviousTime(),
                     this->pModel->getTimeSolver().getInputTimeStepSize());
    }
    RRVector controlStartingVelocity;
    if (controlSetup.enabled)
    {
        const RStructuralDynamicsState &accepted = this->dynamicsIntegrator->acceptedState();
        controlStartingVelocity = accepted.v;
        RStructuralControlObservation observation;
        observation.acceptedTime = accepted.time/this->scales.getSecond();
        observation.nextTime = (accepted.time+dt)/this->scales.getSecond();
        observation.deterministicSeed = controlSetup.deterministicSeed;
        observation.sensors = this->sampleControlSensors(accepted);
        observation.priorControllerState = this->acceptedControllerState;
        const RTimeSolver frozenTimeSolver = this->pModel->getTimeSolver();
        const RStructuralControlSetup frozenControlSetup = controlSetup;
        const RStructuralControlDecision decision = this->controlCallback(observation);
        if (!sameControlTimeSolver(frozenTimeSolver,this->pModel->getTimeSolver()) ||
            !sameControlSetup(frozenControlSetup,
                this->pModel->getProblemSetup().getStructuralControlSetup()))
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural control callback changed the frozen time solver or control setup.");
        if (this->controlCatalogFingerprint() != this->cachedControlCatalogFingerprint)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Structural control callback changed the frozen sensor/actuator catalog.");
        if (decision.requestedCommands.size() != this->activeControlCommands.size() ||
            decision.nextControllerState.size() > 1024*1024)
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Control provider returned wrong channel count or oversized controller state.");
        const auto conditioned = this->controlConditioner->advance(
            this->acceptedConditionerState,decision.requestedCommands);
        this->candidateConditionerState = conditioned.candidate;
        std::vector<RStructuralControlCommand> applied;
        applied.reserve(conditioned.appliedCommands.size());
        const auto &actuators = this->pModel->getProblemSetup().getStructuralInteractionSetup().getActuators();
        bool changed = false;
        for (size_t i=0;i<conditioned.appliedCommands.size();i++)
        {
            const auto &request = conditioned.appliedCommands[i];
            const auto &actuator = actuators[i];
            if (request.actuatorId != this->activeControlCommands[i].actuatorId ||
                !std::isfinite(request.value) || std::isnan(actuator.minimum) ||
                std::isnan(actuator.maximum) || std::isnan(actuator.rateLimit) ||
                actuator.minimum > actuator.maximum || actuator.rateLimit < 0.0)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Control provider command or actuator bounds are invalid for '%s'.",
                             actuator.id.toUtf8().constData());
            const double previous = this->controlLedger->accepted()[i].command;
            const double maxDelta = actuator.rateLimit*(dt/this->scales.getSecond());
            const double lower = std::max(actuator.minimum,previous-maxDelta);
            const double upper = std::min(actuator.maximum,previous+maxDelta);
            if (lower > upper)
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Control actuator '%s' has incompatible magnitude and rate limits.",
                             actuator.id.toUtf8().constData());
            const double value = std::clamp(request.value,lower,upper);
            changed = changed || value != this->activeControlCommands[i].value;
            applied.push_back({request.actuatorId,value});
        }
        if (changed)
        {
            RRVector force(this->A.getNRows(),0.0);
            this->controlForceMap->addForce(applied,force);
            this->activeControlForce = std::move(force);
            this->activeControlCommands = applied;
            // A command jump is a force discontinuity at t_n. Refresh a_n
            // before Newmark predicts, otherwise the first interval looks
            // like a ramp and the held-command work ledger is inconsistent.
            this->dynamicsIntegrator->refreshAcceptedForce();
        }
        this->candidateControllerState = decision.nextControllerState;
        this->trialControlRecord = {};
        this->trialControlRecord.acceptedTime = observation.acceptedTime;
        this->trialControlRecord.nextTime = observation.nextTime;
        this->trialControlRecord.sensors = observation.sensors;
        this->trialControlRecord.requestedCommands = decision.requestedCommands;
        this->trialControlRecord.conditionedCommands = conditioned.appliedCommands;
        this->trialControlRecord.controllerStateSha256 =
            QCryptographicHash::hash(decision.nextControllerState,QCryptographicHash::Sha256).toHex();
    }
    this->dynamicsIntegrator->trialStep(dt);
    this->dynamicsIntegrator->acceptStep();
    const RStructuralDynamicsState &state = this->dynamicsIntegrator->acceptedState();
    if (controlSetup.enabled)
    {
        this->controlLedger->trial(*this->controlForceMap,this->activeControlCommands,
                                   controlStartingVelocity,state.v,dt/this->scales.getSecond());
        this->trialControlRecord.appliedChannels = this->controlLedger->candidate();
    }
    if (this->supportResponseObserver)
    {
        this->trialSupportSamples = this->computeSupportSamples(state,this->supportForce);
        if (this->trialSupportSamples.size() != this->acceptedSupportSamples.size())
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Support response layout changed during integration.");
        }
        for (uint index=0;index<this->trialSupportSamples.size();index++)
        {
            const RStructuralSupportSample &previous = this->acceptedSupportSamples[index];
            RStructuralSupportSample &current = this->trialSupportSamples[index];
            if (current.target.nodeId != previous.target.nodeId ||
                current.target.dof != previous.target.dof ||
                current.motionId != previous.motionId)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Support response identity changed during integration.");
            }
            current.cumulativeWork = previous.cumulativeWork+
                0.5*(previous.power+current.power)*(current.time-previous.time);
        }
    }
    this->x = this->expandReducedSolution(state.q,true);
    this->setDisplacement(this->x);
    this->scatterNodalState(this->expandReducedSolution(state.v,false),
                            this->nodeVelocity,this->nodeAngularVelocity);
    this->scatterNodalState(this->expandReducedSolution(state.a,false),
                            this->nodeAcceleration,this->nodeAngularAcceleration);
    this->setPrescribedMotionResults(state.time);
    this->previousDynamicForce = this->currentDynamicForce;
    this->previousDynamicTime = this->currentDynamicTime;
    // The integrator ledger accounts for free-DOF effective-force work, not
    // full-system support reactions, support kinetic energy or support work.
    RLogger::info("Structural dynamics: reduced kinetic=%g strain=%g energy residual=%g (solver units)\n",
                  state.kineticEnergy,state.strainEnergy,state.energyBalanceResidual);
}

void RSolverStress::solveEigenValue()
{
    RLogger::info("Solving eigen-value problem.\n");

    REigenValueSolverConf conf;

    if (this->pModel->getProblemSetup().getModalSetup().getMethod() == R_MODAL_MULTIPLE_MODES)
    {
        conf.setMethod(REigenValueSolverConf::SubspaceIteration);
    }
    else
    {
        conf.setMethod(REigenValueSolverConf::InversePowerIteration);
    }
    conf.setNEigenValues(this->pModel->getProblemSetup().getModalSetup().getNModesToExtract());
    conf.setNIterations(this->pModel->getProblemSetup().getModalSetup().getNIterations());
    conf.setSolverCvgValue(this->pModel->getProblemSetup().getModalSetup().getConvergenceValue());
    conf.setOutputFrequency(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getOutputFrequency());
    conf.setOutputFileName(this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG).getOutputFileName());

    REigenValueSolver solver(conf,this->pModel->getMatrixSolverConf(RMatrixSolverConf::CG));

    try
    {
        RLogger::indent();
        solver.solve(this->M,this->A,this->d,this->ev);
        RLogger::unindent();
    }
    catch (const RError &error)
    {
        RLogger::unindent();
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve eigen-value problem. %s", error.getMessage().toUtf8().constData());
    }
}

void RSolverStress::setDisplacement(const RRVector &v)
{
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        uint position;
        RR3Vector du(0.0,0.0,0.0);
        RR3Vector dr(0.0,0.0,0.0);

        for (uint c=0;c<3;c++)
        {
            if (this->equationBook.getValue(i,translationDof(c),position))
            {
                du[c] = v[position];
            }
            else
            {
                // Constrained direction. The eigen-value problem is homogeneous
                // so prescribed values do not enter the mode shapes.
                du[c] = (this->problemType == R_PROBLEM_STRESS_MODAL) ? 0.0 : this->nodePrescribedDisplacement[i][c];
            }
        }

        for (uint c=0;c<3;c++)
        {
            const RDofType dof = c == 0 ? RDofType::RX
                                : c == 1 ? RDofType::RY
                                         : RDofType::RZ;
            if (this->equationBook.getValue(i,dof,position))
            {
                dr[c] = v[position];
            }
            else
            {
                dr[c] = (this->problemType == R_PROBLEM_STRESS_MODAL) ? 0.0 : this->nodePrescribedRotation[i][c];
            }
        }

        if (this->localRotations[i].isActive())
        {
            this->localRotations[i].rotateResultsVector(du);
        }
        if (this->rotationLocalRotations[i].isActive())
        {
            this->rotationLocalRotations[i].rotateResultsVector(dr);
        }

        if (this->problemType == R_PROBLEM_STRESS_MODAL)
        {
            du[0] += this->nodeInitialDisplacement.x[i];
            du[1] += this->nodeInitialDisplacement.y[i];
            du[2] += this->nodeInitialDisplacement.z[i];
        }

        this->nodeDisplacement.x[i] = du[0];
        this->nodeDisplacement.y[i] = du[1];
        this->nodeDisplacement.z[i] = du[2];
        this->nodeRotation.x[i] = dr[0];
        this->nodeRotation.y[i] = dr[1];
        this->nodeRotation.z[i] = dr[2];
    }
}

void RSolverStress::process()
{
    if (this->problemType == R_PROBLEM_STRESS)
    {
        // Clear any force field from an earlier execution before processing a
        // new stress result. It is republished below only after this solve has
        // completed result recovery for a supported Heat -> Stress -> Modal
        // task sequence.
        this->pSharedData->addData(QStringLiteral("structural-thermal-axial-force-v1"),RRVector());
        this->pSharedData->addData(QStringLiteral("structural-thermal-axial-force-temperature-v1"),RRVector());
    }

    if (this->isSpecializedStaticTrussMode() || this->isStructuralCableDynamicsMode())
    {
        if (!this->hasConverged())
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Specialized truss results cannot be processed before solver convergence.");
        }
        return;
    }

    if (this->problemType == R_PROBLEM_STRESS_MODAL)
    {
        uint modeNum = this->pModel->getProblemSetup().getModalSetup().getMode();

        // The eigen value of K * phi = lambda * M * phi is lambda = omega^2, so
        // the natural frequency is sqrt(lambda) / (2*pi). The modal setup holds
        // a frequency in Hz, not the raw eigen value.
        double eigenValue = this->d[modeNum];
        double frequency = 0.0;

        if (std::isfinite(eigenValue) && eigenValue > 0.0)
        {
            frequency = std::sqrt(eigenValue) / (2.0 * RConstants::pi);
            RLogger::info("Eigen-value = %g, frequency = %g [Hz]\n",eigenValue,frequency);
        }
        else
        {
            RLogger::warning("Mode %u was not resolved - its eigen value is not usable.\n",modeNum+1);
        }

        this->pModel->getProblemSetup().getModalSetup().setFrequency(frequency);

        RRVector v(this->ev.getNColumns(),0.0);
        for (uint j=0;j<this->ev.getNColumns();j++)
        {
            v[j] = this->ev[modeNum][j];
        }
        this->setDisplacement(this->expandReducedSolution(v,false));
    }

    // Initialize force vector
    this->nodeForce.x.resize(this->pModel->getNNodes());
    this->nodeForce.y.resize(this->pModel->getNNodes());
    this->nodeForce.z.resize(this->pModel->getNNodes());

    this->nodeForce.x.fill(0.0);
    this->nodeForce.y.fill(0.0);
    this->nodeForce.z.fill(0.0);

    this->nodeMoment.x.resize(this->pModel->getNNodes(),0.0);
    this->nodeMoment.y.resize(this->pModel->getNNodes(),0.0);
    this->nodeMoment.z.resize(this->pModel->getNNodes(),0.0);
    this->nodeMoment.x.fill(0.0);
    this->nodeMoment.y.fill(0.0);
    this->nodeMoment.z.fill(0.0);

    // Initialize stress vectors
    for (uint i=0;i<6;i++)
    {
        this->elementStress[i].resize(this->pModel->getNElements(),0.0);
        this->elementStress[i].fill(0.0);
    }
    this->elementNormalStress.resize(this->pModel->getNElements(),0.0);
    this->elementShearStress.resize(this->pModel->getNElements(),0.0);
    this->elementVonMisses.resize(this->pModel->getNElements(),0.0);

    // Process line elements.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        RLine &line = this->pModel->getLine(i);
        double lineCrossArea = line.getCrossArea();

        if (lineCrossArea == 0.0)
        {
            continue;
        }

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(line.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = line.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_LINE(element.getType()));
                if (element.getType() == R_ELEMENT_BEAM1 || element.getType() == R_ELEMENT_BEAM2)
                {
                    const double sectionModulusY = line.getSectionModulusY();
                    const double sectionModulusZ = line.getSectionModulusZ();
                    const double torsionalSectionModulus = line.getTorsionalSectionModulus();
                    if (sectionModulusY <= 0.0 || sectionModulusZ <= 0.0 ||
                        torsionalSectionModulus <= 0.0)
                    {
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                     "Beam stress recovery requires positive Wy, Wz and Wt on line '%s'.",
                                     line.getName().toUtf8().constData());
                    }
                    const uint nSpans = element.getType() == R_ELEMENT_BEAM1 ? 1u : 2u;
                    double representativeAxialStress = 0.0;
                    double maximumNormalStress = 0.0;
                    double maximumShearInvariant = 0.0;
                    for (uint span=0;span<nSpans;span++)
                    {
                        uint first = 0;
                        uint second = 1;
                        if (element.getType() == R_ELEMENT_BEAM2)
                        {
                            first = span == 0 ? 0u : 2u;
                            second = span == 0 ? 2u : 1u;
                        }
                        RElement child(R_ELEMENT_BEAM1);
                        child.setNodeId(0,element.getNodeId(first));
                        child.setNodeId(1,element.getNodeId(second));
                        RRMatrix Me,Ke;
                        RRVector thermalForce;
                        buildBeam1Element(line,child,this->pModel->getNodes(),
                                          this->elementElasticityModulus[elementID],
                                          this->elementPoissonRatio[elementID],
                                          this->elementDensity[elementID],
                                          this->elementThermalExpansion[elementID],
                                          this->elementTemperature[elementID]-this->elementEnvironmentTemperature[elementID],
                                          RR3Vector(0.0,0.0,0.0),false,false,Me,Ke,thermalForce);
                        RRVector xe(12,0.0);
                        const uint parent[2] = {first,second};
                        for (uint node=0;node<2;node++)
                        {
                            const uint nodeId = element.getNodeId(parent[node]);
                            xe[6*node+0] = this->nodeDisplacement.x[nodeId];
                            xe[6*node+1] = this->nodeDisplacement.y[nodeId];
                            xe[6*node+2] = this->nodeDisplacement.z[nodeId];
                            xe[6*node+3] = this->nodeRotation.x[nodeId];
                            xe[6*node+4] = this->nodeRotation.y[nodeId];
                            xe[6*node+5] = this->nodeRotation.z[nodeId];
                        }
                        RRVector internalForce;
                        RRMatrix::mlt(Ke,xe,internalForce);
                        for (uint dof=0;dof<internalForce.size();dof++) internalForce[dof] -= thermalForce[dof];
                        #pragma omp critical
                        {
                            for (uint node=0;node<2;node++)
                            {
                                const uint nodeId = element.getNodeId(parent[node]);
                                this->nodeForce.x[nodeId] += internalForce[6*node+0];
                                this->nodeForce.y[nodeId] += internalForce[6*node+1];
                                this->nodeForce.z[nodeId] += internalForce[6*node+2];
                                this->nodeMoment.x[nodeId] += internalForce[6*node+3];
                                this->nodeMoment.y[nodeId] += internalForce[6*node+4];
                                this->nodeMoment.z[nodeId] += internalForce[6*node+5];
                            }
                        }
                        RRMatrix rotation,transformation,transformationT;
                        const double length = buildBeamTransformation(line,child,this->pModel->getNodes(),rotation,transformation);
                        transformationT.transpose(transformation);
                        RRVector localDisplacement,localEndForce;
                        RRMatrix::mlt(transformationT,xe,localDisplacement);
                        RRMatrix::mlt(transformationT,internalForce,localEndForce);
                        const double axialStress = this->elementElasticityModulus[elementID]
                                                 * ((localDisplacement[6]-localDisplacement[0])/length
                                                    - this->elementThermalExpansion[elementID]
                                                    * (this->elementTemperature[elementID]-this->elementEnvironmentTemperature[elementID]));
                        representativeAxialStress += axialStress/double(nSpans);
                        const double bendingMomentY = std::max(std::fabs(localEndForce[4]),std::fabs(localEndForce[10]));
                        const double bendingMomentZ = std::max(std::fabs(localEndForce[5]),std::fabs(localEndForce[11]));
                        const double torsionalMoment = std::max(std::fabs(localEndForce[3]),std::fabs(localEndForce[9]));
                        maximumNormalStress = std::max(maximumNormalStress,
                                                       std::fabs(axialStress)
                                                       + bendingMomentY/sectionModulusY
                                                       + bendingMomentZ/sectionModulusZ);
                        maximumShearInvariant = std::max(maximumShearInvariant,
                                                         std::sqrt(3.0)*torsionalMoment/torsionalSectionModulus);
                    }
                    this->elementStress[0][elementID] = representativeAxialStress;
                    this->elementNormalStress[elementID] = maximumNormalStress;
                    this->elementShearStress[elementID] = maximumShearInvariant;
                    this->elementVonMisses[elementID] = std::sqrt(maximumNormalStress*maximumNormalStress
                                                                  + maximumShearInvariant*maximumShearInvariant);
                    continue;
                }
                if (element.getType() != R_ELEMENT_TRUSS1 && element.getType() != R_ELEMENT_TRUSS2)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Line element type '%s' is not implemented by result recovery.",
                                 RElement::getName(element.getType()).toUtf8().constData());
                }
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                double QeN = 0.0;

                RRMatrix Be(3*element.size(),1);
                RRMatrix BeT(1,3*element.size());

                double E = this->elementElasticityModulus[elementID];
                // Axial stiffness of the truss - used for the nodal force only.
                double De = E * lineCrossArea;

                RRMatrix Rl;
                RRVector tl;
                element.findTransformationMatrix(this->pModel->getNodes(),Rl,tl);
                Rl.invert();

                RRVector lxe(element.size(),0.0);
                for (uint k=0;k<element.size();k++)
                {
                    RR3Vector xg(this->nodeDisplacement.x[element.getNodeId(k)],
                                 this->nodeDisplacement.y[element.getNodeId(k)],
                                 this->nodeDisplacement.z[element.getNodeId(k)]);
                    RR3Vector xl;
                    RRMatrix::mlt(Rl,xg,xl);
                    lxe[k] = xl[0];
                }

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                    RtT.transpose(Rt);

                    Be.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        Be[3*m+0][0] = Rt[3*m+0][m]*dN[m][0]*J[0][0];
                        Be[3*m+1][0] = Rt[3*m+1][m]*dN[m][0]*J[0][0];
                        Be[3*m+2][0] = Rt[3*m+2][m]*dN[m][0]*J[0][0];
                    }
                    BeT.transpose(Be);

                    // Same stiffness as the one assembled in prepare().
                    RRMatrix BeScaled(Be);
                    BeScaled *= De * detJ * shapeFunc.getW();
                    RRMatrix::mlt(BeScaled,BeT,Ke,true);

                    double integValue = 1.0/double(nInp);

                    // Element level stress. The axial strain follows from the
                    // local axial displacements; the stress is E*(eps - alpha*dT)
                    // and must not carry the cross area, which belongs to the
                    // stiffness only.
                    double axialStrain = 0.0;
                    for (uint m=0;m<element.size();m++)
                    {
                        axialStrain += dN[m][0]*J[0][0] * lxe[m];
                    }
                    QeN += E * (axialStrain - this->elementThermalExpansion[elementID] * dT) * integValue;
                }

                // Nodal force is the internal elastic force K*u. There is no
                // acceleration state in the formulation, so no inertia term.
                RRMatrix::mlt(Ke,xe,fe);

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // A truss carries an axial stress only, along its local x axis.
                this->elementStress[0][elementID] = QeN;
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = 0.0;
                this->elementVonMisses[elementID] = std::fabs(QeN);
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }

    // Process surface elements.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        RSurface &surface = this->pModel->getSurface(i);
        double surfaceThickness = surface.getThickness();
        const bool axisymmetric = surface.getStructuralFormulation() == RElementFormulation::Axisymmetric;

        if (!axisymmetric && surfaceThickness == 0.0)
        {
            continue;
        }
        validateSurfaceFormulation(surface);

        std::atomic<bool> abort{false};
        #pragma omp parallel for default(shared)
        for (int64_t j=0;j<int64_t(surface.size());j++)
        {
            if (abort.load(std::memory_order_relaxed))
            {
                continue;
            }
            try
            {
                uint elementID = surface.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_SURFACE(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*2,element.size()*2,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector thermalForce(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                const uint nStrainComponents = axisymmetric ? 4u : 3u;
                RRVector Qe(nStrainComponents,0.0);
                double accumulatedStressWeight = 0.0;
                double QeN = 0.0;
                double QeS = 0.0;
                double QeVM = 0.0;

                RRMatrix B(element.size(),2,0.0);
                RRMatrix Be(element.size()*2,nStrainComponents,0.0);
                RRMatrix BeT(nStrainComponents,element.size()*2,0.0);
                RRMatrix BeD(element.size()*2,nStrainComponents,0.0);
                RRMatrix Ket(element.size()*2,element.size()*2);
                RRMatrix KeRt(element.size()*2,element.size()*2);
                RRVector localThermalForce(element.size()*2,0.0);

                RRMatrix De;

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                buildSurfaceConstitutiveMatrix(surface.getStructuralFormulation(),E,v,De);

                RRVector lxe(element.size()*2);
                if (axisymmetric)
                {
                    for (uint k=0;k<element.size();k++)
                    {
                        lxe[2*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                        lxe[2*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    }
                }
                else
                {
                    RRMatrix Rl;
                    RRVector tl;
                    element.findTransformationMatrix(this->pModel->getNodes(),Rl,tl);
                    Rl.invert();
                    for (uint k=0;k<element.size();k++)
                    {
                        RR3Vector xg(this->nodeDisplacement.x[element.getNodeId(k)],
                                     this->nodeDisplacement.y[element.getNodeId(k)],
                                     this->nodeDisplacement.z[element.getNodeId(k)]);
                        RR3Vector xl;
                        RRMatrix::mlt(Rl,xg,xl);
                        lxe[2*k+0] = xl[0];
                        lxe[2*k+1] = xl[1];
                    }
                }

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt, RtT;
                    double radius = 0.0;
                    double detJ = 0.0;
                    if (axisymmetric)
                    {
                        detJ = findAxisymmetricKinematics(element,this->pModel->getNodes(),shapeFunc,
                                                          B,Rt,radius);
                    }
                    else
                    {
                        detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);
                        B.fill(0.0);
                        for (uint m=0;m<dN.getNRows();m++)
                        {
                            B[m][0] = dN[m][0]*J[0][0] + dN[m][1]*J[0][1];
                            B[m][1] = dN[m][0]*J[1][0] + dN[m][1]*J[1][1];
                        }
                    }
                    RtT.transpose(Rt);
                    const double integrationScale = (axisymmetric
                                                    ? axisymmetricFullAngle*radius
                                                    : surfaceThickness)
                                                  * detJ*shapeFunc.getW();

                    Be.fill(0.0);
                    for (uint m=0;m<element.size();m++)
                    {
                        if (axisymmetric)
                        {
                            Be[2*m+0][0] = B[m][0];
                            Be[2*m+0][2] = N[m]/radius;
                            Be[2*m+0][3] = B[m][1];
                            Be[2*m+1][1] = B[m][1];
                            Be[2*m+1][3] = B[m][0];
                        }
                        else
                        {
                            Be[2*m+0][0] = B[m][0];
                            Be[2*m+1][1] = B[m][1];
                            Be[2*m+0][2] = B[m][1];
                            Be[2*m+1][2] = B[m][0];
                        }
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    Ket *= integrationScale;
                    RRMatrix::mlt(Rt,Ket,KeRt);
                    RRMatrix::mlt(KeRt,RtT,Ke,true);

                    localThermalForce.fill(0.0);
                    const double thermalStrain = this->elementThermalExpansion[elementID] * dT;
                    const uint nThermalNormals = axisymmetric ? 3u : 2u;
                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint q=0;q<nThermalNormals;q++)
                        {
                            localThermalForce[2*m+0] += thermalStrain*integrationScale*BeD[2*m+0][q];
                            localThermalForce[2*m+1] += thermalStrain*integrationScale*BeD[2*m+1][q];
                        }
                    }
                    RRMatrix::mlt(Rt,localThermalForce,thermalForce,true);

                    const double stressWeight = axisymmetric
                                              ? integrationScale
                                              : 1.0/double(nInp);
                    accumulatedStressWeight += stressWeight;
                    RRVector pointStress(nStrainComponents,0.0);

                    // Element level stress.
                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<nStrainComponents;n++)
                        {
                            pointStress[n] += BeD[2*m+0][n] * lxe[2*m+0]
                                            + BeD[2*m+1][n] * lxe[2*m+1];
                        }
                    }
                    for (uint n=0;n<nStrainComponents;n++)
                    {
                        for (uint q=0;q<nThermalNormals;q++)
                        {
                            pointStress[n] -= thermalStrain*De[n][q];
                        }
                        Qe[n] += pointStress[n]*stressWeight;
                    }
                }
                if (axisymmetric)
                {
                    for (uint n=0;n<nStrainComponents;n++)
                    {
                        Qe[n] /= accumulatedStressWeight;
                    }
                }

                // Nodal force is the internal elastic force K*u. There is no
                // acceleration state in the formulation, so no inertia term.
                RRMatrix::mlt(Ke,xe,fe);
                for (uint m=0;m<fe.size();m++)
                {
                    fe[m] -= thermalForce[m];
                }

                if (axisymmetric)
                {
                    QeN = std::sqrt(Qe[0]*Qe[0] + Qe[1]*Qe[1] + Qe[2]*Qe[2]
                                     - Qe[0]*Qe[1] - Qe[1]*Qe[2] - Qe[2]*Qe[0]);
                    QeS = std::sqrt(3.0*Qe[3]*Qe[3]);
                    QeVM = std::sqrt(QeN*QeN+QeS*QeS);
                }
                else
                {
                    QeN = std::sqrt(Qe[0] * Qe[0] + Qe[1] * Qe[1] - Qe[0] * Qe[1]);
                    QeS = std::sqrt(3.0 * Qe[2] * Qe[2]);
                    // Von Mises combines the normal and the shear invariant in
                    // quadrature, not by adding them.
                    QeVM = std::sqrt(QeN * QeN + QeS * QeS);
                }

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // In-plane components, in the local element frame.
                this->elementStress[0][elementID] = Qe[0];
                this->elementStress[1][elementID] = Qe[1];
                if (axisymmetric)
                {
                    this->elementStress[2][elementID] = Qe[2];
                    this->elementStress[5][elementID] = Qe[3];
                }
                else
                {
                    this->elementStress[5][elementID] = Qe[2];
                }
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = QeS;
                this->elementVonMisses[elementID] = QeVM;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }

    // Process volume elements.
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
                uint elementID = volume.get(uint(j));

                if (!this->computableElements[elementID])
                {
                    continue;
                }

                const RElement &element = this->pModel->getElement(elementID);
                R_ERROR_ASSERT(R_ELEMENT_TYPE_IS_VOLUME(element.getType()));
                uint nInp = RElement::getNIntegrationPoints(element.getType());
                RRMatrix Ke(element.size()*3,element.size()*3,0.0);
                RRVector fe(element.size()*3,0.0);
                RRVector xe(element.size()*3,0.0);
                RRVector Qe(6,0.0);

                RRMatrix B(element.size(),3);
                RRMatrix Be(element.size()*3,6);
                RRMatrix BeT(6,element.size()*3);
                RRMatrix BeD(element.size()*3,6);
                RRMatrix Ket(element.size()*3,element.size()*3);

                RRMatrix De(6,6,0.0);

                double E = this->elementElasticityModulus[elementID];
                double v = this->elementPoissonRatio[elementID];

                De[0][0] = 1-v;   De[0][1] = v;     De[0][2] = v;
                De[1][0] = v;     De[1][1] = 1-v;   De[1][2] = v;
                De[2][0] = v;     De[2][1] = v;     De[2][2] = 1-v;
                De[3][3] = De[4][4] = De[5][5] = (1-2*v)/2;
                De *= E/((1+v)*(1-2*v));

                double dT = this->elementTemperature[elementID] - this->elementEnvironmentTemperature[elementID];

                for (uint k=0;k<element.size();k++)
                {

                    xe[3*k+0] = this->nodeDisplacement.x[element.getNodeId(k)];
                    xe[3*k+1] = this->nodeDisplacement.y[element.getNodeId(k)];
                    xe[3*k+2] = this->nodeDisplacement.z[element.getNodeId(k)];
                }

                for (uint k=0;k<nInp;k++)
                {
                    const RElementShapeFunction &shapeFunc = RElement::getShapeFunction(element.getType(),k);
                    const RRVector &N = shapeFunc.getN();
                    const RRMatrix &dN = shapeFunc.getDN();
                    RRMatrix J, Rt;
                    double detJ = element.findJacobian(this->pModel->getNodes(),k,J,Rt);

                    B.fill(0.0);
                    for (uint m=0;m<dN.getNRows();m++)
                    {
                        B[m][0] += (dN[m][0]*J[0][0] + dN[m][1]*J[0][1] + dN[m][2]*J[0][2]);
                        B[m][1] += (dN[m][0]*J[1][0] + dN[m][1]*J[1][1] + dN[m][2]*J[1][2]);
                        B[m][2] += (dN[m][0]*J[2][0] + dN[m][1]*J[2][1] + dN[m][2]*J[2][2]);
                    }

                    for (uint m=0;m<element.size();m++)
                    {
                        Be[3*m+0][0] = B[m][0];   Be[3*m+1][0] = 0.0;       Be[3*m+2][0] = 0.0;
                        Be[3*m+0][1] = 0.0;       Be[3*m+1][1] = B[m][1];   Be[3*m+2][1] = 0.0;
                        Be[3*m+0][2] = 0.0;       Be[3*m+1][2] = 0.0;       Be[3*m+2][2] = B[m][2];
                        Be[3*m+0][3] = 0.0;       Be[3*m+1][3] = B[m][2];   Be[3*m+2][3] = B[m][1];
                        Be[3*m+0][4] = B[m][2];   Be[3*m+1][4] = 0.0;       Be[3*m+2][4] = B[m][0];
                        Be[3*m+0][5] = B[m][1];   Be[3*m+1][5] = B[m][0];   Be[3*m+2][5] = 0.0;
                    }
                    BeT.transpose(Be);

                    RRMatrix::mlt(Be,De,BeD);
                    RRMatrix::mlt(BeD,BeT,Ket);
                    for (uint m=0;m<3*element.size();m++)
                    {
                        for (uint n=0;n<3*element.size();n++)
                        {
                            // Stiffness matrix
                            Ke[m][n] += Ket[m][n] * detJ * shapeFunc.getW();
                        }
                    }

                    double integValue = 1.0/double(nInp);

                    // Element level stress.
                    for (uint m=0;m<element.size();m++)
                    {
                        for (uint n=0;n<6;n++)
                        {
                            Qe[n] += BeD[3*m+0][n] * this->nodeDisplacement.x[element.getNodeId(m)] * integValue
                                  +  BeD[3*m+1][n] * this->nodeDisplacement.y[element.getNodeId(m)] * integValue
                                  +  BeD[3*m+2][n] * this->nodeDisplacement.z[element.getNodeId(m)] * integValue;
                        }
                        for (uint n=0;n<3;n++)
                        {
                            Qe[0] -= BeD[3*m+0][n] * this->elementThermalExpansion[elementID] * dT * integValue;
                            Qe[1] -= BeD[3*m+1][n] * this->elementThermalExpansion[elementID] * dT * integValue;
                            Qe[2] -= BeD[3*m+2][n] * this->elementThermalExpansion[elementID] * dT * integValue;
                        }
                    }
                }

                // Nodal force is the internal elastic force K*u. There is no
                // acceleration state in the formulation, so no inertia term.
                RRMatrix::mlt(Ke,xe,fe);

                double QeN = std::sqrt(Qe[0]*Qe[0] + Qe[1]*Qe[1] + Qe[2]*Qe[2] - (Qe[0]*Qe[1] + Qe[1]*Qe[2] + Qe[2]*Qe[0]));
                double QeS = std::sqrt(3.0 * (Qe[3]*Qe[3] + Qe[4]*Qe[4] + Qe[5]*Qe[5]));
                // Von Mises combines the normal and the shear invariant in
                // quadrature, not by adding them.
                double QeVM = std::sqrt(QeN*QeN + QeS*QeS);

                #pragma omp critical
                {
                    for (uint m=0;m<element.size();m++)
                    {
                        this->nodeForce.x[element.getNodeId(m)] += fe[3*m+0];
                        this->nodeForce.y[element.getNodeId(m)] += fe[3*m+1];
                        this->nodeForce.z[element.getNodeId(m)] += fe[3*m+2];
                    }
                }

                // Writes below are per-element - no synchronization needed.
                // Stress components in global coordinates.
                for (uint n=0;n<6;n++)
                {
                    this->elementStress[n][elementID] = Qe[n];
                }
                this->elementNormalStress[elementID] = QeN;
                this->elementShearStress[elementID] = QeS;
                this->elementVonMisses[elementID] = QeVM;
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
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to process results.");
        }
    }

    const bool structuralDynamics = this->pModel->getProblemSetup()
                                        .getStructuralDynamicsSetup().getEnabled();
    if (this->problemType == R_PROBLEM_STRESS &&
        !this->pModel->getTimeSolver().getEnabled() && !structuralDynamics &&
        hasCurrentThermalStressBeforeModal(this->pModel->getProblemTaskTree()))
    {
        const uint elementCount = this->pModel->getNElements();
        const double forceScale = this->scales.findScaleFactor(R_VARIABLE_FORCE);
        if (!std::isfinite(forceScale) || forceScale <= 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Thermal axial-force publication scale is invalid.");
        }

        const double unsupported = std::numeric_limits<double>::quiet_NaN();
        RRVector axialForce(elementCount,unsupported);
        for (uint elementId=0;elementId<elementCount;elementId++)
        {
            if (!this->computableElements[elementId] ||
                this->pModel->getElement(elementId).getType() == R_ELEMENT_POINT)
            {
                axialForce[elementId] = 0.0;
            }
        }

        for (uint lineId=0;lineId<this->pModel->getNLines();lineId++)
        {
            const RLine &line = this->pModel->getLine(lineId);
            const double area = line.getCrossArea();
            for (uint position=0;position<line.size();position++)
            {
                const uint elementId = line.get(position);
                const RElementType type = this->pModel->getElement(elementId).getType();
                if (type != R_ELEMENT_TRUSS1 && type != R_ELEMENT_BEAM1)
                {
                    continue;
                }
                if (!this->computableElements[elementId])
                {
                    axialForce[elementId] = 0.0;
                    continue;
                }
                if (!std::isfinite(area) || area <= 0.0)
                {
                    axialForce[elementId] = unsupported;
                    continue;
                }
                const double scaledForce = this->elementStress[0][elementId]*area;
                axialForce[elementId] = std::isfinite(scaledForce)
                                      ? scaledForce/forceScale
                                      : unsupported;
            }
        }

        this->pSharedData->addData(QStringLiteral("structural-thermal-axial-force-v1"),axialForce);
        this->pSharedData->addData(QStringLiteral("structural-thermal-axial-force-temperature-v1"),
                                   this->elementTemperature);
    }
}

void RSolverStress::store()
{
    if (this->isStructuralCableDynamicsMode())
    {
        if (!this->structuralCableDynamicsResultReady)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Cable transient results cannot be published before Newton convergence.");
        }
        RStructuralCableModelAdapter::publishDynamics(
            *this->pModel,this->structuralCableDynamicsResult);
        if (this->controlForceMap)
        {
            RStructuralDynamicsRestartData &checkpoint =
                this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
            const RStructuralControlSetup &controlSetup =
                this->pModel->getProblemSetup().getStructuralControlSetup();
            checkpoint.controlAvailable = true;
            checkpoint.controlProtocolVersion = controlSetup.controlProtocolVersion;
            checkpoint.controlFingerprint = this->cachedControlFingerprint;
            checkpoint.controlChannelIds.clear();
            checkpoint.controlLastCommands.clear();
            checkpoint.controlPowers.clear();
            checkpoint.cumulativeControlWork.clear();
            checkpoint.controlControllerState = RStructuralControlStateEnvelope{
                this->candidateControllerState,
                this->controlConditioner->serialize(this->candidateConditionerState)
            }.serialize();
            checkpoint.controlDeterministicSeed = controlSetup.deterministicSeed;
            for (const RStructuralControlChannelSample &sample :
                 this->controlLedger->candidate())
            {
                checkpoint.controlChannelIds.push_back(sample.actuatorId);
                checkpoint.controlLastCommands.push_back(sample.command);
                checkpoint.controlPowers.push_back(sample.endpointPower);
                checkpoint.cumulativeControlWork.push_back(sample.cumulativeWork);
            }
            checkpoint.validate();
        }
        return;
    }
    if (this->isStructuralPiezoMode())
    {
        if (!this->structuralPiezoResultReady)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Piezo results cannot be published before equilibrium convergence.");
        }
        RStructuralPiezoModelAdapter::publish(*this->pModel,this->structuralPiezoResult);
        return;
    }
    if (this->isStructuralCableMode())
    {
        if (!this->structuralCableResultReady)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Cable results cannot be published before equilibrium convergence.");
        }
        RStructuralCableModelAdapter::publish(*this->pModel,this->structuralCableResult);
        return;
    }
    RLogger::info("Storing results\n");
    RLogger::indent();

    bool hasFrameDofs = false;
    for (uint i=0;i<this->pModel->getNElements() && !hasFrameDofs;i++)
    {
        if (this->computableElements[i] &&
            RElementDescriptor::get(this->pModel->getElement(i).getType()).getDofsPerNode().contains(RDofType::RX))
        {
            hasFrameDofs = true;
        }
    }

    // Displacement
    uint displacementPos = this->pModel->findVariable(R_VARIABLE_DISPLACEMENT);
    if (displacementPos == RConstants::eod)
    {
        displacementPos = this->pModel->addVariable(R_VARIABLE_DISPLACEMENT);

        double umin = 0.0;
        double umax = 0.0;
        for (uint i=0;i<this->nodeDisplacement.x.size();i++)
        {
            double u = RR3Vector(this->nodeDisplacement.x[i],
                                 this->nodeDisplacement.y[i],
                                 this->nodeDisplacement.z[i]).length();
            if (i == 0)
            {
                umin = umax = u;
            }
            else
            {
                umin = std::min(umin,u);
                umax = std::max(umax,u);
            }
        }

        this->pModel->getVariable(displacementPos).getVariableData().setMinMaxDisplayValue(umin,umax);
    }
    RVariable &displacement =  this->pModel->getVariable(displacementPos);

    displacement.setApplyType(R_VARIABLE_APPLY_NODE);
    displacement.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        displacement.setValue(0,i,this->nodeDisplacement.x[i]);
        displacement.setValue(1,i,this->nodeDisplacement.y[i]);
        displacement.setValue(2,i,this->nodeDisplacement.z[i]);
    }

    if (hasFrameDofs)
    {
        // Rotation
        uint rotationPos = this->pModel->findVariable(R_VARIABLE_ROTATION);
        if (rotationPos == RConstants::eod)
        {
            rotationPos = this->pModel->addVariable(R_VARIABLE_ROTATION);
        }
        RVariable &rotation = this->pModel->getVariable(rotationPos);
        rotation.setApplyType(R_VARIABLE_APPLY_NODE);
        rotation.resize(3,this->pModel->getNNodes());
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            rotation.setValue(0,i,this->nodeRotation.x[i]);
            rotation.setValue(1,i,this->nodeRotation.y[i]);
            rotation.setValue(2,i,this->nodeRotation.z[i]);
        }
    }

    if (this->problemType == R_PROBLEM_STRESS &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled())
    {
        auto storeVector = [this](RVariableType type, const RSolverCartesianVector<RRVector> &value)
        {
            uint position = this->pModel->findVariable(type);
            if (position == RConstants::eod)
            {
                position = this->pModel->addVariable(type);
            }
            RVariable &variable = this->pModel->getVariable(position);
            variable.setApplyType(R_VARIABLE_APPLY_NODE);
            variable.resize(3,this->pModel->getNNodes());
            for (uint node=0;node<this->pModel->getNNodes();node++)
            {
                variable.setValue(0,node,value.x[node]);
                variable.setValue(1,node,value.y[node]);
                variable.setValue(2,node,value.z[node]);
            }
        };
        storeVector(R_VARIABLE_VELOCITY,this->nodeVelocity);
        storeVector(R_VARIABLE_ACCELERATION,this->nodeAcceleration);
        if (hasFrameDofs)
        {
            storeVector(R_VARIABLE_ANGULAR_VELOCITY,this->nodeAngularVelocity);
            storeVector(R_VARIABLE_ANGULAR_ACCELERATION,this->nodeAngularAcceleration);
        }
        if (this->dynamicsIntegrator)
        {
            const RStructuralDynamicsCheckpoint committed = this->dynamicsIntegrator->checkpoint();
            const RStructuralDynamicsState &state = committed.state;
            RStructuralDynamicsRestartData &restart =
                this->pModel->getProblemSetup().getStructuralDynamicsRestartData();
            restart.completedSteps = this->pModel->getTimeSolver().getCurrentTimeStep()+1;
            restart.time = state.time;
            restart.fixedTimeStep = committed.fixedStep;
            restart.physicalDofCount = this->nPhysicalEquations;
            restart.dofFingerprint = RStructuralDynamicsRestartData::computeDofFingerprint(
                this->orderedPhysicalDofs());
            restart.systemFingerprint = this->dynamicSystemFingerprint();
            restart.q.assign(state.q.begin(),state.q.end());
            restart.v.assign(state.v.begin(),state.v.end());
            restart.a.assign(state.a.begin(),state.a.end());
            restart.externalForce.assign(state.externalForce.begin(),state.externalForce.end());
            restart.cumulativeExternalWork = state.externalWork;
            restart.cumulativeDampingDissipation = state.dampingDissipation;
            restart.initialMechanicalEnergy = committed.initialMechanicalEnergy;
            restart.supportResponseAvailable = bool(this->supportResponseObserver);
            restart.supportFingerprint = this->supportResponseObserver
                ? this->supportResponseFingerprint() : QString();
            restart.supportReactions.clear();
            restart.supportPowers.clear();
            restart.cumulativeSupportWork.clear();
            for (const RStructuralSupportSample &sample : this->trialSupportSamples)
            {
                restart.supportReactions.push_back(sample.reaction);
                restart.supportPowers.push_back(sample.power);
                restart.cumulativeSupportWork.push_back(sample.cumulativeWork);
            }
            restart.controlAvailable = bool(this->controlForceMap);
            restart.controlProtocolVersion = 1;
            restart.controlFingerprint = this->controlForceMap
                ? this->cachedControlFingerprint : QString();
            restart.controlChannelIds.clear();
            restart.controlLastCommands.clear();
            restart.controlPowers.clear();
            restart.cumulativeControlWork.clear();
            restart.controlControllerState = this->controlForceMap
                ? RStructuralControlStateEnvelope{
                      this->candidateControllerState,
                      this->controlConditioner->serialize(this->candidateConditionerState)}.serialize()
                : QByteArray();
            restart.controlDeterministicSeed = this->controlForceMap
                ? this->pModel->getProblemSetup().getStructuralControlSetup().deterministicSeed : 0;
            if (this->controlForceMap)
            {
                for (const auto &sample : this->controlLedger->candidate())
                {
                    restart.controlChannelIds.push_back(sample.actuatorId);
                    restart.controlLastCommands.push_back(sample.command);
                    restart.controlPowers.push_back(sample.endpointPower);
                    restart.cumulativeControlWork.push_back(sample.cumulativeWork);
                }
            }
            restart.available = true;
            if (restart.q.size() != restart.v.size() ||
                restart.q.size() != restart.a.size() ||
                restart.q.size() != restart.externalForce.size() ||
                restart.q.size() > restart.physicalDofCount)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Structural restart internal size mismatch: physical=%u q=%zu v=%zu a=%zu f=%zu.",
                             restart.physicalDofCount,restart.q.size(),restart.v.size(),
                             restart.a.size(),restart.externalForce.size());
            }
            restart.validate();
        }
    }

    // VonMises Stress
    uint vonMisesStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_VON_MISES);
    if (vonMisesStressPos == RConstants::eod)
    {
        vonMisesStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_VON_MISES);

        this->pModel->getVariable(vonMisesStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementVonMisses),
                    RStatistics::findMaximumValue(this->elementVonMisses));
    }
    RVariable &vonMisesStress =  this->pModel->getVariable(vonMisesStressPos);

    vonMisesStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    vonMisesStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        vonMisesStress.setValue(0,i,this->elementVonMisses[i]);
    }

    // Stress components. Volume elements report them in global coordinates,
    // surface and line elements in their own local element frame.
    const RVariableType stressComponentTypes[6] =
    {
        R_VARIABLE_STRESS_X,
        R_VARIABLE_STRESS_Y,
        R_VARIABLE_STRESS_Z,
        R_VARIABLE_STRESS_YZ,
        R_VARIABLE_STRESS_XZ,
        R_VARIABLE_STRESS_XY
    };

    for (uint c=0;c<6;c++)
    {
        uint stressComponentPos = this->pModel->findVariable(stressComponentTypes[c]);
        if (stressComponentPos == RConstants::eod)
        {
            stressComponentPos = this->pModel->addVariable(stressComponentTypes[c]);

            this->pModel->getVariable(stressComponentPos).getVariableData().setMinMaxDisplayValue(
                        RStatistics::findMinimumValue(this->elementStress[c]),
                        RStatistics::findMaximumValue(this->elementStress[c]));
        }
        RVariable &stressComponent = this->pModel->getVariable(stressComponentPos);

        stressComponent.setApplyType(R_VARIABLE_APPLY_ELEMENT);
        stressComponent.resize(1,this->pModel->getNElements());
        for (uint i=0;i<this->pModel->getNElements();i++)
        {
            stressComponent.setValue(0,i,this->elementStress[c][i]);
        }
    }

    // Normal Stress
    uint normalStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_NORMAL);
    if (normalStressPos == RConstants::eod)
    {
        normalStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_NORMAL);

        this->pModel->getVariable(normalStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementNormalStress),
                    RStatistics::findMaximumValue(this->elementNormalStress));
    }
    RVariable &normalStress =  this->pModel->getVariable(normalStressPos);

    normalStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    normalStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        normalStress.setValue(0,i,this->elementNormalStress[i]);
    }

    // Shear Stress
    uint shearStressPos = this->pModel->findVariable(R_VARIABLE_STRESS_SHEAR);
    if (shearStressPos == RConstants::eod)
    {
        shearStressPos = this->pModel->addVariable(R_VARIABLE_STRESS_SHEAR);

        this->pModel->getVariable(shearStressPos).getVariableData().setMinMaxDisplayValue(
                    RStatistics::findMinimumValue(this->elementShearStress),
                    RStatistics::findMaximumValue(this->elementShearStress));
    }
    RVariable &shearStress =  this->pModel->getVariable(shearStressPos);

    shearStress.setApplyType(R_VARIABLE_APPLY_ELEMENT);
    shearStress.resize(1,this->pModel->getNElements());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        shearStress.setValue(0,i,this->elementShearStress[i]);
    }

    // Force
    uint forcePos = this->pModel->findVariable(R_VARIABLE_FORCE);
    if (forcePos == RConstants::eod)
    {
        forcePos = this->pModel->addVariable(R_VARIABLE_FORCE);

        double fmin = 0.0;
        double fmax = 0.0;
        for (uint i=0;i<this->nodeForce.x.size();i++)
        {
            double f = RR3Vector(this->nodeForce.x[i],
                                 this->nodeForce.y[i],
                                 this->nodeForce.z[i]).length();
            if (i == 0)
            {
                fmin = fmax = f;
            }
            else
            {
                fmin = std::min(fmin,f);
                fmax = std::max(fmax,f);
            }
        }

        this->pModel->getVariable(forcePos).getVariableData().setMinMaxDisplayValue(fmin,fmax);
    }
    RVariable &force =  this->pModel->getVariable(forcePos);

    force.setApplyType(R_VARIABLE_APPLY_NODE);
    force.resize(3,this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        force.setValue(0,i,this->nodeForce.x[i]);
        force.setValue(1,i,this->nodeForce.y[i]);
        force.setValue(2,i,this->nodeForce.z[i]);
    }

    if (hasFrameDofs)
    {
        // Moment
        uint momentPos = this->pModel->findVariable(R_VARIABLE_MOMENT);
        if (momentPos == RConstants::eod)
        {
            momentPos = this->pModel->addVariable(R_VARIABLE_MOMENT);
        }
        RVariable &moment = this->pModel->getVariable(momentPos);
        moment.setApplyType(R_VARIABLE_APPLY_NODE);
        moment.resize(3,this->pModel->getNNodes());
        for (uint i=0;i<this->pModel->getNNodes();i++)
        {
            moment.setValue(0,i,this->nodeMoment.x[i]);
            moment.setValue(1,i,this->nodeMoment.y[i]);
            moment.setValue(2,i,this->nodeMoment.z[i]);
        }
    }

    RLogger::unindent();
}

void RSolverStress::statistics()
{
    if (this->isStructuralPiezoMode())
    {
        if (!this->structuralPiezoResultReady)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Piezo statistics require a converged equilibrium.");
        }
        RLogger::info("Piezo static equilibrium converged: free-residual=%g electric-enthalpy=%g\n",
                      this->structuralPiezoResult.residualNorm,
                      this->structuralPiezoResult.enthalpy);
        this->printStats(R_VARIABLE_DISPLACEMENT);
        this->printStats(R_VARIABLE_STRESS_VON_MISES);
        this->printStats(R_VARIABLE_FORCE);
        this->printStats(R_VARIABLE_ELECTRIC_CHARGE);
        return;
    }
    if (this->isStructuralCableMode())
    {
        if (this->isStructuralCableDynamicsMode())
        {
            if (!this->structuralCableDynamicsResultReady)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Cable transient statistics require a converged Newmark step.");
            }
            RLogger::info("Cable transient step converged: iterations=%zu free-residual=%g "
                          "strain-energy=%g kinetic-energy=%g energy-balance-error=%g\n",
                          this->structuralCableDynamicsResult.iterations,
                          this->structuralCableDynamicsResult.residualNorm,
                          this->structuralCableDynamicsResult.strainEnergy,
                          this->structuralCableDynamicsResult.kineticEnergy,
                          this->structuralCableDynamicsResult.energyBalanceError);
            this->printStats(R_VARIABLE_DISPLACEMENT);
            this->printStats(R_VARIABLE_VELOCITY);
            this->printStats(R_VARIABLE_ACCELERATION);
            this->printStats(R_VARIABLE_STRESS_VON_MISES);
            this->printStats(R_VARIABLE_FORCE);
            return;
        }
        if (!this->structuralCableResultReady)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Cable statistics require a converged equilibrium.");
        }
        RLogger::info("Cable equilibrium converged: iterations=%zu free-residual=%g energy=%g\n",
                      this->structuralCableResult.iterations,
                      this->structuralCableResult.residualNorm,
                      this->structuralCableResult.energy);
        this->printStats(R_VARIABLE_DISPLACEMENT);
        this->printStats(R_VARIABLE_STRESS_VON_MISES);
        this->printStats(R_VARIABLE_FORCE);
        return;
    }
    this->printStats(R_VARIABLE_DISPLACEMENT);
    if (this->pModel->findVariable(R_VARIABLE_ROTATION) != RConstants::eod)
    {
        this->printStats(R_VARIABLE_ROTATION);
    }
    this->printStats(R_VARIABLE_STRESS_VON_MISES);
    this->printStats(R_VARIABLE_FORCE);
    if (this->pModel->findVariable(R_VARIABLE_MOMENT) != RConstants::eod)
    {
        this->printStats(R_VARIABLE_MOMENT);
    }
    this->processMonitoringPoints();
}

bool RSolverStress::isComponentEnabled(const RBoundaryCondition &bc, RVariableType variableType)
{
    uint componentPosition = bc.findComponentPosition(variableType);
    if (componentPosition == RConstants::eod)
    {
        return false;
    }
    return bc.getComponent(componentPosition).getEnabled();
}

double RSolverStress::findComponentValue(const RBoundaryCondition &bc, RVariableType variableType) const
{
    uint componentPosition = bc.findComponentPosition(variableType);
    if (componentPosition == RConstants::eod)
    {
        return 0.0;
    }
    return bc.getComponent(componentPosition).get(this->pModel->getTimeSolver().getCurrentTime());
}

void RSolverStress::updateLocalRotations()
{
    // The local frame of a node follows from the constraints acting on it, so
    // it is built in generateLocalConstraints() instead of from the geometry of
    // a single boundary condition.
}

namespace
{

//! One displacement constraint acting on a node: direction . u = value, with
//! the direction given in global coordinates and of unit length.
struct RNodeConstraint
{
    RR3Vector direction;
    double value;
};

//! Return the boundary condition of an entity which carries a local direction,
//! or null when there is none.
const RBoundaryCondition *findLocalDirectionBc(const RElementGroup &rElementGroup, RProblemType problemType)
{
    for (uint i=0;i<rElementGroup.getNBoundaryConditions();i++)
    {
        const RBoundaryCondition &bc = rElementGroup.getBoundaryCondition(i);
        if ((RBoundaryCondition::getProblemTypeMask(bc.getType()) & problemType) && bc.getHasLocalDirection())
        {
            return &bc;
        }
    }
    return nullptr;
}

}

void RSolverStress::generateLocalConstraints()
{
    uint nNodes = this->pModel->getNNodes();

    this->localRotations.resize(nNodes);
    this->nodeConstrainedDirections.resize(nNodes,0);
    this->nodeConstrainedDirections.fill(0);
    this->nodePrescribedDisplacement.resize(nNodes,3,0.0);
    this->nodePrescribedDisplacement.fill(0.0);
    this->rotationLocalRotations.resize(nNodes);
    this->nodeConstrainedRotations.resize(nNodes,0);
    this->nodeConstrainedRotations.fill(0);
    this->nodeHasRotationCondition.resize(nNodes,false);
    this->nodeHasRotationCondition.fill(false);
    this->nodePrescribedRotation.resize(nNodes,3,0.0);
    this->nodePrescribedRotation.fill(0.0);

    // Directions taken from the geometry. A node lying on several entities
    // collects the average of what they give it, as before.
    std::vector<RR3Vector> nodeNormal(nNodes,RR3Vector(0.0,0.0,0.0));
    RBVector nodeNormalSet(nNodes,false);
    std::vector<RR3Vector> nodeLineDirection(nNodes,RR3Vector(0.0,0.0,0.0));
    RBVector nodeLineDirectionSet(nNodes,false);

    // Surfaces - the averaged element normals, or the entered direction.
    for (uint i=0;i<this->pModel->getNSurfaces();i++)
    {
        const RSurface &rSurface = this->pModel->getSurface(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rSurface,this->problemType);
        if (!pBc)
        {
            continue;
        }

        bool useEntered = pBc->getExplicitLocalDirection();
        RR3Vector entered = pBc->getLocalDirection();

        if (useEntered && entered.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of surface entity \'%s\' has zero length.",
                         rSurface.getName().toUtf8().constData());
        }

        for (uint j=0;j<rSurface.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rSurface.get(j));
            RR3Vector direction;

            if (useEntered)
            {
                direction = entered;
            }
            else if (!rElement.findNormal(this->pModel->getNodes(),direction[0],direction[1],direction[2]))
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Could not calculate element normal for element# = %u.",rSurface.get(j));
            }

            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                nodeNormal[nodeId][0] += direction[0];
                nodeNormal[nodeId][1] += direction[1];
                nodeNormal[nodeId][2] += direction[2];
                nodeNormalSet[nodeId] = true;
            }
        }
    }

    // Points - there is no geometry to follow, the entered direction is used.
    for (uint i=0;i<this->pModel->getNPoints();i++)
    {
        const RPoint &rPoint = this->pModel->getPoint(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rPoint,this->problemType);
        if (!pBc)
        {
            continue;
        }

        RR3Vector direction = pBc->getLocalDirection();
        if (direction.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of point entity \'%s\' has zero length.",
                         rPoint.getName().toUtf8().constData());
        }

        for (uint j=0;j<rPoint.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rPoint.get(j));
            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                nodeNormal[nodeId] = direction;
                nodeNormalSet[nodeId] = true;
            }
        }
    }

    // Lines - the element direction, or the entered direction.
    for (uint i=0;i<this->pModel->getNLines();i++)
    {
        const RLine &rLine = this->pModel->getLine(i);
        const RBoundaryCondition *pBc = findLocalDirectionBc(rLine,this->problemType);
        if (!pBc)
        {
            continue;
        }

        bool useEntered = pBc->getExplicitLocalDirection();
        RR3Vector entered = pBc->getLocalDirection();

        if (useEntered && entered.length() < RConstants::eps)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Local direction of line entity \'%s\' has zero length.",
                         rLine.getName().toUtf8().constData());
        }

        for (uint j=0;j<rLine.size();j++)
        {
            const RElement &rElement = this->pModel->getElement(rLine.get(j));
            if (rElement.size() < 2)
            {
                continue;
            }

            RR3Vector direction;
            if (useEntered)
            {
                direction = entered;
            }
            else
            {
                RR3Vector::subtract(this->pModel->getNode(rElement.getNodeId(1)).toVector(),
                                    this->pModel->getNode(rElement.getNodeId(0)).toVector(),
                                    direction);
            }
            if (direction.normalize() < RConstants::eps)
            {
                continue;
            }

            for (uint k=0;k<rElement.size();k++)
            {
                uint nodeId = rElement.getNodeId(k);
                // Keep the accumulated direction consistently oriented, so that
                // a polyline does not average itself away.
                double sign = 1.0;
                if (nodeLineDirectionSet[nodeId] && RR3Vector::dot(nodeLineDirection[nodeId],direction) < 0.0)
                {
                    sign = -1.0;
                }
                nodeLineDirection[nodeId][0] += sign*direction[0];
                nodeLineDirection[nodeId][1] += sign*direction[1];
                nodeLineDirection[nodeId][2] += sign*direction[2];
                nodeLineDirectionSet[nodeId] = true;
            }
        }
    }

    for (uint i=0;i<nNodes;i++)
    {
        if (nodeNormalSet[i])
        {
            nodeNormalSet[i] = (nodeNormal[i].normalize() > RConstants::eps);
        }
        if (nodeLineDirectionSet[i])
        {
            nodeLineDirectionSet[i] = (nodeLineDirection[i].normalize() > RConstants::eps);
        }
    }

    // Collect the constraints of every node, in global coordinates.
    std::vector<std::vector<RNodeConstraint> > nodeConstraints(nNodes);

    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        REntityGroupType entityType = this->pModel->getEntityGroupType(i);
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }

        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (!(RBoundaryCondition::getProblemTypeMask(bc.getType()) & R_PROBLEM_STRESS))
            {
                continue;
            }
            if (bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT &&
                bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT_NORMAL &&
                bc.getType() != R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER)
            {
                continue;
            }

            for (uint k=0;k<pElementGroup->size();k++)
            {
                const RElement &rElement = this->pModel->getElement(pElementGroup->get(k));

                for (uint l=0;l<rElement.size();l++)
                {
                    uint nodeId = rElement.getNodeId(l);
                    std::vector<RNodeConstraint> &constraints = nodeConstraints[nodeId];

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT)
                    {
                        // Each switched on component holds one global direction.
                        const RVariableType components[3] =
                        {
                            R_VARIABLE_DISPLACEMENT_X,
                            R_VARIABLE_DISPLACEMENT_Y,
                            R_VARIABLE_DISPLACEMENT_Z
                        };
                        for (uint c=0;c<3;c++)
                        {
                            if (!RSolverStress::isComponentEnabled(bc,components[c]))
                            {
                                continue;
                            }
                            RNodeConstraint constraint;
                            constraint.direction = RR3Vector(0.0,0.0,0.0);
                            constraint.direction[c] = 1.0;
                            constraint.value = this->findComponentValue(bc,components[c]);
                            constraints.push_back(constraint);
                        }
                        continue;
                    }

                    double value = this->findComponentValue(bc,R_VARIABLE_DISPLACEMENT);

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT_ROLLER &&
                        entityType == R_ENTITY_GROUP_LINE)
                    {
                        // Free to slide along the line, held across it. The
                        // prescribed value has no single direction to act in
                        // here, so both held directions are held at zero.
                        if (!nodeLineDirectionSet[nodeId])
                        {
                            continue;
                        }
                        RRMatrix frame;
                        nodeLineDirection[nodeId].findRotationMatrix(frame);

                        for (uint c=1;c<3;c++)
                        {
                            RNodeConstraint constraint;
                            constraint.direction = RR3Vector(frame[0][c],frame[1][c],frame[2][c]);
                            constraint.value = 0.0;
                            constraints.push_back(constraint);
                        }
                        continue;
                    }

                    if (!nodeNormalSet[nodeId])
                    {
                        continue;
                    }

                    RNodeConstraint constraint;
                    constraint.direction = nodeNormal[nodeId];
                    constraint.value = value;
                    constraints.push_back(constraint);

                    if (bc.getType() == R_BOUNDARY_CONDITION_DISPLACEMENT_NORMAL)
                    {
                        // Held in the two directions across the normal as well.
                        RRMatrix frame;
                        nodeNormal[nodeId].findRotationMatrix(frame);

                        for (uint c=1;c<3;c++)
                        {
                            RNodeConstraint tangential;
                            tangential.direction = RR3Vector(frame[0][c],frame[1][c],frame[2][c]);
                            tangential.value = 0.0;
                            constraints.push_back(tangential);
                        }
                    }
                }
            }
        }
    }

    // Reduce the constraints of each node to an orthonormal set and build the
    // node frame from it.
    for (uint i=0;i<nNodes;i++)
    {
        const std::vector<RNodeConstraint> &constraints = nodeConstraints[i];
        if (constraints.empty())
        {
            this->localRotations[i].deactivate();
            continue;
        }

        RR3Vector axis[3];
        double prescribed[3] = { 0.0, 0.0, 0.0 };
        uint nConstrained = 0;

        double valueScale = 0.0;
        for (uint j=0;j<constraints.size();j++)
        {
            valueScale = std::max(valueScale,std::fabs(constraints[j].value));
        }

        // Every record is examined, including those beyond the third one - a
        // direction which is already held may still disagree about its value.
        for (uint j=0;j<constraints.size();j++)
        {
            RR3Vector w(constraints[j].direction);
            double wValue = constraints[j].value;

            // Gram-Schmidt against what has been accepted, carrying the
            // prescribed values through the same operations.
            for (uint l=0;l<nConstrained;l++)
            {
                double overlap = RR3Vector::dot(constraints[j].direction,axis[l]);
                w[0] -= overlap*axis[l][0];
                w[1] -= overlap*axis[l][1];
                w[2] -= overlap*axis[l][2];
                wValue -= overlap*prescribed[l];
            }

            double norm = w.length();
            if (norm > 1.0e-8 && nConstrained < 3)
            {
                axis[nConstrained] = RR3Vector(w[0]/norm,w[1]/norm,w[2]/norm);
                prescribed[nConstrained] = wValue/norm;
                nConstrained++;
            }
            else if (std::fabs(wValue) > 1.0e-8*std::max(valueScale,1.0e-12))
            {
                throw RError(RError::Type::Application,R_ERROR_REF,
                             "Node %u carries displacement constraints which contradict each other. "
                             "Check that the entities meeting at this node do not prescribe different "
                             "displacements in the same direction.",i);
            }
        }

        // Everything accepted so far is held; what follows only completes the
        // frame and stays free.
        uint nHeld = nConstrained;

        for (uint c=0;c<3 && nConstrained<3;c++)
        {
            RR3Vector candidate(0.0,0.0,0.0);
            candidate[c] = 1.0;

            for (uint l=0;l<nConstrained;l++)
            {
                double overlap = RR3Vector::dot(candidate,axis[l]);
                candidate[0] -= overlap*axis[l][0];
                candidate[1] -= overlap*axis[l][1];
                candidate[2] -= overlap*axis[l][2];
            }

            double norm = candidate.length();
            if (norm > 1.0e-8)
            {
                axis[nConstrained] = RR3Vector(candidate[0]/norm,candidate[1]/norm,candidate[2]/norm);
                prescribed[nConstrained] = 0.0;
                nConstrained++;
            }
        }

        this->nodeConstrainedDirections[i] = nHeld;

        RRMatrix R(3,3,0.0);
        for (uint r=0;r<3;r++)
        {
            for (uint c=0;c<3;c++)
            {
                R[r][c] = axis[c][r];
            }
        }

        for (uint c=0;c<3;c++)
        {
            this->nodePrescribedDisplacement[i][c] = prescribed[c];
        }

        // A frame which is already the global one needs no rotation at all.
        bool isIdentity = true;
        for (uint r=0;r<3 && isIdentity;r++)
        {
            for (uint c=0;c<3;c++)
            {
                double expected = (r == c) ? 1.0 : 0.0;
                if (std::fabs(R[r][c]-expected) > 1.0e-12)
                {
                    isIdentity = false;
                    break;
                }
            }
        }

        if (isIdentity)
        {
            this->localRotations[i].deactivate();
        }
        else
        {
            this->localRotations[i].activate(R);
        }
    }

    // Rotation conditions use explicit global X/Y/Z components. Keep their
    // equation frame separate from translation-only normal/roller frames so a
    // global rotational restraint cannot silently rotate with a surface normal.
    // Nodes without an explicit rotation condition retain the historical frame
    // behaviour used by the implicit fully-fixed clamp fallback.
    for (uint i=0;i<nNodes;i++)
    {
        this->rotationLocalRotations[i] = this->localRotations[i];
    }

    const RVariableType rotationComponents[3] =
    {
        R_VARIABLE_ROTATION_X,
        R_VARIABLE_ROTATION_Y,
        R_VARIABLE_ROTATION_Z
    };
    for (uint i=0;i<this->pModel->getNElementGroups();i++)
    {
        const RElementGroup *pElementGroup = this->pModel->getElementGroupPtr(i);
        if (!pElementGroup)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Element group could not be found (%u of %u).",i,this->pModel->getNElementGroups());
        }
        for (uint j=0;j<pElementGroup->getNBoundaryConditions();j++)
        {
            const RBoundaryCondition &bc = pElementGroup->getBoundaryCondition(j);
            if (bc.getType() != R_BOUNDARY_CONDITION_ROTATION)
            {
                continue;
            }
            for (uint k=0;k<pElementGroup->size();k++)
            {
                const RElement &element = this->pModel->getElement(pElementGroup->get(k));
                for (uint l=0;l<element.size();l++)
                {
                    const uint nodeId = element.getNodeId(l);
                    this->nodeHasRotationCondition[nodeId] = true;
                    this->rotationLocalRotations[nodeId].deactivate();
                    for (uint c=0;c<3;c++)
                    {
                        if (!RSolverStress::isComponentEnabled(bc,rotationComponents[c]))
                        {
                            continue;
                        }
                        const double value = this->findComponentValue(bc,rotationComponents[c]);
                        const uint bit = 1u << c;
                        if ((this->nodeConstrainedRotations[nodeId] & bit) != 0u &&
                            std::fabs(this->nodePrescribedRotation[nodeId][c]-value) >
                                1.0e-8*std::max({std::fabs(value),
                                                std::fabs(this->nodePrescribedRotation[nodeId][c]),
                                                1.0e-12}))
                        {
                            throw RError(RError::Type::Application,R_ERROR_REF,
                                         "Node %u carries rotation constraints which contradict each other.",nodeId);
                        }
                        this->nodeConstrainedRotations[nodeId] |= bit;
                        this->nodePrescribedRotation[nodeId][c] = value;
                    }
                }
            }
        }
    }
}

void RSolverStress::generateEquationBook()
{
    this->equationBook.resize(this->pModel->getNNodes());
    for (uint i=0;i<this->pModel->getNElements();i++)
    {
        if (this->computableElements[i])
        {
            const RElement &rElement = this->pModel->getElement(i);
            const QList<RDofType> &dofs = RElementDescriptor::get(rElement.getType()).getDofsPerNode();
            for (uint j=0;j<rElement.size();j++)
            {
                for (RDofType dof : dofs)
                {
                    this->equationBook.enable(rElement.getNodeId(j),dof);
                }
            }
        }
    }

    // The held translational directions of a node are the leading directions
    // of its local constraint frame.
    for (uint i=0;i<this->pModel->getNNodes();i++)
    {
        for (uint j=0;j<this->nodeConstrainedDirections[i];j++)
        {
            this->equationBook.disable(i,translationDof(j));
        }
        if (this->nodeHasRotationCondition[i])
        {
            for (uint c=0;c<3;c++)
            {
                if ((this->nodeConstrainedRotations[i] & (1u << c)) != 0u)
                {
                    const RDofType dof = c == 0 ? RDofType::RX
                                        : c == 1 ? RDofType::RY
                                                 : RDofType::RZ;
                    this->equationBook.disable(i,dof);
                }
            }
        }
        // Old models have no rotation condition. Preserve their established
        // meaning: a node fixed in all three translations is a clamped frame
        // node. Adding an explicit Rotation condition opts out of this fallback.
        else if (this->nodeConstrainedDirections[i] == 3)
        {
            this->equationBook.disable(i,RDofType::RX);
            this->equationBook.disable(i,RDofType::RY);
            this->equationBook.disable(i,RDofType::RZ);
        }
    }
}

void RSolverStress::applyLinearConstraints()
{
    this->nPhysicalEquations = this->equationBook.getNEnabled();
    this->constraintNullspace.clear();
    this->constraintParticular.clear();
    this->constraintFreeEquations.clear();
    const RStructuralInteractionSetup &interactionSetup = this->pModel->getProblemSetup()
        .getStructuralInteractionSetup();
    const std::vector<RLinearConstraint> constraints = interactionSetup.buildLinearConstraints(*this->pModel);
    if (constraints.empty())
    {
        return;
    }

    const uint nPhysicalEquations = this->equationBook.getNEnabled();
    struct ReducedEquation
    {
        RRVector coefficients;
        double rhs = 0.0;
        uint pivot = 0;
        QString id;
    };
    std::vector<ReducedEquation> independent;

    for (uint constraintIndex=0;constraintIndex<constraints.size();constraintIndex++)
    {
        const RLinearConstraint &constraint = constraints[constraintIndex];
        // Explicit equations are model data in physical units. Generated
        // couplings are constructed from already scaled node coordinates.
        const bool explicitEquation = constraintIndex < interactionSetup.getLinearConstraints().size();
        ReducedEquation reduced;
        reduced.coefficients.resize(nPhysicalEquations,0.0);
        reduced.coefficients.fill(0.0);
        reduced.rhs = constraint.rhs;
        reduced.id = constraint.id;

        for (const RLinearConstraintTerm &term : constraint.terms)
        {
            const uint nodeId = term.reference.nodeId;
            uint component = 0;
            const bool translation = findTranslationComponent(term.reference.dof,component);
            const bool rotation = !translation && findRotationComponent(term.reference.dof,component);
            if (!translation && !rotation)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Linear constraint '%s' contains a non-mechanical DOF.",
                             constraint.id.toUtf8().constData());
            }

            const RLocalRotation &frame = translation ? this->localRotations[nodeId]
                                                       : this->rotationLocalRotations[nodeId];
            for (uint localComponent=0;localComponent<3;localComponent++)
            {
                const double transform = frame.isActive()
                                       ? frame.getR()[component][localComponent]
                                       : (component == localComponent ? 1.0 : 0.0);
                const double dofScale = translation
                    ? this->scales.findScaleFactor(R_VARIABLE_DISPLACEMENT) : 1.0;
                const double coefficient = term.coefficient*transform /
                    (explicitEquation ? dofScale : 1.0);
                if (std::fabs(coefficient) <= 1.0e-14)
                {
                    continue;
                }
                const RDofType localDof = translation ? translationDof(localComponent)
                    : (localComponent == 0 ? RDofType::RX : localComponent == 1 ? RDofType::RY : RDofType::RZ);
                uint equation = 0;
                if (this->equationBook.getValue(nodeId,localDof,equation))
                {
                    reduced.coefficients[equation] += coefficient;
                    continue;
                }

                bool prescribed = false;
                double prescribedValue = 0.0;
                if (translation && localComponent < this->nodeConstrainedDirections[nodeId])
                {
                    prescribed = true;
                    prescribedValue = this->nodePrescribedDisplacement[nodeId][localComponent];
                }
                else if (rotation)
                {
                    prescribed = (this->nodeHasRotationCondition[nodeId] &&
                                  (this->nodeConstrainedRotations[nodeId] & (1u << localComponent)) != 0u) ||
                                 (!this->nodeHasRotationCondition[nodeId] && this->nodeConstrainedDirections[nodeId] == 3);
                    prescribedValue = this->nodePrescribedRotation[nodeId][localComponent];
                }
                if (!prescribed)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Linear constraint '%s' references inactive DOF '%s' on node %u.",
                                 constraint.id.toUtf8().constData(),
                                 RDofDescriptor::get(localDof).getId().toUtf8().constData(),nodeId);
                }
                reduced.rhs -= coefficient*prescribedValue;
            }
        }

        // Incremental reduced-row echelon form identifies redundant equations
        // and contradictory equations before an indefinite saddle system is built.
        for (const ReducedEquation &basis : independent)
        {
            const double factor = reduced.coefficients[basis.pivot];
            if (std::fabs(factor) <= 1.0e-12)
            {
                continue;
            }
            for (uint column=0;column<nPhysicalEquations;column++)
            {
                reduced.coefficients[column] -= factor*basis.coefficients[column];
            }
            reduced.rhs -= factor*basis.rhs;
        }

        double rowScale = 0.0;
        uint pivot = nPhysicalEquations;
        for (uint column=0;column<nPhysicalEquations;column++)
        {
            rowScale = std::max(rowScale,std::fabs(reduced.coefficients[column]));
            if (pivot == nPhysicalEquations && std::fabs(reduced.coefficients[column]) > 1.0e-10)
            {
                pivot = column;
            }
        }
        if (pivot == nPhysicalEquations)
        {
            if (std::fabs(reduced.rhs) > 1.0e-10*std::max(1.0,rowScale))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Linear constraint '%s' contradicts prescribed or preceding constraints.",
                             constraint.id.toUtf8().constData());
            }
            RLogger::warning("Skipping redundant linear constraint '%s'.\n",constraint.id.toUtf8().constData());
            continue;
        }

        const double divisor = reduced.coefficients[pivot];
        for (uint column=0;column<nPhysicalEquations;column++)
        {
            reduced.coefficients[column] /= divisor;
        }
        reduced.rhs /= divisor;
        reduced.pivot = pivot;

        for (ReducedEquation &basis : independent)
        {
            const double factor = basis.coefficients[pivot];
            if (std::fabs(factor) <= 1.0e-12)
            {
                continue;
            }
            for (uint column=0;column<nPhysicalEquations;column++)
            {
                basis.coefficients[column] -= factor*reduced.coefficients[column];
            }
            basis.rhs -= factor*reduced.rhs;
        }
        independent.push_back(reduced);
        std::sort(independent.begin(),independent.end(),[](const ReducedEquation &a, const ReducedEquation &b)
        {
            return a.pivot < b.pivot;
        });
    }

    this->nLinearConstraintEquations = uint(independent.size());
    const bool eliminate = this->problemType == R_PROBLEM_STRESS_MODAL ||
        (this->problemType == R_PROBLEM_STRESS &&
         this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled());
    if (eliminate)
    {
        if (this->problemType == R_PROBLEM_STRESS_MODAL)
        {
            for (const ReducedEquation &row : independent)
            {
                if (std::fabs(row.rhs) > 1.0e-12)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Modal constraint '%s' is nonhomogeneous; solve about an explicit equilibrium state first.",
                                 row.id.toUtf8().constData());
                }
            }
        }
        const uint nReduced = nPhysicalEquations-this->nLinearConstraintEquations;
        if (nReduced == 0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "All structural degrees of freedom were eliminated by linear constraints.");
        }
        this->constraintNullspace.resize(nPhysicalEquations,nReduced,0.0);
        this->constraintNullspace.fill(0.0);
        this->constraintParticular.resize(nPhysicalEquations,0.0);
        this->constraintParticular.fill(0.0);
        std::vector<bool> dependent(nPhysicalEquations,false);
        for (const ReducedEquation &row : independent)
        {
            dependent[row.pivot] = true;
            this->constraintParticular[row.pivot] = row.rhs;
        }
        uint reducedColumn = 0;
        for (uint column=0;column<nPhysicalEquations;column++)
        {
            if (dependent[column])
            {
                continue;
            }
            this->constraintFreeEquations.push_back(column);
            this->constraintNullspace[column][reducedColumn] = 1.0;
            for (const ReducedEquation &row : independent)
            {
                this->constraintNullspace[row.pivot][reducedColumn] = -row.coefficients[column];
            }
            reducedColumn++;
        }

        std::vector<std::vector<std::pair<uint,double>>> sparseNullspace(nPhysicalEquations);
        for (uint row=0;row<nPhysicalEquations;row++)
        {
            for (uint column=0;column<nReduced;column++)
            {
                const double value = this->constraintNullspace[row][column];
                if (std::fabs(value) > 1.0e-14)
                {
                    sparseNullspace[row].push_back({column,value});
                }
            }
        }
        auto reduceMatrix = [&](const RSparseMatrix &original)
        {
            RSparseMatrix reduced;
            reduced.setNRows(nReduced);
            for (uint row=0;row<nPhysicalEquations;row++)
            {
                for (uint position=0;position<original.getNColumns(row);position++)
                {
                    const uint column = original.getVector(row).getIndex(position);
                    const double value = original.getValue(row,position);
                    for (const auto &left : sparseNullspace[row])
                    {
                        for (const auto &right : sparseNullspace[column])
                        {
                            reduced.addValue(left.first,right.first,left.second*value*right.second);
                        }
                    }
                }
            }
            return reduced;
        };
        RRVector adjustedForce(this->b);
        if (this->problemType != R_PROBLEM_STRESS_MODAL)
        {
            RRVector offsetStiffness;
            RSparseMatrix::mlt(this->A,this->constraintParticular,offsetStiffness);
            for (uint row=0;row<nPhysicalEquations;row++)
            {
                adjustedForce[row] -= offsetStiffness[row];
            }
        }
        RRVector reducedForce(nReduced,0.0);
        for (uint row=0;row<nPhysicalEquations;row++)
        {
            for (const auto &term : sparseNullspace[row])
            {
                reducedForce[term.first] += term.second*adjustedForce[row];
            }
        }
        this->A = reduceMatrix(this->A);
        this->M = reduceMatrix(this->M);
        this->b = reducedForce;
        this->x.resize(nReduced,0.0);
        this->nLinearConstraintEquations = 0;
        return;
    }
    const uint totalEquations = nPhysicalEquations + this->nLinearConstraintEquations;
    double constraintScale = 1.0;
    for (uint row=0;row<nPhysicalEquations;row++)
    {
        constraintScale = std::max(constraintScale,std::fabs(this->A.findValue(row,row)));
    }
    this->A.setNRows(totalEquations);
    this->b.resize(totalEquations,0.0);
    this->x.resize(totalEquations,0.0);
    for (uint i=0;i<this->nLinearConstraintEquations;i++)
    {
        const uint multiplier = nPhysicalEquations+i;
        for (uint column=0;column<nPhysicalEquations;column++)
        {
            const double coefficient = constraintScale*independent[i].coefficients[column];
            if (std::fabs(coefficient) <= 1.0e-14)
            {
                continue;
            }
            this->A.addValue(multiplier,column,coefficient);
            this->A.addValue(column,multiplier,coefficient);
        }
        this->b[multiplier] = constraintScale*independent[i].rhs;
    }
}

RRVector RSolverStress::expandReducedSolution(const RRVector &reduced, bool includeParticular) const
{
    if (this->constraintNullspace.getNRows() == 0)
    {
        return reduced;
    }
    RRVector physical(this->nPhysicalEquations,0.0);
    for (uint row=0;row<this->nPhysicalEquations;row++)
    {
        physical[row] = includeParticular ? this->constraintParticular[row] : 0.0;
        for (uint column=0;column<reduced.size();column++)
        {
            physical[row] += this->constraintNullspace[row][column]*reduced[column];
        }
    }
    return physical;
}

RRVector RSolverStress::projectPhysicalState(const RRVector &physical, bool subtractParticular) const
{
    if (this->constraintNullspace.getNRows() == 0)
    {
        return physical;
    }
    RRVector reduced(this->constraintNullspace.getNColumns(),0.0);
    for (uint column=0;column<this->constraintFreeEquations.size();column++)
    {
        const uint row = this->constraintFreeEquations[column];
        reduced[column] = physical[row] - (subtractParticular ? this->constraintParticular[row] : 0.0);
    }
    return reduced;
}

void RSolverStress::assemblyMatrix(uint elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp, RSparseMatrix &Mp, RSparseMatrix &Kfc, RSparseMatrix &Mfc)
{
    double alpha = this->pModel->getTimeSolver().getTimeMarchApproximationCoefficient();
    double dt = this->pModel->getTimeSolver().getCurrentTimeStepSize();

    const RElement &rElement = this->pModel->getElement(elementID);
    const QList<RDofType> &elementDofs = RElementDescriptor::get(rElement.getType()).getDofsPerNode();
    const uint nDofs = uint(elementDofs.size());
    const uint nElementDofs = nDofs*rElement.size();

    R_ERROR_ASSERT(Me.getNRows() == nElementDofs && Me.getNColumns() == nElementDofs);
    R_ERROR_ASSERT(Ke.getNRows() == nElementDofs && Ke.getNColumns() == nElementDofs);
    R_ERROR_ASSERT(fe.getNRows() == nElementDofs);

    RRMatrix Ae(nElementDofs,nElementDofs);
    RRMatrix Be(nElementDofs,nElementDofs);
    RRVector be(nElementDofs);

    Ae.fill(0.0);
    Be.fill(0.0);
    be.fill(0.0);

    const bool structuralDynamics = this->problemType == R_PROBLEM_STRESS &&
        this->pModel->getProblemSetup().getStructuralDynamicsSetup().getEnabled();
    if (this->pModel->getTimeSolver().getEnabled() && !structuralDynamics)
    {
        for (uint m=0;m<rElement.size();m++)
        {
            for (uint i=0;i<nDofs;i++)
            {
                const uint row = nDofs*m+i;
                be[row] = dt*fe[row];
                for (uint n=0;n<rElement.size();n++)
                {
                    for (uint j=0;j<nDofs;j++)
                    {
                        const uint column = nDofs*n+j;
                        Ae[row][column] = Me[row][column] + alpha*dt*Ke[row][column];

                        uint component = 0;
                        double oldValue = 0.0;
                        if (findTranslationComponent(elementDofs[j],component))
                        {
                            const uint nodeId = rElement.getNodeId(n);
                            oldValue = component == 0 ? this->nodeDisplacement.x[nodeId]
                                     : component == 1 ? this->nodeDisplacement.y[nodeId]
                                                      : this->nodeDisplacement.z[nodeId];
                        }
                        else if (findRotationComponent(elementDofs[j],component))
                        {
                            const uint nodeId = rElement.getNodeId(n);
                            oldValue = component == 0 ? this->nodeRotation.x[nodeId]
                                     : component == 1 ? this->nodeRotation.y[nodeId]
                                                      : this->nodeRotation.z[nodeId];
                        }
                        be[row] += (Me[row][column] - (1.0-alpha)*dt*Ke[row][column])*oldValue;
                    }
                }
            }
        }
    }
    else if (this->problemType == R_PROBLEM_STRESS_MODAL || structuralDynamics)
    {
        Ae = Ke;
        Be = Me;
        be = fe;
    }
    else
    {
        Ae = Ke;
        be = fe;
    }
    this->applyLocalRotations(elementID,Ae);
    this->applyLocalRotations(elementID,Be);
    this->applyLocalRotations(elementID,be);
    // Save the actual applied load before held-column elimination rewrites be.
    const RRVector rawElementForce(be);

    // Apply explicit boundary conditions. Ae and be are in the frame of each
    // node by now, so the prescribed values have to be taken in that same
    // frame - which is what nodePrescribedDisplacement holds. The eigen-value
    // problem is homogeneous, so nothing is prescribed there.
    if (this->problemType != R_PROBLEM_STRESS_MODAL)
    {
        for (uint m=0;m<rElement.size();m++)
        {
            uint position;
            uint nodeID = rElement.getNodeId(m);
            for (uint c=0;c<nDofs;c++)
            {
                if (this->equationBook.getValue(nodeID,elementDofs[c],position))
                {
                    continue;
                }
                uint component = 0;
                double u = 0.0;
                if (findTranslationComponent(elementDofs[c],component))
                {
                    u = this->nodePrescribedDisplacement[nodeID][component];
                }
                else if (findRotationComponent(elementDofs[c],component))
                {
                    u = this->nodePrescribedRotation[nodeID][component];
                }
                if (std::fabs(u) <= 0.0)
                {
                    continue;
                }
                for (uint n=0;n<nElementDofs;n++)
                {
                    be[n] -= Ae[n][nDofs*m+c]*u;
                }
            }
        }
    }

    // Assembly final matrix system
    for (uint m=0;m<rElement.size();m++)
    {
        for (uint i=0;i<nDofs;i++)
        {
            uint mp = 0;

            if (structuralDynamics && this->supportResponseObserver)
            {
                const uint node = rElement.getNodeId(m);
                const uint dof = uint(elementDofs[i]);
                const int heldRow = dof < 6 ? this->heldIndexByNodeDof[6*node+dof] : -1;
                if (heldRow >= 0)
                {
                    const uint row = uint(heldRow);
                    const uint thread = uint(omp_get_thread_num());
                    this->supportForceBuffers[thread][row] += rawElementForce[nDofs*m+i];
                    for (uint n=0;n<rElement.size();n++)
                    {
                        for (uint j=0;j<nDofs;j++)
                        {
                            const uint columnNode = rElement.getNodeId(n);
                            const RDofType columnDof = elementDofs[j];
                            uint column = 0;
                            if (!this->equationBook.getValue(columnNode,columnDof,column))
                            {
                                const uint dofIndex = uint(columnDof);
                                const int heldColumn = dofIndex < 6
                                    ? this->heldIndexByNodeDof[6*columnNode+dofIndex] : -1;
                                if (heldColumn < 0)
                                {
                                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                                 "Support reaction row references an inactive non-support DOF.");
                                }
                                column = this->equationBook.getNEnabled()+uint(heldColumn);
                            }
                            this->supportStiffnessBuffers[thread].addValue(
                                row,column,Ae[nDofs*m+i][nDofs*n+j]);
                            this->supportMassBuffers[thread].addValue(
                                row,column,Be[nDofs*m+i][nDofs*n+j]);
                        }
                    }
                }
            }

            if (this->equationBook.getValue(rElement.getNodeId(m),elementDofs[i],mp))
            {
                bp[mp] += be[nDofs*m+i];

                for (uint n=0;n<rElement.size();n++)
                {
                    for (uint j=0;j<nDofs;j++)
                    {
                        uint np = 0;

                        if (this->equationBook.getValue(rElement.getNodeId(n),elementDofs[j],np))
                        {
                            Ap.addValue(mp,np,Ae[nDofs*m+i][nDofs*n+j]);
                            if (this->problemType == R_PROBLEM_STRESS_MODAL || structuralDynamics)
                            {
                                // Be carries the mass in the node frames.
                                Mp.addValue(mp,np,Be[nDofs*m+i][nDofs*n+j]);
                            }
                        }
                        else if (structuralDynamics && !this->prescribedMotions.empty())
                        {
                            uint component = 0;
                            if (findTranslationComponent(elementDofs[j],component))
                            {
                                const uint node = rElement.getNodeId(n);
                                const int motion = this->motionIndexByNodeDof[3*node+component];
                                if (motion >= 0)
                                {
                                    Kfc.addValue(mp,uint(motion),Ae[nDofs*m+i][nDofs*n+j]);
                                    Mfc.addValue(mp,uint(motion),Be[nDofs*m+i][nDofs*n+j]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void RSolverStress::applyLocalRotations(unsigned int elementID, RRMatrix &Ae)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    const QList<RDofType> &dofs = RElementDescriptor::get(rElement.getType()).getDofsPerNode();
    const uint nDofs = uint(dofs.size());
    RRMatrix T;

    bool first = true;

    for (uint i=0;i<rElement.size();i++)
    {
        uint nodeId = rElement.getNodeId(i);
        const bool translationActive = this->localRotations[nodeId].isActive();
        const bool rotationActive = this->rotationLocalRotations[nodeId].isActive();
        if (translationActive || rotationActive)
        {
            if (first)
            {
                T.setIdentity(Ae.getNRows());
                first = false;
            }
            const int translationOffset = dofs.indexOf(RDofType::UX);
            if (translationActive && translationOffset >= 0 && dofs.indexOf(RDofType::UY) == translationOffset+1 && dofs.indexOf(RDofType::UZ) == translationOffset+2)
            {
                T.setBlock(this->localRotations[nodeId].getR(),nDofs*i+uint(translationOffset),nDofs*i+uint(translationOffset));
            }
            const int rotationOffset = dofs.indexOf(RDofType::RX);
            if (rotationActive && rotationOffset >= 0 && dofs.indexOf(RDofType::RY) == rotationOffset+1 && dofs.indexOf(RDofType::RZ) == rotationOffset+2)
            {
                T.setBlock(this->rotationLocalRotations[nodeId].getR(),nDofs*i+uint(rotationOffset),nDofs*i+uint(rotationOffset));
            }
        }
    }
    if (!first)
    {
        RRMatrix Tt(T);
        Tt.transpose();

        RRMatrix Aetmp;

        RRMatrix::mlt(Tt,Ae,Aetmp);
        RRMatrix::mlt(Aetmp,T,Ae);
    }
}

void RSolverStress::applyLocalRotations(unsigned int elementID, RRVector &fe)
{
    const RElement &rElement = this->pModel->getElement(elementID);
    const QList<RDofType> &dofs = RElementDescriptor::get(rElement.getType()).getDofsPerNode();
    const uint nDofs = uint(dofs.size());
    RRMatrix T;

    bool first = true;

    for (uint i=0;i<rElement.size();i++)
    {
        uint nodeId = rElement.getNodeId(i);
        const bool translationActive = this->localRotations[nodeId].isActive();
        const bool rotationActive = this->rotationLocalRotations[nodeId].isActive();
        if (translationActive || rotationActive)
        {
            if (first)
            {
                T.setIdentity(fe.getNRows());
                first = false;
            }
            const int translationOffset = dofs.indexOf(RDofType::UX);
            if (translationActive && translationOffset >= 0 && dofs.indexOf(RDofType::UY) == translationOffset+1 && dofs.indexOf(RDofType::UZ) == translationOffset+2)
            {
                T.setBlock(this->localRotations[nodeId].getR(),nDofs*i+uint(translationOffset),nDofs*i+uint(translationOffset));
            }
            const int rotationOffset = dofs.indexOf(RDofType::RX);
            if (rotationActive && rotationOffset >= 0 && dofs.indexOf(RDofType::RY) == rotationOffset+1 && dofs.indexOf(RDofType::RZ) == rotationOffset+2)
            {
                T.setBlock(this->rotationLocalRotations[nodeId].getR(),nDofs*i+uint(rotationOffset),nDofs*i+uint(rotationOffset));
            }
        }
    }
    if (!first)
    {
        // T maps local to global, so with u_global = T * u_local the system
        // K_local = T^T * K_global * T has to be paired with the load vector
        // f_local = T^T * f_global. Rotating the load with T instead leaves the
        // solution satisfying equilibrium with a rotated load.
        RRMatrix Tt(T);
        Tt.transpose();

        RRVector fetmp;

        RRMatrix::mlt(Tt,fe,fetmp);
        fe = fetmp;
    }
}
