#ifndef RSOLVERSTRESS_H
#define RSOLVERSTRESS_H

#include <rbl_rvector.h>
#include <rml_structural_interaction.h>
#include <functional>
#include <memory>
#include <QString>

#include "requationbook.h"
#include "rstructural_support_response.h"
#include "rstructural_control_protocol.h"
#include "rstructural_control_conditioner.h"
#include "rstructural_cable_model_adapter.h"
#include "rstructural_cable_dynamics.h"
#include "rstructural_piezo_model_adapter.h"
#include "rsolvergeneric.h"

class RStructuralDynamicsIntegrator;
struct RStructuralDynamicsState;

class RSolverStress : public RSolverGeneric
{

    protected:

        //! Element modulus of elasticity vector.
        RRVector elementElasticityModulus;
        //! Element Poisson ratio.
        RRVector elementPoissonRatio;
        //! Element density vector.
        RRVector elementDensity;
        //! Element thermal expansion coefficient vector.
        RRVector elementThermalExpansion;
        //! Element environment temperature vector.
        RRVector elementEnvironmentTemperature;
        //! Number of constrained directions of each node, 0 to 3. They are the
        //! first directions of the node local frame held in localRotations.
        RUVector nodeConstrainedDirections;
        //! Prescribed displacement of each node, expressed in the frame of that
        //! node - the local frame where one is active, global otherwise. Only
        //! the first nodeConstrainedDirections components carry a meaning.
        RRMatrix nodePrescribedDisplacement;
        //! Rotational equation frame. It normally follows the translational
        //! frame for backward compatibility, but stays global when an explicit
        //! global rotation condition is present on the node.
        std::vector<RLocalRotation> rotationLocalRotations;
        //! Bit mask of explicitly constrained global rotations (RX/RY/RZ).
        RUVector nodeConstrainedRotations;
        //! Whether a rotation boundary condition was explicitly attached.
        RBVector nodeHasRotationCondition;
        //! Prescribed global nodal rotations.
        RRMatrix nodePrescribedRotation;
        //! Node displacement vector.
        RSolverCartesianVector<RRVector> nodeDisplacement;
        //! Node rotation vector used by frame elements.
        RSolverCartesianVector<RRVector> nodeRotation;
        //! Translational velocity/acceleration in global coordinates.
        RSolverCartesianVector<RRVector> nodeVelocity;
        RSolverCartesianVector<RRVector> nodeAcceleration;
        //! Rotational velocity/acceleration in global coordinates (frame nodes).
        RSolverCartesianVector<RRVector> nodeAngularVelocity;
        RSolverCartesianVector<RRVector> nodeAngularAcceleration;
        //! Node initial displacement vector.
        RSolverCartesianVector<RRVector> nodeInitialDisplacement;
        //! Node force vector.
        RSolverCartesianVector<RRVector> nodeForce;
        //! Node moment vector used by frame elements.
        RSolverCartesianVector<RRVector> nodeMoment;
        //! Node pressure.
        RRVector nodePressure;
        //! Element stress components.
        //! Volume elements store them in global coordinates, surface and line
        //! elements in their own local element frame.
        //! Order: xx, yy, zz, yz, xz, xy.
        RRVector elementStress[6];
        //! Element normal stress.
        RRVector elementNormalStress;
        //! Element shear stress.
        RRVector elementShearStress;
        //! Element VonMisses stress.
        RRVector elementVonMisses;

        //! Eigen values.
        RRVector d;
        //! Eigen vectors.
        RRMatrix ev;

        //! Structural node/DOF to equation mapping.
        REquationBook equationBook;
        //! Number of Lagrange multiplier equations appended for MPCs.
        uint nLinearConstraintEquations;
        //! q = constraintParticular + constraintNullspace * reduced q.
        RRMatrix constraintNullspace;
        RRVector constraintParticular;
        std::vector<uint> constraintFreeEquations;
        uint nPhysicalEquations = 0;
        std::unique_ptr<RStructuralDynamicsIntegrator> dynamicsIntegrator;
        RRVector previousDynamicForce;
        RRVector currentDynamicForce;
        RRVector initialDynamicForce;
        bool preparingInitialForce = false;
        //! Linear dynamics reuses matrices; reject unsupported time variation.
        RSparseMatrix dynamicMassReference;
        RSparseMatrix dynamicStiffnessReference;
        //! Free-to-prescribed support blocks (rows are free equations).
        RSparseMatrix motionStiffness;
        RSparseMatrix motionMass;
        RSparseMatrix dynamicMotionStiffnessReference;
        RSparseMatrix dynamicMotionMassReference;
        std::vector<RStructuralPrescribedMotion> prescribedMotions;
        std::vector<int> motionIndexByNodeDof;
        QString dynamicMotionDefinitionReference;
        RRMatrix dynamicConstraintNullspaceReference;
        RRVector dynamicConstraintParticularReference;
        std::vector<double> dynamicBoundaryReference;
        bool dynamicBoundaryReferenceReady = false;
        double previousDynamicTime = 0.0;
        double currentDynamicTime = 0.0;
        struct DynamicsStepSnapshot;
        std::unique_ptr<DynamicsStepSnapshot> dynamicsStepSnapshot;
        std::function<void(const std::vector<RStructuralSupportSample> &)> supportResponseObserver;
        std::vector<RStructuralSupportSample> acceptedSupportSamples;
        std::vector<RStructuralSupportSample> trialSupportSamples;
        QString controllerId;
        QString controllerRevision;
        std::function<RStructuralControlDecision(const RStructuralControlObservation &)> controlCallback;
        std::function<void(const RStructuralControlStepRecord &)> controlStepObserver;
        std::unique_ptr<RStructuralControlForceMap> controlForceMap;
        std::unique_ptr<RStructuralControlLedger> controlLedger;
        std::unique_ptr<RStructuralControlConditioner> controlConditioner;
        RStructuralControlConditioner::State acceptedConditionerState;
        RStructuralControlConditioner::State candidateConditionerState;
        QString cachedControlFingerprint;
        QString cachedControlCatalogFingerprint;
        std::vector<RStructuralControlCommand> activeControlCommands;
        RRVector activeControlForce;
        QByteArray acceptedControllerState;
        QByteArray candidateControllerState;
        RStructuralControlStepRecord trialControlRecord;
        QString pendingControlObserverFailure;
        //! Last converged result from the dedicated static cable equilibrium path.
        RStructuralCableEquilibriumResult structuralCableResult;
        bool structuralCableResultReady = false;
        //! Last candidate/accepted result from the cable Newmark FE path.
        RStructuralCableDynamicsResult structuralCableDynamicsResult;
        bool structuralCableDynamicsResultReady = false;
        bool structuralCableModeAtConstruction = false;
        bool structuralCableDynamicsModeAtConstruction = false;
        RStructuralPiezoModelResult structuralPiezoResult;
        bool structuralPiezoResultReady = false;
        bool structuralPiezoModeAtConstruction = false;
        struct ControlSensorMap
        {
            QString id;
            RStructuralSensorQuantity quantity;
            QString unit;
            std::vector<std::pair<uint,double>> weights;
            double displacementOffset = 0.0;
        };
        std::vector<ControlSensorMap> controlSensors;
        std::vector<RStructuralDofRef> heldDofs;
        std::vector<int> heldIndexByNodeDof;
        std::vector<bool> pointSupportMask;
        RSparseMatrix supportStiffness;
        RSparseMatrix supportMass;
        RSparseMatrix supportStiffnessReference;
        RSparseMatrix supportMassReference;
        RRVector supportForce;
        RRVector initialSupportForce;
        std::vector<RSparseMatrix> supportStiffnessBuffers;
        std::vector<RSparseMatrix> supportMassBuffers;
        std::vector<RRVector> supportForceBuffers;

    public:

        //! Constructor.
        explicit RSolverStress(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData, bool modalAnalysis);

        //! Destructor.
        ~RSolverStress();

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Number of modes the last eigen-value solve actually produced.
        uint getNComputedModes() const override;

        void setSupportResponseObserver(std::function<void(const std::vector<RStructuralSupportSample> &)> observer);
        std::vector<RStructuralSupportSample> getAcceptedSupportSamples() const;
        void setControlCallback(const QString &id, const QString &revision,
                                std::function<RStructuralControlDecision(const RStructuralControlObservation &)> callback);
        void setControlStepObserver(std::function<void(const RStructuralControlStepRecord &)> observer);
        std::vector<RStructuralControlChannelSample> getAcceptedControlSamples() const;
        QString getPendingControlObserverFailure() const;

    protected:

        //! Initialize solver.
        void initialize() override;

        //! Preserve/restore the committed integrator and force history while
        //! RSolverGeneric transacts the externally visible RModel state.
        void beginStepTransaction() override;
        void rollbackStepTransaction() override;
        void commitStepTransaction() override;
        void validateRunConfiguration() const override;

        //! Update scales.
        void updateScales() override;

        //! Specialized static modes validate TRUSS1 topology and bypass the
        //! legacy material-driven computability scan.
        void findComputableElements(RProblemType problemType) override;

        //! Recover previously computed results.
        void recover() override;

        //! Prepare solver.
        void prepare() override;

        //! Run matrix solver.
        void solve() override;

        //! Solve stress-strain problem.
        void solveStressStrain();

        //! Advance the committed second-order structural state one fixed step.
        void solveStructuralDynamics();

        //! Read/write the physical equation-book state in node frames.
        RRVector gatherNodalState(const RSolverCartesianVector<RRVector> &translation,
                                  const RSolverCartesianVector<RRVector> &rotation) const;
        void scatterNodalState(const RRVector &physical,
                               RSolverCartesianVector<RRVector> &translation,
                               RSolverCartesianVector<RRVector> &rotation) const;

        //! Fail explicitly on moving prescribed DOFs (base-motion equations are future work).
        void validateDynamicBoundaryState();

        //! Solve eigen-value problem.
        void solveEigenValue();

        //! Set displacemen.
        void setDisplacement(const RRVector &v);

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        bool isStructuralCableMode() const;
        bool isStructuralCableDynamicsMode() const;
        bool isStructuralPiezoMode() const;
        bool isSpecializedStaticTrussMode() const;

        //! Local rotations are built from the constraints themselves, see
        //! generateLocalConstraints(), so the generic geometric pass is not used.
        void updateLocalRotations() override;

        //! Collect every displacement constraint acting on each node, reduce
        //! them to an orthonormal set, and from that build the node local frame,
        //! the number of constrained directions and the prescribed values.
        void generateLocalConstraints();

        //! Generate the active structural equation book.
        void generateEquationBook();

        //! Validate and append independent linear constraints using Lagrange multipliers.
        void applyLinearConstraints();

        //! Expand an eliminated-coordinate solution to the physical equation book.
        RRVector expandReducedSolution(const RRVector &reduced, bool includeParticular) const;

        //! Project physical equation-book values to the independent coordinates.
        RRVector projectPhysicalState(const RRVector &physical, bool subtractParticular) const;

        //! Physical equation order is part of the persisted dynamics checkpoint contract.
        std::vector<RStructuralDofRef> orderedPhysicalDofs() const;

        //! Stable digest of the reduced operator, constraint map and held supports.
        QString dynamicSystemFingerprint() const;
        QString supportResponseFingerprint() const;
        QString motionDefinitionFingerprint() const;
        void preparePrescribedMotions();
        RStructuralMotionState scaledMotionState(const RStructuralPrescribedMotion &motion, double time) const;
        void motionEffectiveForce(double time, RRVector &force) const;
        void setPrescribedMotionResults(double time);
        void prepareSupportResponse();
        void prepareControlMapping();
        void prepareCableControlMapping();
        void solveStructuralCableDynamics();
        QString controlMappingFingerprint() const;
        QString controlCatalogFingerprint() const;
        std::vector<RStructuralControlSensorSample> sampleControlSensors(
            const RStructuralDynamicsState &state) const;
        std::vector<RStructuralSupportSample> computeSupportSamples(
            const RStructuralDynamicsState &state, const RRVector &heldForce,
            bool initialOperator = false) const;

        //! Return true if given component of the boundary condition is present
        //! and switched on.
        static bool isComponentEnabled(const RBoundaryCondition &bc, RVariableType variableType);

        //! Return the value of given component of the boundary condition at the
        //! current solver time, or zero when the component is absent.
        double findComponentValue(const RBoundaryCondition &bc, RVariableType variableType) const;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp, RSparseMatrix &Mp, RSparseMatrix &Kfc, RSparseMatrix &Mfc);

        //! Apply local rotations to matrix.
        void applyLocalRotations(unsigned int elementID, RRMatrix &Ae);

        //! Apply local rotations to vector.
        void applyLocalRotations(unsigned int elementID, RRVector &fe);

};

#endif // RSOLVERSTRESS_H
