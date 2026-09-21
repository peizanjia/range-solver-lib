#ifndef RSOLVERSTRESS_H
#define RSOLVERSTRESS_H

#include <rbl_rvector.h>

#include "requationbook.h"
#include "rsolvergeneric.h"

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

    public:

        //! Constructor.
        explicit RSolverStress(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData, bool modalAnalysis);

        //! Destructor.
        ~RSolverStress();

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Number of modes the last eigen-value solve actually produced.
        uint getNComputedModes() const override;

    protected:

        //! Initialize solver.
        void initialize() override;

        //! Update scales.
        void updateScales() override;

        //! Recover previously computed results.
        void recover() override;

        //! Prepare solver.
        void prepare() override;

        //! Run matrix solver.
        void solve() override;

        //! Solve stress-strain problem.
        void solveStressStrain();

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

        //! Local rotations are built from the constraints themselves, see
        //! generateLocalConstraints(), so the generic geometric pass is not used.
        void updateLocalRotations() override;

        //! Collect every displacement constraint acting on each node, reduce
        //! them to an orthonormal set, and from that build the node local frame,
        //! the number of constrained directions and the prescribed values.
        void generateLocalConstraints();

        //! Generate the active structural equation book.
        void generateEquationBook();

        //! Return true if given component of the boundary condition is present
        //! and switched on.
        static bool isComponentEnabled(const RBoundaryCondition &bc, RVariableType variableType);

        //! Return the value of given component of the boundary condition at the
        //! current solver time, or zero when the component is absent.
        double findComponentValue(const RBoundaryCondition &bc, RVariableType variableType) const;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp, RSparseMatrix &Mp);

        //! Apply local rotations to matrix.
        void applyLocalRotations(unsigned int elementID, RRMatrix &Ae);

        //! Apply local rotations to vector.
        void applyLocalRotations(unsigned int elementID, RRVector &fe);

};

#endif // RSOLVERSTRESS_H
