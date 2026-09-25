#ifndef RSOLVERGENERIC_H
#define RSOLVERGENERIC_H

#include <QString>

#include <vector>

#include <rbl_book.h>
#include <rml_sparse_matrix.h>

#include "rlocalrotation.h"
#include "rscales.h"
#include "rsolvershareddata.h"

template <typename T>
struct RSolverCartesianVector
{
    public:

        T x, y, z;

};

class RSolverGeneric
{

    protected:

        //! Vacuum permittivity.
        static const double e0;
        //! Stefan–Boltzmann constant.
        static const double sigma;

    protected:

        //! Indicator whether the mesh has changed.
        bool meshChanged;

        //! Problem type.
        RProblemType problemType;
        //! Pointer to model object.
        RModel *pModel;
        //! Model file name.
        QString modelFileName;
        //! Convergence file.
        QString convergenceFileName;
        //! Matrix M (modal analysis).
        RSparseMatrix M;
        //! Matrix A.
        RSparseMatrix A;
        //! Vector x.
        RRVector x;
        //! Vector b.
        RRVector b;
        //! Node book.
        RBook nodeBook;
        //! Local rotations.
        std::vector<RLocalRotation> localRotations;
        //! Element temperature vector.
        RRVector elementTemperature;
        //! Scales.
        RScales scales;
        //! Pointer to shared data container.
        RSolverSharedData *pSharedData;
        //! First run indicator.
        bool firstRun;
        //! Parent task transaction currently owns model-result publication.
        bool resultWritesDeferred;
        //! Propagate configured monitoring-file write failures to a parent
        //! transient task transaction instead of logging and continuing.
        bool monitoringOutputFailuresAreFatal;
        //! Task iteration counter.
        uint taskIteration;
        //! Computable elements array.
        RBVector computableElements;
        //! Includable elements array (not included in coefficient matrix).
        RBVector includableElements;
        //! Surface element inward orientation (if normal is pointing inside computable volume element).
        RBVector inwardElements;

    public:

        //! Constructor.
        RSolverGeneric(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Copy constructor (deleted - solvers should not be copied).
        RSolverGeneric(const RSolverGeneric &genericSolver) = delete;

        //! Assignment operator (deleted - solvers should not be copied).
        RSolverGeneric & operator =(const RSolverGeneric &genericSolver) = delete;

        //! Destructor.
        virtual ~RSolverGeneric();

        //! Run solver.
        void run(bool firstExecution, uint taskIteration);

        //! Defer model-file result publication while a parent task transaction
        //! is collecting results from multiple sub-solvers.
        void setResultWritesDeferred(bool deferred);

        //! Make configured monitoring-file write failures fatal for the run.
        void setMonitoringOutputFailuresAreFatal(bool fatal);

        //! Publish the result accumulated while deferred writes were enabled.
        void publishDeferredResults();

        //! Get mesh changed.
        bool getMeshChanged() const;

        //! Set mesh changed.
        void setMeshChanged(bool meshChanged);

        //! Check if solver has converged.
        virtual bool hasConverged() const = 0;

        //! Number of modes the last eigen-value solve actually produced. The
        //! eigen-value solver can return fewer than were requested, so the
        //! result loop must not assume it got all of them.
        virtual uint getNComputedModes() const;

        //! Update old records.
        static void updateOldRecords(const RTimeSolver &rTimeSolver, const QString &modelFileName);

    protected:

        //! Find element sizes.
        RRVector findElementSizes() const;

        //! Find element scale.
        double findElementScale(bool onlyVolumes = false) const;

        //! Find mesh scale.
        double findMeshScale() const;

        //! Update local rotations.
        virtual void updateLocalRotations();

        //! Clear shared data.
        virtual void clearSharedData();

        //! Store shared data.
        virtual void storeSharedData();

        //! Recover shared data.
        virtual void recoverSharedData();

        //! Initialize.
        virtual void initialize() = 0;

        //! Update scales.
        virtual void updateScales() = 0;

        //! Recover previously computed results.
        virtual void recover() = 0;

        //! Prepare solver.
        virtual void prepare() = 0;

        //! Run matrix solver.
        virtual void solve() = 0;

        //! Process solver results.
        virtual void process() = 0;

        //! Store solver results.
        virtual void store() = 0;

        //! Process statistics.
        virtual void statistics() = 0;

        //! Write results.
        void writeResults();

        //! A structural-dynamics step is committed only after result output
        //! succeeds. Derived solvers can preserve their private integrator
        //! state alongside the model/shared-data transaction in run().
        virtual void beginStepTransaction();
        virtual void rollbackStepTransaction();
        virtual void commitStepTransaction();

        //! Validate mutable model settings before the run changes scales,
        //! results or solver state. Most solvers have no additional gate.
        virtual void validateRunConfiguration() const;

        //! Apply displacement if possible.
        void applyDisplacement();

        //! Remove displacement if possible.
        void removeDisplacement();

        //! Generate node book.
        void generateNodeBook(RProblemType problemType);

        //! Rebuild node book from disabled positions.
        static void rebuildNodeBook(RBook &nodeBook, const RBVector &disabledPositions);

        //! Generate material element vector.
        void generateMaterialVecor(RMaterialProperty::Type materialPropertyType, RRVector &materialPropertyValues) const;

        //! Generate variable element vector from applied environment/initial/boundary condition.
        void generateVariableVector(RVariableType variableType,
                                    RRVector &variableValues,
                                    RBVector &setValues,
                                    bool boundaryConditions,
                                    bool initialConditions,
                                    bool environmentConditions,
                                    bool onlyExplicitBcs = false) const;

        //! Generate element vector holding the heat source density prescribed by
        //! the Heat boundary condition.
        //! The condition prescribes the total heat input [W] for the whole entity,
        //! while the assembly integrates the element value over the element measure.
        //! The total is therefore spread over the measure of the entity it is
        //! applied to - its volume, area or length - which turns it into [W/m^3],
        //! [W/m^2] or [W/m] respectively. A point entity carries no measure, so the
        //! total is split evenly between its elements.
        void generateHeatVector(RRVector &heatValues, RBVector &setValues) const;

        //! Return summed measure - volume, area or length - of all computable
        //! elements of the group. Point groups return the number of such elements.
        double findElementGroupMeasure(const RElementGroup &rElementGroup) const;

        //! Recover variable.
        void recoverVariable(RVariableType variableType,
                             RVariableApplyType applyType,
                             uint expectedSize,
                             uint vectorPosition,
                             RRVector &variableValues,
                             double defaultValue) const;

        //! Sunc values with shared data.
        void syncShared(const QString &keyName, RRVector &values);

        //! Find computable elements based on assigned material properties or has boundary condition assigned.
        virtual void findComputableElements(RProblemType problemType);

        //! Find includable elements (not computed in coefficient matrix) - all points are computable.
        void findIncludableElements();

        //! Find inward surface elements (normal is pointing inside computable element).
        void findInwardElements();

        //! Process monitoring points.
        void processMonitoringPoints() const;

        //! Print results statistics.
        void printStats(RVariableType variableType) const;

};

#endif // RSOLVERGENERIC_H
