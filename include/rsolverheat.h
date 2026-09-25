#ifndef RSOLVERHEAT_H
#define RSOLVERHEAT_H

#include <rbl_r3vector.h>

#include <utility>

#include "rsolvergeneric.h"

class RSolverHeat : public RSolverGeneric
{

    protected:

        //! Element heat capacity vector.
        RRVector elementCapacity;
        //! Element thermal conduction vector.
        RRVector elementConduction;
        //! Element density vector.
        RRVector elementDensity;
        //! Node temperature.
        RRVector nodeTemperature;
        //! Element heat vector.
        RRVector elementHeat;
        //! Element heat rate per unit area - surface elements only.
        RRVector elementHeatRateArea;
        //! Element heat rate per unit volume - volume elements only.
        RRVector elementHeatRateVolume;
        //! Element radiation heat vector.
        RRVector elementRadiativeHeat;
        //! Element joule heat.
        RRVector elementJouleHeat;
        //! Element heat flux vector.
        std::vector<RR3Vector> elementHeatFlux;
        //! Volume element on the fluid side of each surface element, or
        //! RConstants::eod where the surface does not border a fluid domain.
        RUVector fluidElements;
        //! Fluid node temperature recovered from the fluid heat solver.
        //! Empty when no fluid heat solve has run.
        RRVector fluidNodeTemperature;
        //! Fluid node velocity magnitude recovered from the fluid heat solver.
        //! Empty when no fluid heat solve has run.
        RRVector fluidNodeVelocity;
        //! Element heat transfer coefficient.
        RRVector elementHeatTransferCoefficient;
        //! True only after the most recent run has passed a fresh linear residual check.
        mutable bool solveCertified;

    public:

        //! Constructor.
        explicit RSolverHeat(RModel *pModel, const QString &modelFileName, const QString &convergenceFileName, RSolverSharedData &sharedData);

        //! Destructor.
        ~RSolverHeat();

        //! Check if solver has converged.
        bool hasConverged() const override;

        //! Additive W inputs on uniquely identified POINT elements for one
        //! trial heat solve. The caller owns the enclosing macrostep rollback.
        void setControlledPointHeat(const std::vector<std::pair<uint,double>> &powers);

    protected:

        std::vector<std::pair<uint,double>> controlledPointHeat;

        void findComputableElements(RProblemType problemType) override;

        //! Recover previously computed results from the shared data container.
        void recoverSharedData() override;

        //! Return the temperature the fluid heat solver computed on the fluid side
        //! of the given surface element. False when no such result covers it - the
        //! surface borders no fluid, or no fluid heat solve has run yet.
        bool findFluidTemperature(unsigned int elementId, double &fluidTemperature) const;

        //! Return the mean velocity the fluid solver computed on the fluid side of
        //! the given surface element. False when no such result covers it, or when
        //! the fluid is at rest and the correlation would collapse.
        bool findFluidVelocity(unsigned int elementId, double &fluidVelocity) const;

        //! Find the fluid volume element attached to every surface element.
        //! The Forced convection condition does not prescribe a fluid temperature -
        //! it takes the one the fluid heat solver computed on the other side of the
        //! wall, so the wall has to know which element holds it.
        void findFluidElements();

        //! Find temperature scale.
        double findTemperatureScale() const;

        //! Initialize solver.
        void initialize() override;

        //! Clear the previous solve status and validate the linear tolerance at run entry.
        void validateRunConfiguration() const override;

        //! Update scales.
        void updateScales() override;

        //! Recover previously computed results.
        void recover() override;

        //! Prepare solver.
        void prepare() override;

        //! Run matrix solver.
        void solve() override;

        //! Process solver results.
        void process() override;

        //! Store solver results.
        void store() override;

        //! Process statistics.
        void statistics() override;

        //! Assembly matrix
        void assemblyMatrix(unsigned int elementID, const RRMatrix &Me, const RRMatrix &Ke, const RRVector &fe, RSparseMatrix &Ap, RRVector &bp);

        //! Get simple convection BC values.
        bool getSimpleConvection(const RElementGroup &elementGroup, double &htc, double &htt);

        //! Throw when a value configured on a convection condition would leave the
        //! correlation with nothing to work with - the dimensionless groups divide by
        //! the viscosity, the thermal conductivity and the hydraulic diameter, and a
        //! zero density, heat capacity or mean velocity collapses them just as surely.
        void checkConvectionInput(double value,
                                  RVariableType variableType,
                                  RBoundaryConditionType boundaryConditionType,
                                  const RElementGroup &elementGroup) const;

        //! Report a Forced convection condition the solver cannot make use of -
        //! one on a surface bordering no fluid, or one whose properties give no
        //! heat transfer at all.
        void reportForcedConvection(const RElementGroup &elementGroup);

        //! Get forced convection BC values. The fluid temperature comes from the
        //! fluid heat solver results, and only where those are missing from the
        //! value configured on the condition itself.
        bool getForcedConvection(const RElementGroup &elementGroup, unsigned int elementId, double &htc, double &htt);

        //! Get natural convection BC values.
        bool getNaturalConvection(const RElementGroup &elementGroup, unsigned int elementId, double &htc, double &htt);

};

#endif // RSOLVERHEAT_H
