#ifndef RMODAL_DAMPING_H
#define RMODAL_DAMPING_H

#include <vector>

#include <rbl_rmatrix.h>
#include <rml_matrix_solver_conf.h>
#include <rml_sparse_matrix.h>

//! A verified, small-system modal damping operator in reduced solver coordinates.
struct RModalDampingResult
{
    RSparseMatrix damping;
    RRVector eigenvalues;
    //! Rows are mass-normalized eigenvectors; only the configured modes are kept.
    RRMatrix modes;
};

class RModalDamping
{
public:
    static constexpr uint MaxEquations = 512;
    static constexpr uint MaxModes = 16;

    struct Options
    {
        uint eigenIterations = 500;
        double eigenConvergence = 1.0e-14;
    };

    //! C = sum_i 2*zeta_i*sqrt(lambda_i)*(M*phi_i)*(M*phi_i)^T.
    //! The requested ratios correspond to the lowest modes; all other modes
    //! remain undamped. A guard eigenpair is required at a partial-mode boundary.
    static RModalDampingResult synthesize(const RSparseMatrix &mass,
                                          const RSparseMatrix &stiffness,
                                          const std::vector<double> &ratios,
                                          const RMatrixSolverConf &linearSolver);
    static RModalDampingResult synthesize(const RSparseMatrix &mass,
                                          const RSparseMatrix &stiffness,
                                          const std::vector<double> &ratios,
                                          const RMatrixSolverConf &linearSolver,
                                          Options options);
};

#endif // RMODAL_DAMPING_H
