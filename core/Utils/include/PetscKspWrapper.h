#ifndef UTILS_PETSC_KSP_WRAPPER_H
#define UTILS_PETSC_KSP_WRAPPER_H

#ifdef USE_PETSC

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <petscksp.h>
#include <string>
#include <vector>

class PetscKspSolverWrapper {
private:
    std::vector<PetscInt> m_petsc_i;
    std::vector<PetscInt> m_petsc_j;
    KSP m_ksp = nullptr;
    Mat m_petsc_K = nullptr;
    Vec m_petsc_b = nullptr;
    Vec m_petsc_x = nullptr;
    MatNullSpace m_near_nullspace = nullptr;

    bool m_is_wrapped = false;
    bool m_use_rigid_body_modes = false;
    double* m_last_value_ptr = nullptr;

    int m_last_iterations = 0;
    double m_last_residual = 0.0;
    KSPConvergedReason m_last_reason = KSP_CONVERGED_ITERATING;

    Eigen::MatrixXd m_node_coords_2d;
    std::string m_zone_label;

    void wrap_matrix(Eigen::SparseMatrix<double>& K);
    void attach_rigid_body_nullspace();
    void eigen_to_petsc(const Eigen::VectorXd& src, Vec dst) const;
    void petsc_to_eigen(Vec src, Eigen::VectorXd& dst) const;

public:
    explicit PetscKspSolverWrapper(std::string zone_label = "Solve");
    ~PetscKspSolverWrapper();

    PetscKspSolverWrapper(const PetscKspSolverWrapper&) = delete;
    PetscKspSolverWrapper& operator=(const PetscKspSolverWrapper&) = delete;

    void enable_rigid_body_modes(const Eigen::MatrixXd& node_coords_2d);
    void configure(Eigen::SparseMatrix<double>& K, double tol, int max_iterations, const std::string& prec_type = "AMG");
    Eigen::VectorXd solve(const Eigen::VectorXd& F, const Eigen::VectorXd& guess);
    void destroy();

    bool converged() const { return m_last_reason > 0; }
    int iterations() const { return m_last_iterations; }
    double residual() const { return m_last_residual; }
    KSPConvergedReason reason() const { return m_last_reason; }
};

#endif // USE_PETSC
#endif // UTILS_PETSC_KSP_WRAPPER_H