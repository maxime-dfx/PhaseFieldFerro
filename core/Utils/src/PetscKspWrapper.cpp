#include "Utils/include/PetscKspWrapper.h"
#include "Utils/include/Logger.h"
#include "Utils/include/Profiling.h"
#include <cstring>

#ifdef USE_PETSC

PetscKspSolverWrapper::PetscKspSolverWrapper(std::string zone_label)
    : m_zone_label(std::move(zone_label)) {}

PetscKspSolverWrapper::~PetscKspSolverWrapper() {
    destroy();
}

void PetscKspSolverWrapper::destroy() {
    if (m_near_nullspace) { MatNullSpaceDestroy(&m_near_nullspace); m_near_nullspace = nullptr; }
    if (m_petsc_x) { VecDestroy(&m_petsc_x); m_petsc_x = nullptr; }
    if (m_petsc_b) { VecDestroy(&m_petsc_b); m_petsc_b = nullptr; }
    if (m_petsc_K) { MatDestroy(&m_petsc_K); m_petsc_K = nullptr; }
    if (m_ksp)     { KSPDestroy(&m_ksp);     m_ksp = nullptr; }
    m_is_wrapped = false;
    m_last_value_ptr = nullptr;
}

void PetscKspSolverWrapper::enable_rigid_body_modes(const Eigen::MatrixXd& node_coords_2d) {
    m_use_rigid_body_modes = true;
    m_node_coords_2d = node_coords_2d;
}

void PetscKspSolverWrapper::eigen_to_petsc(const Eigen::VectorXd& src, Vec dst) const {
    PetscScalar* arr = nullptr;
    VecGetArray(dst, &arr);
    std::memcpy(arr, src.data(), static_cast<size_t>(src.size()) * sizeof(PetscScalar));
    VecRestoreArray(dst, &arr);
}

void PetscKspSolverWrapper::petsc_to_eigen(Vec src, Eigen::VectorXd& dst) const {
    const PetscScalar* arr = nullptr;
    VecGetArrayRead(src, &arr);
    std::memcpy(dst.data(), arr, static_cast<size_t>(dst.size()) * sizeof(PetscScalar));
    VecRestoreArrayRead(src, &arr);
}

void PetscKspSolverWrapper::attach_rigid_body_nullspace() {
    if (!m_use_rigid_body_modes || m_node_coords_2d.rows() == 0) return;

    const PetscInt n_nodes = static_cast<PetscInt>(m_node_coords_2d.rows());
    const PetscInt block_size = 2;
    Vec coords_vec;
    VecCreateSeq(PETSC_COMM_SELF, n_nodes * block_size, &coords_vec);
    PetscScalar* arr = nullptr;
    VecGetArray(coords_vec, &arr);
    for (PetscInt i = 0; i < n_nodes; ++i) {
        arr[2 * i]     = m_node_coords_2d(i, 0);
        arr[2 * i + 1] = m_node_coords_2d(i, 1);
    }
    VecRestoreArray(coords_vec, &arr);

    if (m_near_nullspace) { MatNullSpaceDestroy(&m_near_nullspace); m_near_nullspace = nullptr; }
    MatNullSpaceCreateRigidBody(coords_vec, &m_near_nullspace);
    MatSetNearNullSpace(m_petsc_K, m_near_nullspace);
    VecDestroy(&coords_vec);
}

void PetscKspSolverWrapper::wrap_matrix(Eigen::SparseMatrix<double>& K) {
    PROFILE_ZONE_NC("PETSc_Mat_Wrap", PROFILE_COLOR_SOLVE);
    K.makeCompressed();
    const PetscInt n = static_cast<PetscInt>(K.rows());
    
    const int* eigen_i = K.outerIndexPtr();
    const int* eigen_j = K.innerIndexPtr();
    PetscScalar* v_ptr = reinterpret_cast<PetscScalar*>(K.valuePtr());

    const bool pattern_changed = !m_is_wrapped || v_ptr != m_last_value_ptr;

    if (pattern_changed) {
        if (m_is_wrapped) {
            if (m_petsc_x) { VecDestroy(&m_petsc_x); m_petsc_x = nullptr; }
            if (m_petsc_b) { VecDestroy(&m_petsc_b); m_petsc_b = nullptr; }
            if (m_petsc_K) { MatDestroy(&m_petsc_K); m_petsc_K = nullptr; }
        }

        size_t nnz = K.nonZeros();
        size_t outer_size = K.rows() + 1;

        m_petsc_i.resize(outer_size);
        for (size_t i = 0; i < outer_size; ++i) {
            m_petsc_i[i] = static_cast<PetscInt>(eigen_i[i]);
        }

        m_petsc_j.resize(nnz);
        for (size_t i = 0; i < nnz; ++i) {
            m_petsc_j[i] = static_cast<PetscInt>(eigen_j[i]);
        }

        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, n, n, m_petsc_i.data(), m_petsc_j.data(), v_ptr, &m_petsc_K);
        
        MatSetOption(m_petsc_K, MAT_SYMMETRIC, PETSC_TRUE);
        
        VecCreateSeq(PETSC_COMM_SELF, n, &m_petsc_b);
        VecCreateSeq(PETSC_COMM_SELF, n, &m_petsc_x);

        m_is_wrapped = true;
        m_last_value_ptr = v_ptr;
        attach_rigid_body_nullspace();
    }
    MatAssemblyBegin(m_petsc_K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_petsc_K, MAT_FINAL_ASSEMBLY);
}

void PetscKspSolverWrapper::configure(Eigen::SparseMatrix<double>& K, double tol, int max_iterations, const std::string& prec_type) {
    PROFILE_ZONE_NC("Solve_Compute", PROFILE_COLOR_SOLVE);
    wrap_matrix(K);
    
    if (!m_ksp) {
        KSPCreate(PETSC_COMM_SELF, &m_ksp);
        
        PC pc;
        KSPGetPC(m_ksp, &pc);
        
        // Routage dynamique basé sur config.toml
        if (prec_type == "AMG" || prec_type == "BOOMERAMG") {
            KSPSetType(m_ksp, KSPCG);
            PCSetType(pc, PCHYPRE);
            PCHYPRESetType(pc, "boomeramg");
        } else if (prec_type == "ILU") {
            KSPSetType(m_ksp, KSPGMRES); // ILU se marie très bien avec GMRES
            PCSetType(pc, PCILU);
        } else if (prec_type == "BJACOBI") {
            KSPSetType(m_ksp, KSPCG);
            PCSetType(pc, PCBJACOBI);
        } else {
            KSPSetType(m_ksp, KSPCG);
            PCSetType(pc, PCJACOBI); // Équivalent de Diagonal
        }

        KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE);
        KSPSetFromOptions(m_ksp);
    }
    
    KSPSetTolerances(m_ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iterations);
    KSPSetOperators(m_ksp, m_petsc_K, m_petsc_K);
}

Eigen::VectorXd PetscKspSolverWrapper::solve(const Eigen::VectorXd& F, const Eigen::VectorXd& guess) {
    PROFILE_ZONE_NC("Solve_Apply", PROFILE_COLOR_SOLVE);
    eigen_to_petsc(F, m_petsc_b);
    eigen_to_petsc(guess, m_petsc_x);

    KSPSolve(m_ksp, m_petsc_b, m_petsc_x);

    Eigen::VectorXd result(guess.size());
    petsc_to_eigen(m_petsc_x, result);

    KSPGetIterationNumber(m_ksp, &m_last_iterations);
    KSPGetResidualNorm(m_ksp, &m_last_residual);
    KSPGetConvergedReason(m_ksp, &m_last_reason);

    Logger::debug("[", m_zone_label, "] n_dof=", F.size(),
         " KSP_iters=", m_last_iterations,
         " residual=", m_last_residual,
         " reason=", static_cast<int>(m_last_reason));
    return result;
}

#endif // USE_PETSC