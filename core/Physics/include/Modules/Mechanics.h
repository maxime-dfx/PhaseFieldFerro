#pragma once
//
// Mechanics.h
// ------------------
// Module Mechanics complet : equation elementaire, 
// gestion des DOFs/BCs et post-traitement E/D.

// --- INCLUSIONS (Uniques pour tout le fichier) ---
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Physics/include/State/FieldState.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Mesh/include/FEM.h"
#include "Mesh/include/IDofMap.h"
#include "Physics/Core/BoundaryManager.h"
#include "IO/include/Datafile.h"
#include "IO/include/ConfigTypes.h"
#include "Utils/include/Profiling.h"
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

class MechanicsEquation : public IElementEquation {
public:
    // --- TAGS POUR TRACY ---
    static constexpr const char* ModuleName       = "Mechanics";
    static constexpr uint32_t    Color            = PROFILE_COLOR_MECHANICS;
    static constexpr const char* ZoneAssemblySlow = "Mechanics::Assembly_Triplets";
    static constexpr const char* ZoneAssemblyFast = "Mechanics::Assembly_CSR";
    static constexpr const char* ZoneSolveCompute = "Mechanics::Solve_Compute";
    static constexpr const char* ZoneSolveApply   = "Mechanics::Solve_Apply";
    static constexpr const char* ZonePostProcess  = "Mechanics::PostProcessing";

    void compute_element_matrices(const ElementContext& ctx, const PhysicsState& state, 
                                  Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
                                  std::vector<int>& global_dofs) const override 
    {
        const MaterialModel& material = state.polycrystal
            ? state.polycrystal->get_material(ctx.elem_idx, ctx.material_id, *state.materials_manager)
            : state.materials_manager->get_material(ctx.material_id);

        const Eigen::VectorXd* Px = state.get_field("Px_prev");
        const Eigen::VectorXd* Py = state.get_field("Py_prev");
        const Eigen::VectorXd* v = state.get_field("v_prev");
        
        int n_nodes = ctx.elem.get_num_nodes();
        global_dofs.clear();
        std::array<int, 8> geom_indices;
        std::array<int, 8> p_dofs;
        
        for (int i=0; i<n_nodes; ++i) { 
            int idx = ctx.mesh.get_node_index(ctx.elem_idx, i);
            geom_indices[i] = idx;
            p_dofs[i] = state.p_dof_map ? state.p_dof_map->dof_for(ctx.elem_idx, i) : idx;
            global_dofs.push_back(2*idx); 
            global_dofs.push_back(2*idx+1); 
        }
        
        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(8);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, 8);
        
        ElementIntegrator::integrate(ctx.elem, ctx.coords, N_buffer, grad_N_buffer, [&](int num_n, double dV, const GaussPoint2D& gp) {
            (void)gp;
            Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
            double v_gp = 0.0;
            
            for (int i = 0; i < num_n; ++i) {
                int g_idx = geom_indices[i];
                int p_idx = p_dofs[i];
                P_gp(0) += N_buffer(i) * (Px ? (*Px)(p_idx) : 0.0);
                P_gp(1) += N_buffer(i) * (Py ? (*Py)(p_idx) : 0.0);
                v_gp    += N_buffer(i) * (v ? (*v)(g_idx) : 0.0);
            }
            
            Eigen::Matrix<double, 3, Eigen::Dynamic, 0, 3, 16> B; B.setZero(3, 2 * num_n);
            for (int i = 0; i < num_n; ++i) {
                B(0, 2*i) = grad_N_buffer(0, i);
                B(1, 2*i+1) = grad_N_buffer(1, i);
                B(2, 2*i) = grad_N_buffer(1, i);
                B(2, 2*i+1) = grad_N_buffer(0, i);
            }
            
            const Eigen::Matrix3d C = material.get_elastic_matrix();
            const Eigen::Vector3d sigma_0 = material.compute_sigma_0(P_gp);
            const double degradation = (v_gp * v_gp) + material.get_eta_k();
            
            int n_local_dofs = 2 * num_n;
            K_local.topLeftCorner(n_local_dofs, n_local_dofs).noalias() += B.transpose() * (C * degradation) * B * dV;
            F_local.head(n_local_dofs).noalias() -= B.transpose() * (sigma_0 * degradation) * dV;
        });
    }
};

// ============================= MechanicsDofMapper =============================

class MechanicsDofMapper : public IPhysicsDofMapper {
private:
    const Mesh& m_mesh;
    const BoundaryManager& m_bc;
    FieldState m_ux, m_uy;
public:
    MechanicsDofMapper(const Mesh& mesh, const BoundaryManager& bc) : m_mesh(mesh), m_bc(bc) {
        int n = m_mesh.get_num_nodes();
        m_ux.resize_and_reset(n); m_uy.resize_and_reset(n);
    }
    void register_fields(PhysicsState& state) override {
        state.fields["ux"] = &m_ux.current; state.fields["uy"] = &m_uy.current;
        state.fields["ux_prev"] = &m_ux.prev_iter; state.fields["uy_prev"] = &m_uy.prev_iter;
    }
    int get_system_size() const override { return m_mesh.get_num_nodes() * 2; }
    std::vector<FlattenedBC> get_bcs() const override {
        std::vector<FlattenedBC> flat;
        const auto& ux_bcs = m_bc.get_bcs("ux");
        const auto& uy_bcs = m_bc.get_bcs("uy");
        for (size_t i = 0; i < ux_bcs.size(); ++i) {
            if (ux_bcs[i].type == BCType::DIRICHLET) flat.push_back({(int)(2*i), ux_bcs[i].value, true});
            if (uy_bcs[i].type == BCType::DIRICHLET) flat.push_back({(int)(2*i+1), uy_bcs[i].value, true});
        }
        return flat;
    }
    Eigen::VectorXd build_guess_vector() const override {
        Eigen::VectorXd g(get_system_size());
        for (int i=0; i<m_mesh.get_num_nodes(); ++i) { g(2*i) = m_ux.current(i); g(2*i+1) = m_uy.current(i); }
        return g;
    }
    void map_solution_to_states(const Eigen::VectorXd& sol) override {
        for (int i=0; i<m_mesh.get_num_nodes(); ++i) { m_ux.current(i) = sol(2*i); m_uy.current(i) = sol(2*i+1); }
    }
    void apply_initial_conditions() override { }
    void save_previous_iteration() override { m_ux.save_iteration(); m_uy.save_iteration(); }
    void save_previous_state() override { m_ux.save_state(); m_uy.save_state(); }
    void restore_previous_state() override { m_ux.restore_state(); m_uy.restore_state(); }
    void update_history() override { m_ux.update_history(); m_uy.update_history(); }
    double calculate_error() const override { return 0.0; }
};

// ============================= StressPostProcessor =============================

class StressPostProcessor : public IPostProcessor {
public:
    FieldState sigma_xx, sigma_yy, sigma_xy, von_mises;

    void register_fields(PhysicsState& state) override {
        state.fields["sigma_xx"] = &sigma_xx.current;
        state.fields["sigma_yy"] = &sigma_yy.current;
        state.fields["sigma_xy"] = &sigma_xy.current;
        state.fields["von_mises"] = &von_mises.current;
    }

    void compute(const PhysicsState& state, const Mesh& mesh) override {
        const int n_nodes = mesh.get_num_nodes();
        if (sigma_xx.current.size() != n_nodes) {
            sigma_xx.resize_and_reset(n_nodes);
            sigma_yy.resize_and_reset(n_nodes);
            sigma_xy.resize_and_reset(n_nodes);
            von_mises.resize_and_reset(n_nodes);
        }

        const Eigen::VectorXd* Px = state.get_field("Px");
        const Eigen::VectorXd* Py = state.get_field("Py");
        const Eigen::VectorXd* ux = state.get_field("ux");
        const Eigen::VectorXd* uy = state.get_field("uy");
        const Eigen::VectorXd* v  = state.get_field("v");
        const IDofMap* p_dof_map = state.p_dof_map;

        if (!ux || !uy || !state.materials_manager) {
            sigma_xx.current.setZero();
            sigma_yy.current.setZero();
            sigma_xy.current.setZero();
            von_mises.current.setZero();
            return;
        }

        // ETAPE 5 (NUMA) : Allocation sans initialisation, puis mise à zéro parallèle
        Eigen::VectorXd num_xx(n_nodes);
        Eigen::VectorXd num_yy(n_nodes);
        Eigen::VectorXd num_xy(n_nodes);
        Eigen::VectorXd weight(n_nodes);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_nodes; ++i) {
            num_xx(i) = 0.0; num_yy(i) = 0.0; num_xy(i) = 0.0; weight(i) = 0.0;
        }

        const auto& color_groups = mesh.get_color_groups();

        // ETAPE 4 : Parallélisation SANS atomiques grâce aux couleurs
        #pragma omp parallel
        {
            Eigen::RowVectorXd N(8);
            Eigen::MatrixXd grad_N(2, 8);

            for (size_t c = 0; c < color_groups.size(); ++c) {
                const auto& group = color_groups[c];

                #pragma omp for schedule(guided)
                for (size_t idx = 0; idx < group.size(); ++idx) {
                    int e = group[idx];
                    
                    const Element& elem = mesh.get_elements()[e];
                    auto coords = mesh.get_element_coords(e);
                    const int mat_id = elem.ref_tag;
                    const MaterialModel& material = state.polycrystal
                        ? state.polycrystal->get_material(e, mat_id, *state.materials_manager)
                        : state.materials_manager->get_material(mat_id); 

                    ElementIntegrator::integrate(elem, coords, N, grad_N, [&](int n_local, double dV, const GaussPoint2D& gp) {
                        (void)gp;
                        Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
                        double v_gp = 0.0;
                        Eigen::Matrix2d eps_gp = Eigen::Matrix2d::Zero();

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const int p_idx = p_dof_map ? p_dof_map->dof_for(e, i) : g_idx;
                            const double ux_val = (*ux)[g_idx];
                            const double uy_val = (*uy)[g_idx];

                            P_gp(0) += N(i) * (Px ? (*Px)[p_idx] : 0.0);
                            P_gp(1) += N(i) * (Py ? (*Py)[p_idx] : 0.0);
                            v_gp    += N(i) * (v  ? (*v)[g_idx]  : 0.0);

                            eps_gp(0,0) += grad_N(0,i) * ux_val;
                            eps_gp(1,1) += grad_N(1,i) * uy_val;
                            eps_gp(0,1) += 0.5 * (grad_N(1,i) * ux_val + grad_N(0,i) * uy_val);
                            eps_gp(1,0) = eps_gp(0,1);
                        }

                        const Eigen::Vector3d eps_voigt(eps_gp(0,0), eps_gp(1,1), 2.0 * eps_gp(0,1));
                        const Eigen::Matrix3d C = material.get_elastic_matrix();
                        const Eigen::Vector3d sigma_0 = material.compute_sigma_0(P_gp);
                        const double degradation = (v_gp * v_gp) + material.get_eta_k();
                        const Eigen::Vector3d sigma = degradation * (C * eps_voigt - sigma_0);

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const double w = N(i) * dV;
                            // Ecriture directe mathématiquement sûre !
                            num_xx(g_idx) += sigma(0) * w;
                            num_yy(g_idx) += sigma(1) * w;
                            num_xy(g_idx) += sigma(2) * w;
                            weight(g_idx) += w;
                        }
                    });
                }
            }
        }

        // Parallélisation de la division nodale et von Mises
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < n_nodes; ++n) {
            const double w = weight(n);
            sigma_xx.current(n) = (w > 1e-14) ? num_xx(n) / w : 0.0;
            sigma_yy.current(n) = (w > 1e-14) ? num_yy(n) / w : 0.0;
            sigma_xy.current(n) = (w > 1e-14) ? num_xy(n) / w : 0.0;

            const double sxx = sigma_xx.current(n);
            const double syy = sigma_yy.current(n);
            const double sxy = sigma_xy.current(n);
            von_mises.current(n) = std::sqrt(std::max(0.0,
                sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy));
        }
    }
};
