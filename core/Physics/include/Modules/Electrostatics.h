#pragma once
//
// Electrostatics.h
// ------------------
// Module Electrostatics complet : equation elementaire, 
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

class ElectrostaticsEquation : public IElementEquation {
public:
    // --- TAGS POUR TRACY ---
    static constexpr const char* ModuleName       = "Electrostatics";
    static constexpr uint32_t    Color            = PROFILE_COLOR_ELECTROSTATICS;
    static constexpr const char* ZoneAssemblySlow = "Electrostatics::Assembly_Triplets";
    static constexpr const char* ZoneAssemblyFast = "Electrostatics::Assembly_CSR";
    static constexpr const char* ZoneSolveCompute = "Electrostatics::Solve_Compute";
    static constexpr const char* ZoneSolveApply   = "Electrostatics::Solve_Apply";
    static constexpr const char* ZonePostProcess  = "Electrostatics::PostProcessing";
    
    void compute_element_matrices(const ElementContext& ctx, const PhysicsState& state, 
                                  Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
                                  std::vector<int>& global_dofs) const override 
    {
        const MaterialModel& material = state.polycrystal
            ? state.polycrystal->get_material(ctx.elem_idx, ctx.material_id, *state.materials_manager)
            : state.materials_manager->get_material(ctx.material_id);

        bool is_impermeable = (state.config->fracture.mode == CrackBCType::IMPERMEABLE);
        const double eta_k = material.get_eta_k();
        

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
            global_dofs.push_back(idx); 
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
                P_gp(0) += N_buffer(i) * (Px ? (*Px)[p_idx] : 0.0);
                P_gp(1) += N_buffer(i) * (Py ? (*Py)[p_idx] : 0.0);
                v_gp    += N_buffer(i) * (v ? (*v)[g_idx] : 0.0);
            }
            
            const double eps_eff = material.compute_effective_permittivity(v_gp, eta_k, is_impermeable);
            const Eigen::Vector2d P_eff = material.compute_effective_polarization(P_gp, v_gp, eta_k, is_impermeable);
            
            auto grad_N = grad_N_buffer.leftCols(num_n);
            K_local.topLeftCorner(num_n, num_n).noalias() += (eps_eff * dV) * (grad_N.transpose() * grad_N);
            F_local.head(num_n).noalias() += dV * (grad_N.transpose() * P_eff);
        });
    }
};

// ============================= ElectrostaticsDofMapper =============================

class ElectrostaticsDofMapper : public IPhysicsDofMapper {
private:
    const Mesh& m_mesh;
    const BoundaryManager& m_bc;
    FieldState m_phi;
public:
    ElectrostaticsDofMapper(const Mesh& mesh, const BoundaryManager& bc) : m_mesh(mesh), m_bc(bc) {
        m_phi.resize_and_reset(mesh.get_num_nodes());
    }
    void register_fields(PhysicsState& state) override {
        state.fields["phi"] = &m_phi.current;
        state.fields["phi_prev"] = &m_phi.prev_iter;
    }
    int get_system_size() const override { return m_mesh.get_num_nodes(); }
    std::vector<FlattenedBC> get_bcs() const override {
        std::vector<FlattenedBC> flat;
        const auto& bcs = m_bc.get_bcs("phi");
        for (size_t i = 0; i < bcs.size(); ++i) {
            if (bcs[i].type == BCType::DIRICHLET) flat.push_back({(int)i, bcs[i].value, true});
        }
        return flat;
    }
    Eigen::VectorXd build_guess_vector() const override { return m_phi.current; }
    void map_solution_to_states(const Eigen::VectorXd& sol) override { m_phi.current = sol; }
    void apply_initial_conditions() override { }
    void save_previous_iteration() override { m_phi.save_iteration(); }
    void save_previous_state() override { m_phi.save_state(); }
    void restore_previous_state() override { m_phi.restore_state(); }
    void update_history() override { m_phi.update_history(); }
    double calculate_error() const override { return 0.0; }
};

// ============================= ElectrostaticsPostProcessor =============================

class ElectrostaticsPostProcessor : public IPostProcessor {
public:
    FieldState Ex, Ey, Dx, Dy;

    void register_fields(PhysicsState& state) override {
        state.fields["Ex"] = &Ex.current;
        state.fields["Ey"] = &Ey.current;
        state.fields["Dx"] = &Dx.current;
        state.fields["Dy"] = &Dy.current;
    }

    void compute(const PhysicsState& state, const Mesh& mesh) override {
        const int n_nodes = mesh.get_num_nodes();
        if (Ex.current.size() != n_nodes) {
            Ex.resize_and_reset(n_nodes);
            Ey.resize_and_reset(n_nodes);
            Dx.resize_and_reset(n_nodes);
            Dy.resize_and_reset(n_nodes);
        }

        const Eigen::VectorXd* Px = state.get_field("Px");
        const Eigen::VectorXd* Py = state.get_field("Py");
        const Eigen::VectorXd* phi = state.get_field("phi");
        const Eigen::VectorXd* v = state.get_field("v");
        const IDofMap* p_dof_map = state.p_dof_map;

        if (!phi || !state.materials_manager) {
            Ex.current.setZero();
            Ey.current.setZero();
            Dx.current.setZero();
            Dy.current.setZero();
            return;
        }

        const bool is_impermeable = state.config && (state.config->fracture.mode == CrackBCType::IMPERMEABLE);

        // ETAPE NUMA : Allocation pure puis initialisation en parallèle
        Eigen::VectorXd num_Ex(n_nodes);
        Eigen::VectorXd num_Ey(n_nodes);
        Eigen::VectorXd num_Dx(n_nodes);
        Eigen::VectorXd num_Dy(n_nodes);
        Eigen::VectorXd weight(n_nodes);

        #pragma omp parallel for schedule(static)
        for (int n = 0; n < n_nodes; ++n) {
            num_Ex(n) = 0.0; num_Ey(n) = 0.0; num_Dx(n) = 0.0; num_Dy(n) = 0.0; weight(n) = 0.0;
        }

        const auto& color_groups = mesh.get_color_groups();

        // ETAPE PARALLÉLISATION : Sans atomiques grâce aux couleurs
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
                        Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
                        double v_gp = 0.0;

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const int p_idx = p_dof_map ? p_dof_map->dof_for(e, i) : g_idx;
                            const double phi_val = (*phi)[g_idx];

                            P_gp(0) += N(i) * (Px ? (*Px)[p_idx] : 0.0);
                            P_gp(1) += N(i) * (Py ? (*Py)[p_idx] : 0.0);
                            v_gp    += N(i) * (v  ? (*v)[g_idx]  : 0.0);

                            E_gp(0) -= grad_N(0, i) * phi_val;
                            E_gp(1) -= grad_N(1, i) * phi_val;
                        }

                        const double eps_eff = material.compute_effective_permittivity(v_gp, material.get_eta_k(), is_impermeable);
                        const Eigen::Vector2d P_eff = material.compute_effective_polarization(P_gp, v_gp, material.get_eta_k(), is_impermeable);
                        const Eigen::Vector2d D_gp = eps_eff * E_gp + P_eff;

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const double w = N(i) * dV;
                            num_Ex(g_idx) += E_gp(0) * w;
                            num_Ey(g_idx) += E_gp(1) * w;
                            num_Dx(g_idx) += D_gp(0) * w;
                            num_Dy(g_idx) += D_gp(1) * w;
                            weight(g_idx) += w;
                        }
                    });
                }
            }
        }

        // Parallélisation de la division finale
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < n_nodes; ++n) {
            const double w = weight(n);
            Ex.current(n) = (w > 1e-14) ? num_Ex(n) / w : 0.0;
            Ey.current(n) = (w > 1e-14) ? num_Ey(n) / w : 0.0;
            Dx.current(n) = (w > 1e-14) ? num_Dx(n) / w : 0.0;
            Dy.current(n) = (w > 1e-14) ? num_Dy(n) / w : 0.0;
        }
    }
};
