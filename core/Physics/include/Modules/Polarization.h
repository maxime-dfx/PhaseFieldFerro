#pragma once
//
// Polarization.h
// ------------------
// Module Polarization complet : equation elementaire, 
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

class PolarizationEquation : public IElementEquation {
public:
    // --- TAGS POUR TRACY ---
    static constexpr const char* ModuleName       = "Polarization";
    static constexpr uint32_t    Color            = PROFILE_COLOR_POLARIZATION;
    static constexpr const char* ZoneAssemblySlow = "Polarization::Assembly_Triplets";
    static constexpr const char* ZoneAssemblyFast = "Polarization::Assembly_CSR";
    static constexpr const char* ZoneSolveCompute = "Polarization::Solve_Compute";
    static constexpr const char* ZoneSolveApply   = "Polarization::Solve_Apply";
    static constexpr const char* ZonePostProcess  = "Polarization::PostProcessing";

    void compute_element_matrices(const ElementContext& ctx, const PhysicsState& state, 
                                  Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
                                  std::vector<int>& global_dofs) const override 
    {
        const MaterialModel& material = state.polycrystal
            ? state.polycrystal->get_material(ctx.elem_idx, ctx.material_id, *state.materials_manager)
            : state.materials_manager->get_material(ctx.material_id);

        bool is_impermeable = (state.config->fracture.mode == CrackBCType::IMPERMEABLE);
        const double mu_p = material.get_mu_p();
        const double a0   = material.get_a0();
        const double eta_k = material.get_eta_k();
        const double dt = state.dt;

        const Eigen::VectorXd* Px_prev = state.get_field("Px_prev");
        const Eigen::VectorXd* Py_prev = state.get_field("Py_prev");
        const Eigen::VectorXd* Px = state.get_field("Px");
        const Eigen::VectorXd* Py = state.get_field("Py");
        const Eigen::VectorXd* v = state.get_field("v_prev");
        const Eigen::VectorXd* ux = state.get_field("ux_prev");
        const Eigen::VectorXd* uy = state.get_field("uy_prev");
        const Eigen::VectorXd* phi = state.get_field("phi_prev");

        int n_nodes = ctx.elem.get_num_nodes();
        global_dofs.clear();
        std::array<int, 8> geom_indices;
        std::array<int, 8> p_dofs;

        for (int i=0; i<n_nodes; ++i) { 
            int idx = ctx.mesh.get_node_index(ctx.elem_idx, i);
            geom_indices[i] = idx;
            p_dofs[i] = state.p_dof_map ? state.p_dof_map->dof_for(ctx.elem_idx, i) : idx;
            global_dofs.push_back(2 * p_dofs[i]); 
            global_dofs.push_back(2 * p_dofs[i] + 1); 
        }

        Eigen::RowVectorXd N_buffer = Eigen::RowVectorXd::Zero(8);
        Eigen::MatrixXd grad_N_buffer = Eigen::MatrixXd::Zero(2, 8);

        ElementIntegrator::integrate(ctx.elem, ctx.coords, N_buffer, grad_N_buffer, [&](int num_n, double dV, const GaussPoint2D& gp) {
            (void)gp;
            double v_gp = 0.0;
            double Px_prev_gp = 0.0, Py_prev_gp = 0.0;
            Eigen::Vector2d P_gp = Eigen::Vector2d::Zero();
            Eigen::Matrix2d eps_gp = Eigen::Matrix2d::Zero();
            Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();

            for (int i = 0; i < num_n; ++i) {
                int g_idx = geom_indices[i]; 
                int p_idx = p_dofs[i];
                double Ni = N_buffer(i);

                P_gp(0) += Ni * (Px ? (*Px)[p_idx] : 0.0);
                P_gp(1) += Ni * (Py ? (*Py)[p_idx] : 0.0);
                
                Px_prev_gp += Ni * (Px_prev ? (*Px_prev)[p_idx] : 0.0);
                Py_prev_gp += Ni * (Py_prev ? (*Py_prev)[p_idx] : 0.0);

                v_gp += Ni * (v ? (*v)[g_idx] : 0.0);

                double ux_val = ux ? (*ux)[g_idx] : 0.0;
                double uy_val = uy ? (*uy)[g_idx] : 0.0;
                eps_gp(0,0) += grad_N_buffer(0, i) * ux_val;
                eps_gp(1,1) += grad_N_buffer(1, i) * uy_val;
                eps_gp(0,1) += 0.5 * (grad_N_buffer(1, i) * ux_val + grad_N_buffer(0, i) * uy_val);
                eps_gp(1,0) = eps_gp(0,1);

                double phi_val = phi ? (*phi)[g_idx] : 0.0;
                E_gp(0) -= grad_N_buffer(0, i) * phi_val;
                E_gp(1) -= grad_N_buffer(1, i) * phi_val;
            }

            const double penalite = (v_gp * v_gp) + eta_k;
            const GinzburgLandauTerms GL = material.compute_GL_terms(P_gp, eps_gp, E_gp, penalite, is_impermeable);
            const double inv_dt = mu_p / dt;

            // CORRECTIF : Boucle explicite pour l'assemblage (Remplace Eigen::Map avec Stride dynamique)
            for (int i = 0; i < num_n; ++i) {
                // 1. Assemblage du second membre (Vecteur F_local)
                double Fx_val = dV * (inv_dt * Px_prev_gp + GL.J_11 * P_gp(0) + GL.J_12 * P_gp(1) - GL.force_px) * N_buffer(i);
                double Fy_val = dV * (inv_dt * Py_prev_gp + GL.J_12 * P_gp(0) + GL.J_22 * P_gp(1) - GL.force_py) * N_buffer(i);
                
                F_local(2 * i)     += Fx_val;
                F_local(2 * i + 1) += Fy_val;

                // 2. Assemblage de la matrice de rigidité (Matrice K_local)
                for (int j = 0; j < num_n; ++j) {
                    double mass_term = (inv_dt * dV) * N_buffer(i) * N_buffer(j);
                    double grad_term = (a0 * penalite * dV) * (grad_N_buffer(0, i) * grad_N_buffer(0, j) + grad_N_buffer(1, i) * grad_N_buffer(1, j));
                    double common_ij = mass_term + grad_term;

                    K_local(2 * i,     2 * j)     += common_ij + GL.J_11 * dV * N_buffer(i) * N_buffer(j);
                    K_local(2 * i,     2 * j + 1) += GL.J_12 * dV * N_buffer(i) * N_buffer(j);
                    K_local(2 * i + 1, 2 * j)     += GL.J_12 * dV * N_buffer(i) * N_buffer(j);
                    K_local(2 * i + 1, 2 * j + 1) += common_ij + GL.J_22 * dV * N_buffer(i) * N_buffer(j);
                }
            }
        });
    }
}; // <-- LA FAMEUSE ACCOLADE MANQUANTE !

// ============================= PolarizationDofMap =============================

class PolarizationDofMap : public IDofMap {
private:
    const Mesh& m_mesh;
    int m_num_p_dofs = 0;
    std::vector<int> m_dof_to_node;
public:
    PolarizationDofMap(const Mesh& mesh, const Polycrystal& polycrystal) : m_mesh(mesh) {
        (void)polycrystal; // non utilise par ce fallback continu, voir note ci-dessus
        m_num_p_dofs = mesh.get_num_nodes();
        m_dof_to_node.resize(m_num_p_dofs);
        for (int n = 0; n < m_num_p_dofs; ++n) m_dof_to_node[n] = n;
    }
    int dof_for(int element_id, int local_node) const override {
        return m_mesh.get_node_index(element_id, local_node);
    }
    int node_for(int dof_id) const override {
        return m_dof_to_node[dof_id];
    }
    int num_dofs() const override { return m_num_p_dofs; }
    int num_p_dofs() const { return m_num_p_dofs; }
};

// ============================= PolarizationDofMapper =============================

class PolarizationDofMapper : public IPhysicsDofMapper {
private:
    const Mesh& m_mesh;
    const BoundaryManager& m_bc;
    const Polycrystal& m_polycrystal;
    const MaterialManager& m_materials;
    const Datafile& m_config; 
    int m_gb_material_id;
    PolarizationDofMap m_dof_map;
    FieldState m_px, m_py;
    double m_P0;
    double m_noise_fraction = 0.01;

    bool element_is_ferroelectric(int e) const {
        if (m_gb_material_id != -1 && m_polycrystal.is_polarization_lock_element(e)) {
            return m_materials.get_material(m_gb_material_id).is_ferroelectric();
        }
        const Element& elem = m_mesh.get_elements()[static_cast<std::size_t>(e)];
        return m_materials.get_material(elem.ref_tag).is_ferroelectric();
    }

    std::vector<bool> compute_node_is_matrix() const {
        int n = m_dof_map.num_p_dofs();
        std::vector<bool> is_matrix(n, false);
        const int num_elements = m_mesh.get_num_elements();
        for (int e = 0; e < num_elements; ++e) {
            if (element_is_ferroelectric(e)) continue;
            const Element& elem = m_mesh.get_elements()[static_cast<std::size_t>(e)];
            for (int i = 0; i < elem.get_num_nodes(); ++i) {
                int p_idx = m_dof_map.dof_for(e, i);
                is_matrix[p_idx] = true;
            }
        }
        return is_matrix;
    }

public:
    PolarizationDofMapper(const Mesh& mesh, const BoundaryManager& bc, const Polycrystal& pc,
                          const Datafile& config, const MaterialManager& materials)
        : m_mesh(mesh), m_bc(bc), m_polycrystal(pc), m_materials(materials),
          m_config(config), m_gb_material_id(config.crystal.grain_boundary_material_id),
          m_dof_map(mesh, pc), m_P0(config.material.P0) {
        int n = m_dof_map.num_p_dofs();
        m_px.resize_and_reset(n); m_py.resize_and_reset(n);
    }

    void register_fields(PhysicsState& state) override {
        state.fields["Px"] = &m_px.current; state.fields["Py"] = &m_py.current;
        state.fields["Px_n"] = &m_px.n; state.fields["Py_n"] = &m_py.n;
        state.fields["Px_prev"] = &m_px.prev_iter; state.fields["Py_prev"] = &m_py.prev_iter;
        state.p_dof_map = &m_dof_map;
    }

    int get_system_size() const override { return m_dof_map.num_p_dofs() * 2; }

    std::vector<FlattenedBC> get_bcs() const override {
        std::vector<bool> is_matrix = compute_node_is_matrix();
        std::vector<FlattenedBC> bcs;
        for (int i = 0; i < m_dof_map.num_p_dofs(); ++i) {
            if (!is_matrix[i]) continue;
            bcs.push_back(FlattenedBC{2 * i, 0.0, true});
            bcs.push_back(FlattenedBC{2 * i + 1, 0.0, true});
        }
        return bcs;
    }

    Eigen::VectorXd build_guess_vector() const override {
        Eigen::VectorXd g(get_system_size());
        for (int i=0; i<m_dof_map.num_p_dofs(); ++i) { g(2*i) = m_px.current(i); g(2*i+1) = m_py.current(i); }
        return g;
    }

    void map_solution_to_states(const Eigen::VectorXd& sol) override {
        for (int i=0; i<m_dof_map.num_p_dofs(); ++i) { m_px.current(i) = sol(2*i); m_py.current(i) = sol(2*i+1); }
    }

    void apply_initial_conditions() override {
        int n = m_dof_map.num_p_dofs();
        std::vector<bool> is_matrix = compute_node_is_matrix();
        
        // On lit les valeurs du TOML !
        double px_val = m_config.polarization.val_x_0; 
        double py_val = m_config.polarization.val_y_0;
        
        std::mt19937 rng(m_polycrystal.get_seed());
        std::uniform_real_distribution<double> noise(-m_noise_fraction * m_P0, m_noise_fraction * m_P0);
        
        Eigen::VectorXd px0(n), py0(n);
        for (int i = 0; i < n; ++i) {
            if (is_matrix[i]) {
                px0(i) = 0.0;
                py0(i) = 0.0;
            } else {
                px0(i) = px_val;
                py0(i) = py_val;
                // On ajoute du bruit uniquement si RANDOM est demandé
                if (m_config.polarization.type == InitializationType::RANDOM) {
                    px0(i) += noise(rng);
                    py0(i) += noise(rng);
                }
            }
        }
        m_px.set(px0);
        m_py.set(py0);
    }

    void save_previous_iteration() override { m_px.save_iteration(); m_py.save_iteration(); }
    void save_previous_state() override { m_px.save_state(); m_py.save_state(); }
    void restore_previous_state() override { m_px.restore_state(); m_py.restore_state(); }
    void update_history() override { m_px.update_history(); m_py.update_history(); }
    double calculate_error() const override {
        return std::max((m_px.current - m_px.prev_iter).cwiseAbs().maxCoeff(), (m_py.current - m_py.prev_iter).cwiseAbs().maxCoeff());
    }
};