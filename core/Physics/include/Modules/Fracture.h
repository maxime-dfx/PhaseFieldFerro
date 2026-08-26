#pragma once
//
// Fracture.h
// ------------------
// Module Fracture complet : equation elementaire, 
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

class FractureEquation : public IElementEquation {
public:
    // --- TAGS POUR TRACY ---
    static constexpr const char* ModuleName       = "Fracture";
    static constexpr uint32_t    Color            = PROFILE_COLOR_FRACTURE;
    static constexpr const char* ZoneAssemblySlow = "Fracture::Assembly_Triplets";
    static constexpr const char* ZoneAssemblyFast = "Fracture::Assembly_CSR";
    static constexpr const char* ZoneSolveCompute = "Fracture::Solve_Compute";
    static constexpr const char* ZoneSolveApply   = "Fracture::Solve_Apply";
    static constexpr const char* ZonePostProcess  = "Fracture::PostProcessing";
    
    void compute_element_matrices(const ElementContext& ctx, const PhysicsState& state, 
                                  Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
                                  std::vector<int>& global_dofs) const override 
    {
        const MaterialModel& material = state.polycrystal
            ? state.polycrystal->get_material(ctx.elem_idx, ctx.material_id, *state.materials_manager)
            : state.materials_manager->get_material(ctx.material_id);
        bool is_impermeable = (state.config->fracture.mode == CrackBCType::IMPERMEABLE);
    
        double mu_v = material.get_mu_v();
        double kappa = material.get_kappa();
        const double dt = state.dt;
        double Gc_elem = material.get_Gc();
        if (state.polycrystal && state.polycrystal->is_fracture_band_element(ctx.elem_idx)) {
            int gb_id = state.config->crystal.grain_boundary_material_id;
            if (gb_id != -1) {
                const auto& gb_mat = state.materials_manager->get_material(gb_id);
                Gc_elem = gb_mat.get_Gc();
                mu_v = gb_mat.get_mu_v();   
                kappa = gb_mat.get_kappa(); 
            }
        }
        const Eigen::VectorXd* Px_prev = state.get_field("Px_prev");
        const Eigen::VectorXd* Py_prev = state.get_field("Py_prev");
        const Eigen::VectorXd* ux_prev = state.get_field("ux_prev");
        const Eigen::VectorXd* uy_prev = state.get_field("uy_prev");
        const Eigen::VectorXd* phi_prev = state.get_field("phi_prev");
        const Eigen::VectorXd* v_prev = state.get_field("v_prev");
        
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
            Eigen::Matrix2d grad_P_gp = Eigen::Matrix2d::Zero();
            Eigen::Matrix2d eps_gp = Eigen::Matrix2d::Zero();
            Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
            double v_n_gp = 0.0;
            
            for (int i = 0; i < num_n; ++i) {
                int g_idx = geom_indices[i];
                int p_idx = p_dofs[i];
                
                double px_val = Px_prev ? (*Px_prev)[p_idx] : 0.0;
                double py_val = Py_prev ? (*Py_prev)[p_idx] : 0.0;
                
                P_gp(0) += N_buffer(i) * px_val; 
                P_gp(1) += N_buffer(i) * py_val;
                
                grad_P_gp(0, 0) += grad_N_buffer(0,i) * px_val; grad_P_gp(0, 1) += grad_N_buffer(1,i) * px_val;
                grad_P_gp(1, 0) += grad_N_buffer(0,i) * py_val; grad_P_gp(1, 1) += grad_N_buffer(1,i) * py_val;
                
                v_n_gp += N_buffer(i) * (v_prev ? (*v_prev)[g_idx] : 0.0);
                
                double ux_val = ux_prev ? (*ux_prev)[g_idx] : 0.0;
                double uy_val = uy_prev ? (*uy_prev)[g_idx] : 0.0;
                eps_gp(0,0) += grad_N_buffer(0, i) * ux_val;
                eps_gp(1,1) += grad_N_buffer(1, i) * uy_val;
                eps_gp(0,1) += 0.5 * (grad_N_buffer(1, i) * ux_val + grad_N_buffer(0, i) * uy_val);
                eps_gp(1,0) = eps_gp(0,1);
                
                double phi_val = phi_prev ? (*phi_prev)[g_idx] : 0.0;
                E_gp(0) -= grad_N_buffer(0, i) * phi_val;
                E_gp(1) -= grad_N_buffer(1, i) * phi_val;
            }
            
            if (!is_impermeable) E_gp.setZero();
            const double H_drive = material.compute_H_drive(grad_P_gp, P_gp, eps_gp, E_gp, is_impermeable);
            
            double mass_coeff = (mu_v / dt) + (Gc_elem / (2.0 * kappa)) + 2.0 * H_drive;
            const double mass_coeff_floor = 1e-6 * (mu_v / dt);
            if (mass_coeff < mass_coeff_floor) {
                Logger::debug("[FractureEquation] mass_coeff_floor active : mass_coeff=",
                               mass_coeff, " < floor=", mass_coeff_floor,
                               " (H_drive=", H_drive, ", elem=", ctx.elem_idx, ")");
            }
            mass_coeff = std::max(mass_coeff, mass_coeff_floor);
            const double diff_coeff = 2.0 * Gc_elem * kappa;
            const double rhs_coeff = (mu_v / dt) * v_n_gp + (Gc_elem / (2.0 * kappa));
            
            auto N = N_buffer.head(num_n);
            auto grad_N = grad_N_buffer.leftCols(num_n);
            K_local.topLeftCorner(num_n, num_n).noalias() += (mass_coeff * dV) * (N.transpose() * N);
            K_local.topLeftCorner(num_n, num_n).noalias() += (diff_coeff * dV) * (grad_N.transpose() * grad_N);
            F_local.head(num_n).noalias() += (rhs_coeff * dV) * N.transpose();
        });
    }
};

// ============================= FractureDofMapper =============================

class FractureDofMapper : public IPhysicsDofMapper {
private:
    const Mesh& m_mesh;
    const BoundaryManager& m_bc;
    const Datafile& m_config; 
    FieldState m_v;
    double m_alpha;

public:
    FractureDofMapper(const Mesh& mesh, const BoundaryManager& bc, const Datafile& config)
        : m_mesh(mesh), m_bc(bc), m_config(config), m_alpha(config.fracture.alpha_irreversibility) {
        m_v.resize_and_reset(mesh.get_num_nodes());
    }

    void register_fields(PhysicsState& state) override {
        state.fields["v"] = &m_v.current;
        state.fields["v_n"] = &m_v.n;
        state.fields["v_prev"] = &m_v.prev_iter;
    }

    int get_system_size() const override { return m_mesh.get_num_nodes(); }

    std::vector<FlattenedBC> get_bcs() const override {
        std::vector<FlattenedBC> flat;
        const int n = m_mesh.get_num_nodes();
        flat.reserve(static_cast<size_t>(n) / 8);
        for (int i = 0; i < n; ++i) {
            if (m_v.n(i) <= m_alpha) {
                flat.push_back({i, 0.0, true});
            }
        }
        return flat;
    }

    Eigen::VectorXd build_guess_vector() const override { return m_v.current; }

    void map_solution_to_states(const Eigen::VectorXd& sol) override {
        for (int i=0; i<sol.size(); ++i) m_v.current(i) = std::max(0.0, std::min(sol(i), m_v.n(i)));
    }

    void apply_initial_conditions() override { 
        m_v.current.setOnes(); 
        m_v.n.setOnes(); 

        // Application de la pré-fissure spatiale demandée par le TOML !
        if (m_config.fracture.enable_precrack) {
            const auto& pc = m_config.fracture.precrack;
            for (int i = 0; i < m_mesh.get_num_nodes(); ++i) {
                auto coords = m_mesh.get_node_coords(i);
                double x = coords[0], y = coords[1];
                
                double d = 1e9;
                if (pc.shape == PrecrackShape::SEGMENT) {
                    if (x < pc.x0) d = std::sqrt((x-pc.x0)*(x-pc.x0) + (y-pc.y0)*(y-pc.y0));
                    else if (x > pc.x0 + pc.length) d = std::sqrt((x-pc.x0-pc.length)*(x-pc.x0-pc.length) + (y-pc.y0)*(y-pc.y0));
                    else d = std::abs(y - pc.y0);
                }
                
                if (d < 1e8) {
                    if (pc.smooth) {
                        double v_val = 1.0 - std::exp(-d / pc.smoothing_length);
                        m_v.current(i) = std::min(m_v.current(i), v_val);
                        m_v.n(i) = std::min(m_v.n(i), v_val);
                    } else if (d <= pc.half_width) {
                        m_v.current(i) = 0.0;
                        m_v.n(i) = 0.0;
                    }
                }
            }
        }
    }

    void save_previous_iteration() override { m_v.save_iteration(); }
    void save_previous_state() override { m_v.save_state(); }
    void restore_previous_state() override { m_v.restore_state(); }
    void update_history() override { m_v.update_history(); }
    double calculate_error() const override { return (m_v.current - m_v.prev_iter).cwiseAbs().maxCoeff(); }
};
