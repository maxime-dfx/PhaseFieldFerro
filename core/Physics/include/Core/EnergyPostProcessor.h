#pragma once
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Physics/include/State/FieldState.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Microstructure/Polycrystal.h"
#include "Mesh/include/FEM.h"
#include "Mesh/include/IDofMap.h"
#include "IO/include/Datafile.h"
#include "IO/include/ConfigTypes.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Diagnostic transverse (lit Px/Py, ux/uy, phi, v depuis le blackboard,
// exactement comme SystemMetricsCalculator dont il reprend les formules de
// densite d'energie), donc rattache a Core plutot qu'a un module particulier
// : ce n'est ni un worker de Polarization/Mechanics/Electrostatics/Fracture,
// ni une equation resolue, juste une projection L2 aux noeuds (meme
// convention masse lumpee que StressPostProcessor) des densites d'energie
// de SystemMetricsCalculator::compute, pour visualisation VTK spatiale en
// complement du total scalaire deja ecrit dans energies.csv.
class EnergyPostProcessor : public IPostProcessor {
public:
    FieldState energy_gradient;   // U_energy(grad_P), degradee par la fracture
    FieldState energy_elastic;    // W_energy(P, strain), degradee par la fracture
    FieldState energy_landau;     // chi_energy(P) (terme de Landau / anisotropie)
    FieldState energy_electric;   // -P.E - 0.5*eps0*|E|^2 (degradee si IMPERMEABLE)
    FieldState energy_surface;    // densite de surface AT-type (Gc/kappa), fracture
    FieldState energy_bulk;       // somme gradient+elastic+landau+electric

    void register_fields(PhysicsState& state) override {
        state.fields["energy_gradient"] = &energy_gradient.current;
        state.fields["energy_elastic"]  = &energy_elastic.current;
        state.fields["energy_landau"]   = &energy_landau.current;
        state.fields["energy_electric"] = &energy_electric.current;
        state.fields["energy_surface"]  = &energy_surface.current;
        state.fields["energy_bulk"]     = &energy_bulk.current;
    }

    void compute(const PhysicsState& state, const Mesh& mesh) override {
        const int n_nodes = mesh.get_num_nodes();
        if (energy_gradient.current.size() != n_nodes) {
            energy_gradient.resize_and_reset(n_nodes);
            energy_elastic.resize_and_reset(n_nodes);
            energy_landau.resize_and_reset(n_nodes);
            energy_electric.resize_and_reset(n_nodes);
            energy_surface.resize_and_reset(n_nodes);
            energy_bulk.resize_and_reset(n_nodes);
        }

        if (!state.materials_manager || !state.config) {
            energy_gradient.current.setZero();
            energy_elastic.current.setZero();
            energy_landau.current.setZero();
            energy_electric.current.setZero();
            energy_surface.current.setZero();
            energy_bulk.current.setZero();
            return;
        }

        const Eigen::VectorXd* Px = state.get_field("Px");
        const Eigen::VectorXd* Py = state.get_field("Py");
        const Eigen::VectorXd* ux = state.get_field("ux");
        const Eigen::VectorXd* uy = state.get_field("uy");
        const Eigen::VectorXd* phi = state.get_field("phi");
        const Eigen::VectorXd* v = state.get_field("v");
        const IDofMap* p_dof_map = state.p_dof_map;

        const bool is_impermeable = (state.config->fracture.mode == CrackBCType::IMPERMEABLE);
        // Gc/kappa passaient par la constante globale Datafile::material ;
        // ils sont maintenant lus par élément plus bas (via `material`,
        // résolu depuis MaterialManager) pour refléter la ténacité propre
        // à chaque phase du composite.

        // ETAPE NUMA : Allocation pure puis initialisation en parallèle
        Eigen::VectorXd num_U(n_nodes);
        Eigen::VectorXd num_W(n_nodes);
        Eigen::VectorXd num_chi(n_nodes);
        Eigen::VectorXd num_elec(n_nodes);
        Eigen::VectorXd num_surf(n_nodes);
        Eigen::VectorXd weight(n_nodes);

        #pragma omp parallel for schedule(static)
        for (int n = 0; n < n_nodes; ++n) {
            num_U(n) = 0.0; num_W(n) = 0.0; num_chi(n) = 0.0; 
            num_elec(n) = 0.0; num_surf(n) = 0.0; weight(n) = 0.0;
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
                        Eigen::Matrix2d grad_P_gp = Eigen::Matrix2d::Zero();
                        Eigen::Matrix2d strain_gp = Eigen::Matrix2d::Zero();
                        Eigen::Vector2d E_gp = Eigen::Vector2d::Zero();
                        Eigen::Vector2d grad_v_gp = Eigen::Vector2d::Zero();
                        double v_gp = 0.0;

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const int p_idx = p_dof_map ? p_dof_map->dof_for(e, i) : g_idx;

                            const double px_val = Px ? (*Px)[p_idx] : 0.0;
                            const double py_val = Py ? (*Py)[p_idx] : 0.0;
                            const double ux_val = ux ? (*ux)[g_idx] : 0.0;
                            const double uy_val = uy ? (*uy)[g_idx] : 0.0;
                            const double phi_val = phi ? (*phi)[g_idx] : 0.0;
                            const double v_val_node = v ? (*v)[g_idx] : 0.0;

                            P_gp(0) += N(i) * px_val;
                            P_gp(1) += N(i) * py_val;
                            grad_P_gp(0,0) += px_val * grad_N(0,i); grad_P_gp(0,1) += px_val * grad_N(1,i);
                            grad_P_gp(1,0) += py_val * grad_N(0,i); grad_P_gp(1,1) += py_val * grad_N(1,i);

                            strain_gp(0,0) += ux_val * grad_N(0,i);
                            strain_gp(1,1) += uy_val * grad_N(1,i);
                            strain_gp(0,1) += 0.5 * (ux_val * grad_N(1,i) + uy_val * grad_N(0,i));
                            strain_gp(1,0) = strain_gp(0,1);

                            E_gp(0) -= grad_N(0,i) * phi_val;
                            E_gp(1) -= grad_N(1,i) * phi_val;

                            v_gp += N(i) * v_val_node;
                            grad_v_gp(0) += v_val_node * grad_N(0,i);
                            grad_v_gp(1) += v_val_node * grad_N(1,i);
                        }

                        const double deg = (v_gp * v_gp) + material.get_eta_k();

                        const double U_density = deg * material.U_energy(grad_P_gp);
                        const double W_density = deg * material.W_energy(P_gp, strain_gp);
                        const double chi_density = material.chi_energy(P_gp);
                        const double elec_density = is_impermeable
                            ? deg * (-P_gp.dot(E_gp) - 0.5 * material.get_eps0() * E_gp.dot(E_gp))
                            : (-P_gp.dot(E_gp) - 0.5 * material.get_eps0() * E_gp.dot(E_gp));
                        const double surf_density = material.get_Gc() * (((1.0 - v_gp) * (1.0 - v_gp)) / (4.0 * material.get_kappa()) + material.get_kappa() * grad_v_gp.dot(grad_v_gp));

                        for (int i = 0; i < n_local; ++i) {
                            const int g_idx = mesh.get_node_index(e, i);
                            const double w = N(i) * dV;
                            num_U(g_idx) += U_density * w;
                            num_W(g_idx) += W_density * w;
                            num_chi(g_idx) += chi_density * w;
                            num_elec(g_idx) += elec_density * w;
                            num_surf(g_idx) += surf_density * w;
                            weight(g_idx) += w;
                        }
                    });
                }
            }
        }

        // Parallélisation de la division finale et du vrac (bulk)
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < n_nodes; ++n) {
            const double w = weight(n);
            energy_gradient.current(n) = (w > 1e-14) ? num_U(n) / w : 0.0;
            energy_elastic.current(n)  = (w > 1e-14) ? num_W(n) / w : 0.0;
            energy_landau.current(n)   = (w > 1e-14) ? num_chi(n) / w : 0.0;
            energy_electric.current(n) = (w > 1e-14) ? num_elec(n) / w : 0.0;
            energy_surface.current(n)  = (w > 1e-14) ? num_surf(n) / w : 0.0;
            
            energy_bulk.current(n) = energy_gradient.current(n) + energy_elastic.current(n)
                                    + energy_landau.current(n) + energy_electric.current(n);
        }
    }
};