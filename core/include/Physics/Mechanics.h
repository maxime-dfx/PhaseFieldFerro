#ifndef PHYSICS_MECHANICS_H
#define PHYSICS_MECHANICS_H

#include <Eigen/Core>
#include <vector>
#include <array>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "BC/BoundaryManager.h"
#include <Eigen/CholmodSupport>


class Polarization;
class Fracture;
class MaterialModel;
class Element;
struct GaussPoint2D;

class Mechanics {
private:
    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager;

    Eigen::VectorXd ux_current, uy_current;
    Eigen::VectorXd ux_prev_iter, uy_prev_iter;
    Eigen::VectorXd ux_backup, uy_backup;
    Eigen::VectorXd sigma_xx_current, sigma_yy_current, sigma_xy_current;


    // --- Nouvelles fonctions privées spécialisées (Travailleurs) ---
    void appliquer_conditions_initiales();

    Eigen::Matrix2d calculer_deformation_triangle(
        const Element& elem, 
        const std::array<std::array<double, 2>, 8>& coords) const;

    Eigen::Matrix2d calculer_deformation_quadrangle(
        const Element& elem, 
        const GaussPoint2D& gp, 
        const std::array<std::array<double, 2>, 8>& coords) const;
        
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> m_cholmod_solver;
    bool m_is_pattern_analyzed = false;
public:
    Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    void update_u(double time, const Polarization& polarization, const Fracture& fracture, const MaterialModel& material);
    void map_global_vector_to_components(const Eigen::VectorXd& U_new);

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    Eigen::Matrix2d get_strain_at_gp(
        const Element& elem, 
        const GaussPoint2D& gp, 
        const std::array<std::array<double, 2>, 8>& coords) const;

    Eigen::Matrix2d get_stress_at_gp(
        const Element& elem, 
        const GaussPoint2D& gp,
        const std::array<std::array<double, 2>, 8>& coords,
        const Polarization& polarization,
        const Fracture& fracture,
        const MaterialModel& material) const;
    
    const Eigen::VectorXd& get_ux() const { return ux_current; }
    const Eigen::VectorXd& get_uy() const { return uy_current; }
    void set_ux(const Eigen::VectorXd& ux) { ux_current = ux; }
    void set_uy(const Eigen::VectorXd& uy) { uy_current = uy; }

    // Projection nodale du tenseur des contraintes (L2), à appeler après update_u()
    void compute_stress_field(const Polarization& polarization, const Fracture& fracture, const MaterialModel& material);

    const Eigen::VectorXd& get_sigma_xx() const { return sigma_xx_current; }
    const Eigen::VectorXd& get_sigma_yy() const { return sigma_yy_current; }
    const Eigen::VectorXd& get_sigma_xy() const { return sigma_xy_current; }

};

#endif // PHYSICS_MECHANICS_H