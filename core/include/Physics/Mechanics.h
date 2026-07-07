#pragma once
#include <Eigen/Dense>
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"
#include "Core/ElementIntegrator.h" // Pour GaussPoint2D

class Polarization; 
class Fracture;
class Math;

class Mechanics {
private:
    Eigen::VectorXd ux_current, uy_current;
    Eigen::VectorXd ux_prev_iter, uy_prev_iter;
    Eigen::VectorXd ux_backup, uy_backup;

    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager; 

    void map_global_vector_to_components(const Eigen::VectorXd& U_new);

public:
    Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
    
    void update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math);
    
    const Eigen::VectorXd& get_ux() const { return ux_current; }
    const Eigen::VectorXd& get_uy() const { return uy_current; }

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    Eigen::Matrix2d get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const;
};