#ifndef PHYSICS_ELECTROSTATICS_H
#define PHYSICS_ELECTROSTATICS_H

#include <Eigen/Core>
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"

class Polarization;
class Fracture;
class Math;
class Element;
struct GaussPoint2D;

class Electrostatics {
private:
    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager;

    Eigen::VectorXd phi_current;
    Eigen::VectorXd phi_prev_iter;
    Eigen::VectorXd phi_n;
    Eigen::VectorXd phi_backup;

    Eigen::VectorXd Ex_current;
    Eigen::VectorXd Ey_current;

public:
    Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    void update_phi(double time, const Polarization& polarization, const Fracture& fracture, const Math& math);
    void compute_electric_field();

    double get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const;
    double get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const;

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    const Eigen::VectorXd& get_Ex() const { return Ex_current; }
    const Eigen::VectorXd& get_Ey() const { return Ey_current; }
    const Eigen::VectorXd& get_phi() const { return phi_current; }
};

#endif // PHYSICS_ELECTROSTATICS_H