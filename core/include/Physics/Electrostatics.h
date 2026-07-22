#ifndef PHYSICS_ELECTROSTATICS_H
#define PHYSICS_ELECTROSTATICS_H

#include <Eigen/Core>
#include <vector>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "BC/BoundaryManager.h"
#include <Eigen/CholmodSupport>

class Polarization;
class Fracture;
class MaterialModel;
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

    // --- Fonctions de routage et d'initialisation ---
    void appliquer_conditions_initiales();

    double calculer_champ_triangle(const Element& elem, const Eigen::VectorXd& champ_nodal) const;
    double calculer_champ_quadrangle(const Element& elem, const GaussPoint2D& gp, const Eigen::VectorXd& champ_nodal) const;

    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> m_cholmod_solver;
    bool m_is_pattern_analyzed = false;
public:
    Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    void update_phi(double time, const Polarization& polarization, const Fracture& fracture, const MaterialModel& material);
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
    void set_phi(const Eigen::VectorXd& val) { phi_n = val; }
};

#endif // PHYSICS_ELECTROSTATICS_H