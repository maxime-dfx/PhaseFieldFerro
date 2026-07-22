#ifndef PHYSICS_FRACTURE_H
#define PHYSICS_FRACTURE_H

#include <Eigen/Core>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include <Eigen/CholmodSupport>


class Polarization;
class Mechanics;
class Electrostatics;
class MaterialModel;

class Fracture {
private:
    const Datafile& config;
    const Mesh& mesh;

    Eigen::VectorXd v_current;
    Eigen::VectorXd v_prev_iter;
    Eigen::VectorXd v_n;
    Eigen::VectorXd v_backup;

    // --- Fonction d'initialisation isolée ---
    void appliquer_conditions_initiales();
    double calculer_distance_signee_precrack(double x, double y, const PrecrackConfig& pc) const;

    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> m_cholmod_solver;
    bool m_is_pattern_analyzed = false;

public:
    Fracture(const Datafile& config, const Mesh& mesh);

    void update_precrack_geometry(double time);
    void update_v(double dt_relax, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const MaterialModel& material);
    void enforce_physical_bounds(const Eigen::VectorXd& v_new);

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    const Eigen::VectorXd& get_v() const { return v_current; }
    void set_v(const Eigen::VectorXd& val) { v_n = val; }
};

#endif // PHYSICS_FRACTURE_H