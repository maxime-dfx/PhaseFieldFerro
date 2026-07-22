#ifndef PHYSICS_POLARIZATION_H
#define PHYSICS_POLARIZATION_H

#include <Eigen/Core>
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "BC/BoundaryManager.h"
#include <Eigen/CholmodSupport>


class Fracture;
class Mechanics;
class Electrostatics;
class MaterialModel;

class Polarization {
private:
    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager;

    Eigen::VectorXd Px_current, Py_current;
    Eigen::VectorXd Px_prev_iter, Py_prev_iter;
    Eigen::VectorXd Px_n, Py_n;
    Eigen::VectorXd Px_backup, Py_backup;

    // --- Fonction d'initialisation isolée ---
    void appliquer_conditions_initiales();

    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> m_cholmod_solver;
    bool m_is_pattern_analyzed = false;

public:
    Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    void update_P(double time, double dt_relax,  const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const MaterialModel& material);
    void map_global_vector_to_components(const Eigen::VectorXd& P_new);

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    const Eigen::VectorXd& get_Px() const { return Px_current; }
    const Eigen::VectorXd& get_Py() const { return Py_current; }
    void set_Px(const Eigen::VectorXd& val) { Px_n = val; }
    void set_Py(const Eigen::VectorXd& val) { Py_n = val; }
};

#endif // PHYSICS_POLARIZATION_H