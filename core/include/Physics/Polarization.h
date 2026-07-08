#pragma once

#include <Eigen/Dense>
#include "IO/Datafile.h"
#include "Core/Mesh.h"
#include "Core/BoundaryManager.h"

// Forward declarations pour eviter les inclusions circulaires
class Fracture;
class Mechanics;
class Electrostatics;
class Math;

class Polarization {
public:
    Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    // ==========================================================================
    // NOUVELLE METHODE MONOLITHIQUE : resout Px et Py EN UNE SEULE FOIS,
    // dans un unique systeme lineaire de taille 2*n_dof (DOFs interleaves
    // comme pour Mechanics::update_u), au lieu de deux resolutions
    // sequentielles update_Px() puis update_Py().
    // ==========================================================================
    void update_polarization(double time, const Fracture& fracture, const Mechanics& mechanics,
                              const Electrostatics& electrostatics, const Math& math);

    // --- Anciennes methodes conservees pour compatibilite ascendante ---
    // (peuvent etre supprimees une fois la migration validee)
    void update_Px(double time, const Fracture& fracture, const Mechanics& mechanics,
                    const Electrostatics& electrostatics, const Math& math) {
        update_polarization_component(time, fracture, mechanics, electrostatics, math, 0);
    }
    void update_Py(double time, const Fracture& fracture, const Mechanics& mechanics,
                    const Electrostatics& electrostatics, const Math& math) {
        update_polarization_component(time, fracture, mechanics, electrostatics, math, 1);
    }

    const Eigen::VectorXd& get_Px() const { return Px_current; }
    const Eigen::VectorXd& get_Py() const { return Py_current; }

    double calculate_error() const {
        return (Px_current - Px_prev_iter).norm() + (Py_current - Py_prev_iter).norm();
    }

    void save_previous_iteration() {
        Px_prev_iter = Px_current;
        Py_prev_iter = Py_current;
    }

    void save_previous_state() {
        Px_backup = Px_current;
        Py_backup = Py_current;
    }

    void restore_previous_state() {
        Px_current = Px_backup;
        Py_current = Py_backup;
    }

    void freeze_time_step() {
        Px_n = Px_current;
        Py_n = Py_current;
    }

    void update_history() {
        // Rien pour le moment, comme dans Mechanics
    }

private:
    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager;

    Eigen::VectorXd Px_current, Py_current;
    Eigen::VectorXd Px_prev_iter, Py_prev_iter;
    Eigen::VectorXd Px_n, Py_n;
    Eigen::VectorXd Px_backup, Py_backup;

    // Conservee pour les anciennes methodes update_Px/update_Py (compatibilite)
    void update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics,
                                        const Electrostatics& electrostatics, const Math& math, int component);

    void map_global_vector_to_components(const Eigen::VectorXd& P_new);
};