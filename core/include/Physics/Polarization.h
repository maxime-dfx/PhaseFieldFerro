#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Core/BoundaryManager.h"

class Math;
class Fracture;
class Mechanics;
class Electrostatics;

class Polarization {
private:
    Eigen::VectorXd Px_current, Py_current;
    Eigen::VectorXd Px_prev_iter, Py_prev_iter;
    Eigen::VectorXd Px_n, Py_n;
    Eigen::VectorXd Px_backup, Py_backup;

    const Datafile& config;
    const Mesh& mesh;
    const BoundaryManager& bc_manager;

public:
    Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);

    void update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, int component);
    
    void update_Px(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
        update_polarization_component(time, fracture, mechanics, electrostatics, math, 0);
    }

    void update_Py(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
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

    void freeze_time_step() {
        Px_n = Px_current;
        Py_n = Py_current;
    }

    void save_previous_state() {
        if (Px_backup.size() != Px_current.size()) {
            Px_backup.resize(Px_current.size());
            Py_backup.resize(Py_current.size());
        }
        Px_backup = Px_current;
        Py_backup = Py_current;
    }

    void restore_previous_state() {
        Px_current = Px_backup;
        Py_current = Py_backup;
    }

    void update_history() {
        Px_n = Px_current;
        Py_n = Py_current;
    }
};