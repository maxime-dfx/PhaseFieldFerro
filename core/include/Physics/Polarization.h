#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Core/BoundaryManager.h"

class Math;
class Fracture;
class Mechanics;
class Electrostatics;

class Polarization {
    private:
        Eigen::VectorXd Px_current;
        Eigen::VectorXd Py_current;
        Eigen::VectorXd Px_prev_iter;
        Eigen::VectorXd Py_prev_iter;
        Eigen::VectorXd Px_n;
        Eigen::VectorXd Py_n;
        
        Eigen::VectorXd b;
        Eigen::SparseMatrix<double> local_A;
        Eigen::VectorXd local_b;
        Eigen::VectorXd Px_local;
        Eigen::VectorXd Py_local;
<<<<<<< HEAD

        Eigen::VectorXd Px_backup;
        Eigen::VectorXd Py_backup;
=======
        Eigen::VectorXd v_local;
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

        double Px_0;
        double Py_0;
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager;
<<<<<<< HEAD
        void assemble_component_system(double dt, const Fracture& fracture, const Mechanics& mechanics, 
                                    const Electrostatics& electrostatics, const Math& math, int component,
                                    std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b, Eigen::VectorXd& diag_check);
        void handle_floating_dofs(std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& b, const Eigen::VectorXd& diag_check, int component);
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da

    public:
        // Initialization of fracture based on the configuration and mesh
        Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
        void apply_boundary_conditions(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const std::vector<NodeBC>& bcs);

        void update_polarization_component(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math, int component);
        
        void update_Px(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
            update_polarization_component(time, fracture, mechanics, electrostatics, math, 0);
        }

        void update_Py(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
            update_polarization_component(time, fracture, mechanics, electrostatics, math, 1);
        }
        const Eigen::VectorXd& get_Px() const { return Px_current; }
        const Eigen::VectorXd& get_Py() const { return Py_current; }

        double calculate_error() {
            return (Px_current - Px_prev_iter).norm();
        }

        void save_previous_iteration() {
            Px_prev_iter = Px_current;
            Py_prev_iter = Py_current;
        }

        void freeze_time_step() {
            Px_n = Px_current;
            Py_n = Py_current;
        }
<<<<<<< HEAD
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
=======
>>>>>>> 1b56e6de1054186eb666ba43dbeb72efb8eda2da
};