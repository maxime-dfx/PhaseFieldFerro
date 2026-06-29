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
        Eigen::VectorXd Px_current;
        Eigen::VectorXd Py_current;
        Eigen::VectorXd Px_prev_iter;
        Eigen::VectorXd Py_prev_iter;
        
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        Eigen::MatrixXd local_A;
        Eigen::VectorXd local_b;
        Eigen::VectorXd Px_local;
        Eigen::VectorXd Py_local;
        Eigen::VectorXd v_local;

        double Px_0;
        double Py_0;
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager;

    public:
        // Initialization of fracture based on the configuration and mesh
        Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
        void apply_boundary_conditions(Eigen::MatrixXd& A, Eigen::VectorXd& b, const std::vector<NodeBC>& bcs);

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
};