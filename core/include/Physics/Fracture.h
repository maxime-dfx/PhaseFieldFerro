#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Core/Mesh.h"
#include "IO/Datafile.h"

class Math;
class Polarization; 
class Mechanics;
class Electrostatics;

class Fracture {
    private:
        Eigen::VectorXd v_current;     
        Eigen::VectorXd v_prev_iter;
        Eigen::VectorXd v_n;
        const Datafile& config;
        const Mesh& mesh;       

        // Solveur persistant : analyzePattern() une seule fois. Le pattern reste
        // stable meme avec le blocage des DOFs irreversibles (coeffRef sur la
        // diagonale, deja presente via les termes de masse/diffusion).
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
        bool pattern_analyzed_ = false;

    public:
        Fracture(const Datafile& config, const Mesh& mesh);
        
        void update_v(double dt_relax, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math);
        
        const Eigen::VectorXd& get_v() const { return v_current; }
        const Eigen::VectorXd& get_v_prev() const { return v_prev_iter; }

        double calculate_error() {
            return (v_current - v_prev_iter).norm();
        }

        void save_previous_iteration() {
            v_prev_iter = v_current;
        }

        void freeze_time_step() {
            v_n = v_current;
        }
        void set_state(const Eigen::VectorXd& v, const Eigen::VectorXd& v_n_in) {
            v_current = v;
            v_n = v_n_in;
            v_prev_iter = v;
        }

};