#pragma once
#include <vector>
#include <Eigen/Dense>
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

    public:
        Fracture(const Datafile& config, const Mesh& mesh);
        
        void update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math);
        
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
};