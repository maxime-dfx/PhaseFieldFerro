#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Core/Mesh.h"
#include "IO/Datafile.h"

struct GaussPoint2D;
class Math;
class Polarization; 
class Mechanics;
class Fracture;

class Electrostatics {
    private:
        Eigen::VectorXd Ex_current;     
        Eigen::VectorXd Ey_current;
        Eigen::VectorXd Ex_prev_iter;
        Eigen::VectorXd Ey_prev_iter;
        double Ex_0;
        double Ey_0;
        const Datafile& config;
        const Mesh& mesh;

    public:
        // Initialization of fracture based on the configuration and mesh
        Electrostatics(const Datafile& config, const Mesh& mesh);
        void update_Ex(double time, const Polarization& polarization, const Mechanics& mechanics, const Fracture& fracture);
        void update_Ey(double time, const Polarization& polarization, const Mechanics& mechanics, const Fracture& fracture);
        const Eigen::VectorXd& get_Ex() const { return Ex_current; }
        const Eigen::VectorXd& get_Ey() const { return Ey_current; }

        double calculate_error() {
            return (Ex_current - Ex_prev_iter).norm();
        }

        void save_previous_iteration() {
            Ex_prev_iter = Ex_current;
            Ey_prev_iter = Ey_current;
        }

        double get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const;
        double get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const;
};