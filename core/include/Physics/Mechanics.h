#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Core/BoundaryManager.h"

struct GaussPoint2D;
class Math;
class Polarization; 
class Fracture;
class Electrostatics;

class Mechanics {
    private:
        Eigen::VectorXd ux_current;     
        Eigen::VectorXd uy_current;
        Eigen::VectorXd ux_prev_iter;
        Eigen::VectorXd uy_prev_iter;
        Eigen::SparseMatrix<double> K_global;
        Eigen::VectorXd F_global;

        double ux_0;
        double uy_0;
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager; 

    public:
        Mechanics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
        
        void update_u(double time, const Polarization& polarization, const Fracture& fracture, const Math& math);
        
        const Eigen::VectorXd& get_ux() const { return ux_current; }
        const Eigen::VectorXd& get_uy() const { return uy_current; }

        double calculate_error() {
            double err_x = (ux_current - ux_prev_iter).norm();
            double err_y = (uy_current - uy_prev_iter).norm();
            return err_x + err_y;
        }

        void save_previous_iteration() {
            ux_prev_iter = ux_current;
            uy_prev_iter = uy_current;
        }

        Eigen::Matrix2d get_strain_at_gp(const Element& elem, const GaussPoint2D& gp, const std::vector<std::array<double, 2>>& coords) const;
        
    private:
        void apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, 
                                       const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y);
};