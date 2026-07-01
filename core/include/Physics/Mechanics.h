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
        Eigen::VectorXd ux_backup;
        Eigen::VectorXd uy_backup;

        double ux_0;
        double uy_0;
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager; 
        
        void assemble_mechanical_system(const Polarization& polarization, const Fracture& fracture, const Math& math,
                                        std::vector<Eigen::Triplet<double>>& triplets);
        void map_global_vector_to_components(const Eigen::VectorXd& U_new);

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
        void save_previous_state() {
            if (ux_backup.size() != ux_current.size()) {
                ux_backup.resize(ux_current.size());
                uy_backup.resize(uy_current.size());
            }
            ux_backup = ux_current;
            uy_backup = uy_current;
        }
        void restore_previous_state() {
            ux_current = ux_backup;
            uy_current = uy_backup;
        }
        void update_history() {
        }
        
    private:
        void apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, 
                                       const std::vector<NodeBC>& bcs_x, const std::vector<NodeBC>& bcs_y);

};