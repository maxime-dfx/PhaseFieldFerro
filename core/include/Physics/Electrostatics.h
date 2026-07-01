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
class Mechanics;
class Fracture;

class Electrostatics {
    private:
        Eigen::VectorXd phi_current;    
        Eigen::VectorXd Ex_current;     
        Eigen::VectorXd Ey_current;
        Eigen::VectorXd Ex_prev_iter;
        Eigen::VectorXd Ey_prev_iter;
        
        double Ex_0;
        double Ey_0;
        
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager; 

        void compute_nodal_electric_field();
        void apply_boundary_conditions(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, const std::vector<NodeBC>& bcs);
        void assemble_electrostatic_system(const Polarization& polarization, const Fracture& fracture, const Math& math,
                                        std::vector<Eigen::Triplet<double>>& triplets, Eigen::VectorXd& F_global);

    public:
        Electrostatics(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
        void update_electric_potential(double time, const Polarization& polarization, const Fracture& fracture, const Math& math);
        
        const Eigen::VectorXd& get_phi() const { return phi_current; }
        const Eigen::VectorXd& get_Ex() const { return Ex_current; }
        const Eigen::VectorXd& get_Ey() const { return Ey_current; }

        double calculate_error() {
            return (Ex_current - Ex_prev_iter).norm() + (Ey_current - Ey_prev_iter).norm();
        }

        void save_previous_iteration() {
            Ex_prev_iter = Ex_current;
            Ey_prev_iter = Ey_current;
        }

        double get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const;
        double get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const;

};