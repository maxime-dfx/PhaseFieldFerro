// #pragma once
// #include <Eigen/Sparse>
// #include <Eigen/SparseLU>
// #include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <cmath>
// #include <algorithm>
// #include "MultiphysicsAssembler.hpp"
// #include "BoundaryManager.hpp"
// #include "Mesh.hpp"

// using namespace Eigen;

// class StaggeredSolver {
// private:
//     MultiphysicsAssembler& assembler_;
//     BoundaryManager& bc_manager_; 
    
//     VectorXd u_;   
//     VectorXd phi_; 
//     VectorXd p_;   
//     VectorXd v_;   

//     VectorXd p_n_; 
//     VectorXd v_n_;

//     double dt_;
//     double tol_;        
//     int max_iter_;      

//     SparseLU<SparseMatrix<double>> linear_solver_;

// public:
//     StaggeredSolver(MultiphysicsAssembler& assembler, BoundaryManager& bc_manager, double dt, int total_nodes) 
//         : assembler_(assembler), bc_manager_(bc_manager), dt_(dt), tol_(1e-4), max_iter_(50) 
//     {
//         u_   = VectorXd::Zero(total_nodes * 2);
//         phi_ = VectorXd::Zero(total_nodes * 1);
//         p_   = VectorXd::Zero(total_nodes * 2);
//         v_   = VectorXd::Ones(total_nodes * 1); 

//         p_n_ = p_;
//         v_n_ = v_;
//     }

//     void run_simulation(double t_max, const Mesh& mesh, const std::string& output_dir) {
        
//         const std::vector<Element>& elements = mesh.get_elements();
        
//         double current_time = 0.0;
//         int step = 0;

//         while (current_time < t_max) {
//             std::cout << "--- Pas de temps " << step << " | t = " << current_time << " ---" << std::endl;

//             bc_manager_.update_time(current_time);

//             p_n_ = p_;
//             v_n_ = v_;

//             bool converged = false;
//             int iter = 0;

//             while (!converged && iter < max_iter_) {
//                 VectorXd p_prev_iter = p_;
//                 VectorXd v_prev_iter = v_;

//                 SparseMatrix<double> K_p;
//                 VectorXd F_p;
//                 assembler_.assemblePolarization(elements, p_n_, v_, phi_, u_, K_p, F_p);
//                 apply_dirichlet_bc(K_p, F_p, "polarization"); 
//                 linear_solver_.compute(K_p);
//                 p_ = linear_solver_.solve(F_p); 

//                 SparseMatrix<double> K_u;
//                 VectorXd F_u;
//                 assembler_.assembleMechanics(elements, p_, v_, K_u, F_u); 
//                 apply_dirichlet_bc(K_u, F_u, "mechanics");  
//                 linear_solver_.compute(K_u);
//                 u_ = linear_solver_.solve(F_u);

//                 SparseMatrix<double> K_phi;
//                 VectorXd F_phi;
//                 assembler_.assembleElectrostatics(elements, p_, v_, K_phi, F_phi);
//                 apply_dirichlet_bc(K_phi, F_phi, "electrostatics"); 
//                 linear_solver_.compute(K_phi);
//                 phi_ = linear_solver_.solve(F_phi);

//                 SparseMatrix<double> K_v;
//                 VectorXd F_v;
//                 assembler_.assembleFracture(elements, v_n_, p_, u_, K_v, F_v); 
//                 linear_solver_.compute(K_v);
//                 VectorXd v_new = linear_solver_.solve(F_v);

//                 for (int i = 0; i < v_.size(); ++i) {
//                     v_(i) = std::min(v_new(i), v_n_(i)); 
//                     v_(i) = std::max(v_(i), 0.0); 
//                 }

//                 double err_p = (p_ - p_prev_iter).norm() / (p_.norm() + 1e-12);
//                 double err_v = (v_ - v_prev_iter).norm() / (v_.norm() + 1e-12);

//                 std::cout << "  Iter " << iter << " | Err(p) = " << err_p << " | Err(v) = " << err_v << std::endl;

//                 if (err_p < tol_ && err_v < tol_) {
//                     converged = true;
//                 }
//                 iter++;
//             }

//             if (!converged) {
//                 std::cerr << "AVERTISSEMENT : Non-convergence au pas de temps " << step << "!" << std::endl;
//             }

//             export_to_vtk(step, current_time, mesh, output_dir);
//             current_time += dt_;
//             step++;
//         }
//     }

// private:
//     void apply_dirichlet_bc(SparseMatrix<double>& K, VectorXd& F, const std::string& field) {
//         const double PENALTY = 1e12; 
        
//         if (field == "electrostatics") {
//             const auto& bcs = bc_manager_.get_phi_bcs();
//             for (size_t i = 0; i < bcs.size(); ++i) {
//                 if (bcs[i].type == BCType::DIRICHLET) {
//                     K.coeffRef(i, i) = PENALTY;
//                     F(i) = PENALTY * bcs[i].value;
//                 }
//             }
//         } 
//         else if (field == "mechanics") {
//             const auto& bcs_ux = bc_manager_.get_ux_bcs();
//             const auto& bcs_uy = bc_manager_.get_uy_bcs();
//             for (size_t i = 0; i < bcs_ux.size(); ++i) {
//                 if (bcs_ux[i].type == BCType::DIRICHLET) {
//                     K.coeffRef(2*i, 2*i) = PENALTY;
//                     F(2*i) = PENALTY * bcs_ux[i].value;
//                 }
//                 if (bcs_uy[i].type == BCType::DIRICHLET) {
//                     K.coeffRef(2*i + 1, 2*i + 1) = PENALTY;
//                     F(2*i + 1) = PENALTY * bcs_uy[i].value;
//                 }
//             }
//         }
//         else if (field == "polarization") {
//             const auto& bcs_px = bc_manager_.get_px_bcs();
//             const auto& bcs_py = bc_manager_.get_py_bcs();
//             for (size_t i = 0; i < bcs_px.size(); ++i) {
//                 if (bcs_px[i].type == BCType::DIRICHLET) {
//                     K.coeffRef(2*i, 2*i) = PENALTY;
//                     F(2*i) = PENALTY * bcs_px[i].value;
//                 }
//                 if (bcs_py[i].type == BCType::DIRICHLET) {
//                     K.coeffRef(2*i + 1, 2*i + 1) = PENALTY;
//                     F(2*i + 1) = PENALTY * bcs_py[i].value;
//                 }
//             }
//         }
//     }

//     void export_to_vtk(int step, double time, const Mesh& mesh, const std::string& output_dir) {
//         std::ostringstream filename;
//         filename << output_dir << "/result_" << std::setw(4) << std::setfill('0') << step << ".vtk";

//         std::ofstream vtk_file(filename.str());
//         if (!vtk_file.is_open()) {
//             std::cerr << "ERREUR : Impossible de créer le fichier " << filename.str() << " (Le dossier existe-t-il ?)" << std::endl;
//             return;
//         }

//         const auto& nodes = mesh.get_nodes();
//         const auto& elements = mesh.get_elements();
//         int num_nodes = nodes.size();
//         int num_elements = elements.size();

//         vtk_file << "# vtk DataFile Version 3.0\n";
//         vtk_file << "Multiphysics Simulation - Step " << step << " - Time " << time << "\n";
//         vtk_file << "ASCII\n";
//         vtk_file << "DATASET UNSTRUCTURED_GRID\n";

//         vtk_file << "POINTS " << num_nodes << " double\n";
//         for (int i = 0; i < num_nodes; ++i) {
//             vtk_file << nodes[i].x << " " << nodes[i].y << " 0.0\n"; 
//         }

//         int connectivity_size = 0;
//         for (const auto& el : elements) {
//             connectivity_size += 1 + el.node_indices.size(); 
//         }

//         vtk_file << "\nCELLS " << num_elements << " " << connectivity_size << "\n";
//         for (const auto& el : elements) {
//             vtk_file << el.node_indices.size();
//             for (int idx : el.node_indices) {
//                 vtk_file << " " << idx;
//             }
//             vtk_file << "\n";
//         }

//         vtk_file << "\nCELL_TYPES " << num_elements << "\n";
//         for (const auto& el : elements) {
//             int n = el.node_indices.size();
//             int vtk_type = (n == 3) ? 5 : ((n == 4) ? 9 : 7); 
//             vtk_file << vtk_type << "\n";
//         }

//         vtk_file << "\nPOINT_DATA " << num_nodes << "\n";

//         vtk_file << "SCALARS damage double 1\n";
//         vtk_file << "LOOKUP_TABLE default\n";
//         for (int i = 0; i < num_nodes; ++i) {
//             vtk_file << v_(i) << "\n";
//         }

//         vtk_file << "SCALARS potential double 1\n";
//         vtk_file << "LOOKUP_TABLE default\n";
//         for (int i = 0; i < num_nodes; ++i) {
//             vtk_file << phi_(i) << "\n";
//         }

//         vtk_file << "VECTORS displacement double\n";
//         for (int i = 0; i < num_nodes; ++i) {
//             vtk_file << u_(2 * i) << " " << u_(2 * i + 1) << " 0.0\n";
//         }

//         vtk_file << "VECTORS polarization double\n";
//         for (int i = 0; i < num_nodes; ++i) {
//             vtk_file << p_(2 * i) << " " << p_(2 * i + 1) << " 0.0\n";
//         }

//         vtk_file.close();
//     }
// };