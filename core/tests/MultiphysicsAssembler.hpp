#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "FerroelectricMaterial.hpp"
#include "Element.hpp" 
using namespace Eigen;
using TripletList = std::vector<Triplet<double>>;

class MultiphysicsAssembler {
private:
    FerroelectricMaterial material_;
    double dt_; 

public:
    MultiphysicsAssembler(const FerroelectricMaterial& mat, double dt) 
        : material_(mat), dt_(dt) {}

    void assemblePolarization(
        const std::vector<Element>& elements,
        const VectorXd& p_n,     
        const VectorXd& v_prev,  
        const VectorXd& phi_prev,
        const VectorXd& u_prev,  
        SparseMatrix<double>& K_global,
        VectorXd& F_global) 
    {
        int num_dofs = p_n.size();
        K_global.resize(num_dofs, num_dofs);
        TripletList triplets;
        triplets.reserve(elements.size() * 36); 
        F_global = VectorXd::Zero(num_dofs);

        for (const auto& el : elements) {
            int nodes = el.node_indices.size();
            MatrixXd K_elem = MatrixXd::Zero(nodes * 2, nodes * 2);
            VectorXd F_elem = VectorXd::Zero(nodes * 2);

            for (int q = 0; q < el.num_gauss_points(); ++q) {
                double w_detJ = el.weight(q);
                VectorXd N = el.N(q);
                MatrixXd grad_N = el.grad_N(q); 
                
                double v_local = interpolate_scalar(N, v_prev, el.node_indices);
                Vector2d p_local = interpolate_vector(N, p_n, el.node_indices);
                Vector3d eps_local = el.B_mech(q) * extract_local(u_prev, el.node_indices, 2);
                Vector2d E_local = -grad_N * extract_local(phi_prev, el.node_indices, 1);

                double deg = material_.degradation(v_local);

                MatrixXd N_mat = expand_N_to_vector_dof(N);             
                MatrixXd B_mat = expand_gradN_to_vector_dof(grad_N);    

                K_elem += ((material_.mu_p() / dt_) * N_mat.transpose() * N_mat + 
                           deg * material_.a0() * B_mat.transpose() * B_mat) * w_detJ;
                
                Vector2d force_meca   = deg * material_.mechanical_driving_force_on_p(p_local, eps_local);
                Vector2d force_landau = -material_.landau_derivative(p_local);
                Vector2d inertia      = (material_.mu_p() / dt_) * p_local;
                Vector2d total_thermo_force = inertia + force_meca + force_landau + E_local;

                F_elem += N_mat.transpose() * total_thermo_force * w_detJ;
            }
            add_to_global_system(triplets, F_global, el.node_indices, K_elem, F_elem, 2);
        }
        K_global.setFromTriplets(triplets.begin(), triplets.end());
    }

    void assembleFracture(
        const std::vector<Element>& elements,
        const VectorXd& v_n,     
        const VectorXd& p_prev,  
        const VectorXd& u_prev,  
        SparseMatrix<double>& K_global,
        VectorXd& F_global) 
    {
        int num_dofs = v_n.size();
        K_global.resize(num_dofs, num_dofs);
        TripletList triplets;
        F_global = VectorXd::Zero(num_dofs);

        double Gc = material_.Gc();
        double kappa = material_.kappa();
        double mu_v = material_.mu_v();

        for (const auto& el : elements) {
            int nodes = el.node_indices.size();
            MatrixXd K_elem = MatrixXd::Zero(nodes, nodes);
            VectorXd F_elem = VectorXd::Zero(nodes);

            for (int q = 0; q < el.num_gauss_points(); ++q) {
                double w_detJ = el.weight(q);
                VectorXd N = el.N(q);
                MatrixXd B = el.grad_N(q); 
                
                double v_local = interpolate_scalar(N, v_n, el.node_indices);
                Vector2d p_local = interpolate_vector(N, p_prev, el.node_indices);
                Matrix2d grad_p_local = interpolate_gradient(B, p_prev, el.node_indices);
                Vector3d eps_local = el.B_mech(q) * extract_local(u_prev, el.node_indices, 2);

                K_elem += ( (mu_v / dt_ + Gc / (2.0 * kappa)) * N * N.transpose() + 
                            (2.0 * Gc * kappa) * B.transpose() * B ) * w_detJ;

                double driving_energy = material_.crack_driving_force(p_local, eps_local, grad_p_local);
                double rhs_scalar = (mu_v / dt_) * v_local + (Gc / (2.0 * kappa)) - 2.0 * v_local * driving_energy;

                F_elem += N * rhs_scalar * w_detJ;
            }
            add_to_global_system(triplets, F_global, el.node_indices, K_elem, F_elem, 1);
        }
        K_global.setFromTriplets(triplets.begin(), triplets.end());
    }
    
    void assembleMechanics(
        const std::vector<Element>& elements,
        const VectorXd& p_prev,  
        const VectorXd& v_prev,  
        SparseMatrix<double>& K_global,
        VectorXd& F_global) 
    {
        int num_dofs = v_prev.size() * 2; 
        K_global.resize(num_dofs, num_dofs);
        TripletList triplets;
        F_global = VectorXd::Zero(num_dofs);
        Matrix3d C = material_.C(); 

        for (const auto& el : elements) {
            int nodes = el.node_indices.size();
            MatrixXd K_elem = MatrixXd::Zero(nodes * 2, nodes * 2);
            VectorXd F_elem = VectorXd::Zero(nodes * 2);

            for (int q = 0; q < el.num_gauss_points(); ++q) {
                double w_detJ = el.weight(q);
                VectorXd N = el.N(q);
                MatrixXd B_u = el.B_mech(q); 
                
                double v_local = interpolate_scalar(N, v_prev, el.node_indices);
                Vector2d p_local = interpolate_vector(N, p_prev, el.node_indices);
                double deg = material_.degradation(v_local);
                Vector3d eps0 = material_.spontaneous_strain(p_local);

                Matrix3d C_deg = deg * C;
                K_elem += B_u.transpose() * C_deg * B_u * w_detJ;
                F_elem += B_u.transpose() * (C_deg * eps0) * w_detJ;
            }
            add_to_global_system(triplets, F_global, el.node_indices, K_elem, F_elem, 2);
        }
        K_global.setFromTriplets(triplets.begin(), triplets.end());
    }

    void assembleElectrostatics(
        const std::vector<Element>& elements,
        const VectorXd& p_prev,  
        const VectorXd& v_prev,  
        SparseMatrix<double>& K_global,
        VectorXd& F_global) 
    {
        int num_dofs = v_prev.size(); 
        K_global.resize(num_dofs, num_dofs);
        TripletList triplets;
        F_global = VectorXd::Zero(num_dofs);
        double eps0_const = material_.eps0();

        for (const auto& el : elements) {
            int nodes = el.node_indices.size();
            MatrixXd K_elem = MatrixXd::Zero(nodes, nodes);
            VectorXd F_elem = VectorXd::Zero(nodes);

            for (int q = 0; q < el.num_gauss_points(); ++q) {
                double w_detJ = el.weight(q);
                VectorXd N = el.N(q);
                MatrixXd B_phi = el.grad_N(q); 
                Vector2d p_local = interpolate_vector(N, p_prev, el.node_indices);

                K_elem += B_phi.transpose() * eps0_const * B_phi * w_detJ;
                F_elem += B_phi.transpose() * p_local * w_detJ;
            }
            add_to_global_system(triplets, F_global, el.node_indices, K_elem, F_elem, 1);
        }
        K_global.setFromTriplets(triplets.begin(), triplets.end());
    }

private:
    inline VectorXd extract_local(const VectorXd& global_vec, const std::vector<int>& indices, int dofs) const {
        VectorXd local(indices.size() * dofs);
        for (size_t i = 0; i < indices.size(); ++i) {
            for (int d = 0; d < dofs; ++d) {
                local(i * dofs + d) = global_vec(indices[i] * dofs + d);
            }
        }
        return local;
    }

    inline double interpolate_scalar(const VectorXd& N, const VectorXd& global_vec, const std::vector<int>& indices) const {
        double val = 0.0;
        for (size_t i = 0; i < indices.size(); ++i) { val += N(i) * global_vec(indices[i]); }
        return val;
    }

    inline Vector2d interpolate_vector(const VectorXd& N, const VectorXd& global_vec, const std::vector<int>& indices) const {
        Vector2d val = Vector2d::Zero();
        for (size_t i = 0; i < indices.size(); ++i) {
            val(0) += N(i) * global_vec(indices[i] * 2);
            val(1) += N(i) * global_vec(indices[i] * 2 + 1);
        }
        return val;
    }

    inline Matrix2d interpolate_gradient(const MatrixXd& B, const VectorXd& global_vec, const std::vector<int>& indices) const {
        Matrix2d grad = Matrix2d::Zero();
        for (size_t i = 0; i < indices.size(); ++i) {
            grad(0, 0) += B(0, i) * global_vec(indices[i] * 2);    
            grad(0, 1) += B(1, i) * global_vec(indices[i] * 2);    
            grad(1, 0) += B(0, i) * global_vec(indices[i] * 2 + 1); 
            grad(1, 1) += B(1, i) * global_vec(indices[i] * 2 + 1); 
        }
        return grad;
    }

    inline MatrixXd expand_N_to_vector_dof(const VectorXd& N) const {
        MatrixXd N_mat = MatrixXd::Zero(2, 2 * N.size());
        for (int i = 0; i < N.size(); ++i) {
            N_mat(0, 2 * i) = N(i);
            N_mat(1, 2 * i + 1) = N(i);
        }
        return N_mat;
    }

    inline MatrixXd expand_gradN_to_vector_dof(const MatrixXd& grad_N) const {
        MatrixXd B_mat = MatrixXd::Zero(4, 2 * grad_N.cols());
        for (int i = 0; i < grad_N.cols(); ++i) {
            B_mat(0, 2 * i)     = grad_N(0, i); 
            B_mat(1, 2 * i)     = grad_N(1, i); 
            B_mat(2, 2 * i + 1) = grad_N(0, i); 
            B_mat(3, 2 * i + 1) = grad_N(1, i); 
        }
        return B_mat;
    }

    inline void add_to_global_system(TripletList& triplets, VectorXd& F_global, const std::vector<int>& indices, 
                                     const MatrixXd& K_elem, const VectorXd& F_elem, int dofs) const {
        for (size_t i = 0; i < indices.size(); ++i) {
            for (int di = 0; di < dofs; ++di) {
                F_global(indices[i] * dofs + di) += F_elem(i * dofs + di);
                for (size_t j = 0; j < indices.size(); ++j) {
                    for (int dj = 0; dj < dofs; ++dj) {
                        triplets.push_back(Triplet<double>(indices[i] * dofs + di, indices[j] * dofs + dj, K_elem(i * dofs + di, j * dofs + dj)));
                    }
                }
            }
        }
    }
};