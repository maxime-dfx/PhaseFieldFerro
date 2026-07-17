#include "IO/Diagnostics.h"

#include "Physics/Polarization.h"
#include "Physics/Mechanics.h"
#include "Physics/Fracture.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "Core/ShapeFunctions.h"
#include "Core/Quadrature.h"

#include <algorithm>
#include <limits>
#include <iomanip>

Diagnostics::Diagnostics(const Datafile& config, const Mesh& mesh)
    : config(config), mesh(mesh) {
    // Initialisation des tailles pour les énergies nodales (Export VTK)
    size_t n_nodes = mesh.get_num_nodes();
    U_nodal.setZero(n_nodes);
    W_nodal.setZero(n_nodes);
    chi_nodal.setZero(n_nodes);
    elec_nodal.setZero(n_nodes);
    surf_nodal.setZero(n_nodes);
}

void Diagnostics::integrate_bulk_terms(const Polarization& polarization,
                                        const Mechanics& mechanics,
                                        const Fracture& fracture,
                                        const Electrostatics& electrostatics,
                                        const Math& math,
                                        EnergyRecord& rec) const {
    double U_sum = 0.0, W_sum = 0.0, chi_sum = 0.0, elec_sum = 0.0;
    const bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    const auto& elements = mesh.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd Px_loc(n_nodes), Py_loc(n_nodes), v_loc(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            Px_loc(i) = polarization.get_Px()[indices[i]];
            Py_loc(i) = polarization.get_Py()[indices[i]];
            v_loc(i)  = fracture.get_v()[indices[i]];
        }

        auto integrate_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N,
                                 double dV, const GaussPoint2D& gp) {
            double v_gp = N.dot(v_loc);
            double phase_factor = (v_gp * v_gp) + math.eta_k;

            Eigen::Vector2d Pi_gp(N.dot(Px_loc), N.dot(Py_loc));
            Eigen::Matrix2d Pij_gp;
            Pij_gp(0, 0) = grad_N.row(0).dot(Px_loc); 
            Pij_gp(0, 1) = grad_N.row(1).dot(Px_loc); 
            Pij_gp(1, 0) = grad_N.row(0).dot(Py_loc); 
            Pij_gp(1, 1) = grad_N.row(1).dot(Py_loc); 

            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            double U   = math.U_energy(Pij_gp);
            double W   = math.W_energy(Pi_gp, eps_gp);
            double chi = math.chi_energy(Pi_gp);

            U_sum   += phase_factor * U * dV;
            W_sum   += phase_factor * W * dV;
            chi_sum += chi * dV; 

            double Ex = electrostatics.get_Ex_at_gp(elem, gp);
            double Ey = electrostatics.get_Ey_at_gp(elem, gp);
            Eigen::Vector2d E_gp(Ex, Ey);
            double elec = -0.5 * math.eps0 * E_gp.squaredNorm() - E_gp.dot(Pi_gp);

            if (is_impermeable) {
                elec_sum += phase_factor * elec * dV;
            } else {
                elec_sum += elec * dV;
            }
        };

        if (n_nodes == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

            double dV = std::abs(detJ) / 2.0;
            Eigen::RowVectorXd N(3); N << N_std[0], N_std[1], N_std[2];
            Eigen::MatrixXd grad_N(2, 3);
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }
            GaussPoint2D gp_fake{xi, eta, 1.0};
            integrate_gp(N, grad_N, dV, gp_fake);

        } else if (n_nodes == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

                double dV = gp.weight * std::abs(detJ);
                Eigen::RowVectorXd N(4); N << N_std[0], N_std[1], N_std[2], N_std[3];
                Eigen::MatrixXd grad_N(2, 4);
                for (int i = 0; i < 4; ++i) {
                    grad_N(0, i) = dN_xy[0][i];
                    grad_N(1, i) = dN_xy[1][i];
                }
                integrate_gp(N, grad_N, dV, gp);
            }
        }
    }

    rec.U_total = U_sum;
    rec.W_total = W_sum;
    rec.chi_total = chi_sum;
    rec.electric_total = elec_sum;
    rec.bulk_enthalpy = U_sum + W_sum + chi_sum + elec_sum;
}

void Diagnostics::integrate_surface_term(const Fracture& fracture, EnergyRecord& rec) const {
    double surf_sum = 0.0;
    const double Gc = config.material.Gc;
    const double kappa = config.material.kappa;
    const auto& elements = mesh.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd v_loc(n_nodes);
        for (int i = 0; i < n_nodes; ++i) v_loc(i) = fracture.get_v()[indices[i]];

        auto integrate_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N, double dV) {
            double v_gp = N.dot(v_loc);
            Eigen::Vector2d grad_v_gp = grad_N * v_loc;
            double density = (1.0 - v_gp) * (1.0 - v_gp) / (4.0 * kappa) + kappa * grad_v_gp.squaredNorm();
            surf_sum += Gc * density * dV;
        };

        if (n_nodes == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

            double dV = std::abs(detJ) / 2.0;
            Eigen::RowVectorXd N(3); N << N_std[0], N_std[1], N_std[2];
            Eigen::MatrixXd grad_N(2, 3);
            for (int i = 0; i < 3; ++i) {
                grad_N(0, i) = dN_xy[i][0];
                grad_N(1, i) = dN_xy[i][1];
            }
            integrate_gp(N, grad_N, dV);

        } else if (n_nodes == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

                double dV = gp.weight * std::abs(detJ);
                Eigen::RowVectorXd N(4); N << N_std[0], N_std[1], N_std[2], N_std[3];
                Eigen::MatrixXd grad_N(2, 4);
                for (int i = 0; i < 4; ++i) {
                    grad_N(0, i) = dN_xy[0][i];
                    grad_N(1, i) = dN_xy[1][i];
                }
                integrate_gp(N, grad_N, dV);
            }
        }
    }

    rec.surface_energy = surf_sum;
}

// =========================================================================
// LA FONCTION ORIGINALE QUI AVAIT DISPARU EST ICI !
// =========================================================================
void Diagnostics::compute_nodal_extrema(const Polarization& polarization,
                                         const Mechanics& mechanics,
                                         const Fracture& fracture,
                                         const Electrostatics& electrostatics,
                                         EnergyRecord& rec) const {
    const auto& v = fracture.get_v();
    const auto& phi = electrostatics.get_phi();
    const auto& Px = polarization.get_Px();
    const auto& Py = polarization.get_Py();
    const auto& ux = mechanics.get_ux();
    const auto& uy = mechanics.get_uy();

    rec.v_min = v.minCoeff();   rec.v_max = v.maxCoeff();
    rec.phi_min = phi.minCoeff(); rec.phi_max = phi.maxCoeff();
    rec.Px_min = Px.minCoeff(); rec.Px_max = Px.maxCoeff();
    rec.Py_min = Py.minCoeff(); rec.Py_max = Py.maxCoeff();
    rec.ux_max_abs = ux.cwiseAbs().maxCoeff();
    rec.uy_max_abs = uy.cwiseAbs().maxCoeff();

    const double alpha = 2e-2; 
    double area_cracked = 0.0;
    const auto& elements = mesh.get_elements();
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        double v_avg = 0.0;
        for (int i = 0; i < n_nodes; ++i) v_avg += v[indices[i]];
        v_avg /= n_nodes;

        if (v_avg < alpha) {
            double area = 0.0;
            if (n_nodes == 3) {
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
                area = std::abs(detJ) / 2.0;
            } else if (n_nodes == 4) {
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(0.0, 0.0);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                area = 4.0 * std::abs(detJ); 
            }
            area_cracked += area;
        }
    }
    rec.crack_length_proxy = area_cracked;
}

EnergyRecord Diagnostics::record(double load_step, double time,
                                  const Polarization& polarization,
                                  const Mechanics& mechanics,
                                  const Fracture& fracture,
                                  const Electrostatics& electrostatics,
                                  const Math& math) {
    EnergyRecord rec;
    rec.load_step = load_step;
    rec.time = time;

    integrate_bulk_terms(polarization, mechanics, fracture, electrostatics, math, rec);
    integrate_surface_term(fracture, rec);
    compute_nodal_extrema(polarization, mechanics, fracture, electrostatics, rec);

    rec.total_energy = rec.bulk_enthalpy + rec.surface_energy;

    history.push_back(rec);
    return rec;
}

void Diagnostics::write_csv(const std::string& path) const {
    std::ofstream file(path);
    file << "load_step,time,U_total,W_total,chi_total,electric_total,"
         << "bulk_enthalpy,surface_energy,total_energy,"
         << "v_min,v_max,crack_length_proxy,"
         << "phi_min,phi_max,Px_min,Px_max,Py_min,Py_max,ux_max_abs,uy_max_abs\n";

    file << std::setprecision(10);
    for (const auto& r : history) {
        file << r.load_step << "," << r.time << ","
             << r.U_total << "," << r.W_total << "," << r.chi_total << "," << r.electric_total << ","
             << r.bulk_enthalpy << "," << r.surface_energy << "," << r.total_energy << ","
             << r.v_min << "," << r.v_max << "," << r.crack_length_proxy << ","
             << r.phi_min << "," << r.phi_max << ","
             << r.Px_min << "," << r.Px_max << "," << r.Py_min << "," << r.Py_max << ","
             << r.ux_max_abs << "," << r.uy_max_abs << "\n";
    }
}

void Diagnostics::append_csv(const std::string& path) {
    if (history.empty()) return;
    std::ofstream file(path, std::ios::app);

    if (!append_header_written) {
        std::ifstream check(path);
        bool empty = check.peek() == std::ifstream::traits_type::eof();
        check.close();
        if (empty) {
            file << "load_step,time,U_total,W_total,chi_total,electric_total,"
                 << "bulk_enthalpy,surface_energy,total_energy,"
                 << "v_min,v_max,crack_length_proxy,"
                 << "phi_min,phi_max,Px_min,Px_max,Py_min,Py_max,ux_max_abs,uy_max_abs\n";
        }
        append_header_written = true;
    }

    const auto& r = history.back();
    file << std::setprecision(10)
         << r.load_step << "," << r.time << ","
         << r.U_total << "," << r.W_total << "," << r.chi_total << "," << r.electric_total << ","
         << r.bulk_enthalpy << "," << r.surface_energy << "," << r.total_energy << ","
         << r.v_min << "," << r.v_max << "," << r.crack_length_proxy << ","
         << r.phi_min << "," << r.phi_max << ","
         << r.Px_min << "," << r.Px_max << "," << r.Py_min << "," << r.Py_max << ","
         << r.ux_max_abs << "," << r.uy_max_abs << "\n";
}

void Diagnostics::compute_nodal_energies(const Polarization& polarization,
                                         const Mechanics& mechanics,
                                         const Fracture& fracture,
                                         const Electrostatics& electrostatics,
                                         const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    U_nodal.setZero(n_nodes);
    W_nodal.setZero(n_nodes);
    chi_nodal.setZero(n_nodes);
    elec_nodal.setZero(n_nodes);
    surf_nodal.setZero(n_nodes);
    Eigen::VectorXd M_lumped = Eigen::VectorXd::Zero(n_nodes);

    const bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);
    const double Gc = config.material.Gc;
    const double kappa = config.material.kappa;
    const auto& elements = mesh.get_elements();

    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& elem = elements[elem_idx];
        int n_elem_nodes = elem.get_num_nodes();
        auto coords = mesh.get_element_coords(elem_idx);
        const auto& indices = elem.get_node_indices();

        Eigen::VectorXd Px_loc(n_elem_nodes), Py_loc(n_elem_nodes), v_loc(n_elem_nodes);
        for (int i = 0; i < n_elem_nodes; ++i) {
            Px_loc(i) = polarization.get_Px()[indices[i]];
            Py_loc(i) = polarization.get_Py()[indices[i]];
            v_loc(i)  = fracture.get_v()[indices[i]];
        }

        auto process_gp = [&](const Eigen::RowVectorXd& N, const Eigen::MatrixXd& grad_N,
                              double dV, const GaussPoint2D& gp) {
            double v_gp = N.dot(v_loc);
            double phase_factor = (v_gp * v_gp) + math.eta_k;

            Eigen::Vector2d Pi_gp(N.dot(Px_loc), N.dot(Py_loc));
            Eigen::Matrix2d Pij_gp;
            Pij_gp(0, 0) = grad_N.row(0).dot(Px_loc); Pij_gp(0, 1) = grad_N.row(1).dot(Px_loc);
            Pij_gp(1, 0) = grad_N.row(0).dot(Py_loc); Pij_gp(1, 1) = grad_N.row(1).dot(Py_loc);

            Eigen::Matrix2d eps_gp = mechanics.get_strain_at_gp(elem, gp, coords);

            double U_gp = phase_factor * math.U_energy(Pij_gp);
            double W_gp = phase_factor * math.W_energy(Pi_gp, eps_gp);
            double chi_gp = math.chi_energy(Pi_gp); 

            double Ex = electrostatics.get_Ex_at_gp(elem, gp);
            double Ey = electrostatics.get_Ey_at_gp(elem, gp);
            Eigen::Vector2d E_gp(Ex, Ey);
            double elec_gp = -0.5 * math.eps0 * E_gp.squaredNorm() - E_gp.dot(Pi_gp);
            if (is_impermeable) elec_gp *= phase_factor;

            Eigen::Vector2d grad_v_gp = grad_N * v_loc;
            double density = (1.0 - v_gp) * (1.0 - v_gp) / (4.0 * kappa) + kappa * grad_v_gp.squaredNorm();
            double surf_gp = Gc * density;

            for (int i = 0; i < n_elem_nodes; ++i) {
                double Ni_dV = N[i] * dV;
                U_nodal(indices[i]) += U_gp * Ni_dV;
                W_nodal(indices[i]) += W_gp * Ni_dV;
                chi_nodal(indices[i]) += chi_gp * Ni_dV;
                elec_nodal(indices[i]) += elec_gp * Ni_dV;
                surf_nodal(indices[i]) += surf_gp * Ni_dV;
                M_lumped(indices[i]) += Ni_dV;
            }
        };

        if (n_elem_nodes == 3) {
            double xi = 1.0 / 3.0, eta = 1.0 / 3.0;
            auto N_std = ShapeFunctions::get_shape_functions_tri(xi, eta);
            auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives_tri(coords);
            if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

            double dV = std::abs(detJ) / 2.0;
            Eigen::RowVectorXd N(3); N << N_std[0], N_std[1], N_std[2];
            Eigen::MatrixXd grad_N(2, 3);
            for (int i = 0; i < 3; ++i) { grad_N(0, i) = dN_xy[i][0]; grad_N(1, i) = dN_xy[i][1]; }
            
            GaussPoint2D gp_fake{xi, eta, 1.0};
            process_gp(N, grad_N, dV, gp_fake);

        } else if (n_elem_nodes == 4) {
            const auto& gauss_points = Quadrature::get_gauss_2x2();
            for (const auto& gp : gauss_points) {
                auto N_std = ShapeFunctions::get_shape_functions(gp.xi, gp.eta);
                auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
                auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
                if (!std::isfinite(detJ) || std::abs(detJ) <= 1e-12) continue;

                double dV = gp.weight * std::abs(detJ);
                Eigen::RowVectorXd N(4); N << N_std[0], N_std[1], N_std[2], N_std[3];
                Eigen::MatrixXd grad_N(2, 4);
                for (int i = 0; i < 4; ++i) { grad_N(0, i) = dN_xy[0][i]; grad_N(1, i) = dN_xy[1][i]; }
                process_gp(N, grad_N, dV, gp);
            }
        }
    }

    for (size_t i = 0; i < n_nodes; ++i) {
        if (M_lumped(i) > 1e-12) {
            U_nodal(i) /= M_lumped(i);
            W_nodal(i) /= M_lumped(i);
            chi_nodal(i) /= M_lumped(i);
            elec_nodal(i) /= M_lumped(i);
            surf_nodal(i) /= M_lumped(i);
        }
    }
}