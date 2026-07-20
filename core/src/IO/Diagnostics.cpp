#include "IO/Diagnostics.h"
#include "Physics/Polarization.h"
#include "Physics/Mechanics.h"
#include "Physics/Fracture.h"
#include "Physics/Electrostatics.h"
#include "Physics/Math.h"
#include "FEM/ElementIntegrator.h"
#include <iostream>
#include <iomanip>

Diagnostics::Diagnostics(const Datafile& config, const Mesh& mesh)
    : config(config), mesh(mesh) {
    U_nodal = Eigen::VectorXd::Zero(mesh.get_num_nodes());
    W_nodal = Eigen::VectorXd::Zero(mesh.get_num_nodes());
    chi_nodal = Eigen::VectorXd::Zero(mesh.get_num_nodes());
    elec_nodal = Eigen::VectorXd::Zero(mesh.get_num_nodes());
    surf_nodal = Eigen::VectorXd::Zero(mesh.get_num_nodes());
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

    rec.total_energy = rec.bulk_enthalpy + rec.surface_energy;
    
    compute_nodal_extrema(polarization, mechanics, fracture, electrostatics, rec);
    
    history.push_back(rec);
    return rec;
}

void Diagnostics::integrate_bulk_terms(const Polarization& polarization, const Mechanics& mechanics, 
                                       const Fracture& fracture, const Electrostatics& electrostatics, 
                                       const Math& math, EnergyRecord& rec) const {
    double U_sum = 0.0, W_sum = 0.0, chi_sum = 0.0, elec_sum = 0.0;
    
    Eigen::RowVectorXd N(8);
    Eigen::MatrixXd grad_N(2, 8);
    
    bool is_impermeable = (config.fracture.mode == CrackBCType::IMPERMEABLE);

    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        const Element& elem = mesh.get_elements()[e];
        auto coords = mesh.get_element_coords(e);
        
        ElementIntegrator::integrate(elem, coords, N, grad_N, 
            [&](int n_nodes, double dV, const GaussPoint2D& gp) {
                Eigen::Vector2d P = Eigen::Vector2d::Zero();
                Eigen::Matrix2d grad_P = Eigen::Matrix2d::Zero();
                double v_val = 0.0;

                for (int i = 0; i < n_nodes; ++i) {
                    int node_id = elem.get_node_index(i);
                    double Px = polarization.get_Px()(node_id);
                    double Py = polarization.get_Py()(node_id);
                    
                    P(0) += Px * N(i);
                    P(1) += Py * N(i);
                    grad_P(0,0) += Px * grad_N(0, i);
                    grad_P(0,1) += Px * grad_N(1, i);
                    grad_P(1,0) += Py * grad_N(0, i);
                    grad_P(1,1) += Py * grad_N(1, i);
                    v_val += fracture.get_v()(node_id) * N(i);
                }

                Eigen::Matrix2d strain = mechanics.get_strain_at_gp(elem, gp, coords);
                Eigen::Vector2d E(electrostatics.get_Ex_at_gp(elem, gp), electrostatics.get_Ey_at_gp(elem, gp));

                double deg = (v_val * v_val + math.eta_k);
                
                U_sum   += deg * math.U_energy(grad_P) * dV;
                W_sum   += deg * math.W_energy(P, strain) * dV;
                chi_sum += math.chi_energy(P) * dV;

                if (is_impermeable) {
                    elec_sum += deg * (-P.dot(E) - 0.5 * math.eps0 * E.dot(E)) * dV;
                } else {
                    elec_sum += (-P.dot(E) - 0.5 * math.eps0 * E.dot(E)) * dV;
                }
            });
    }

    rec.U_total = U_sum;
    rec.W_total = W_sum;
    rec.chi_total = chi_sum;
    rec.electric_total = elec_sum;
    rec.bulk_enthalpy = U_sum + W_sum + chi_sum + elec_sum;
}

void Diagnostics::integrate_surface_term(const Fracture& fracture, EnergyRecord& rec) const {
    double Gc = config.material.Gc;
    double kappa = config.material.kappa;
    double surf_sum = 0.0;
    
    Eigen::RowVectorXd N(8);
    Eigen::MatrixXd grad_N(2, 8);

    for (int e = 0; e < mesh.get_num_elements(); ++e) {
        const Element& elem = mesh.get_elements()[e];
        auto coords = mesh.get_element_coords(e);
        
        ElementIntegrator::integrate(elem, coords, N, grad_N, 
            [&](int n_nodes, double dV, const GaussPoint2D& gp) {
                double v_val = 0.0;
                Eigen::Vector2d grad_v = Eigen::Vector2d::Zero();
                for (int i = 0; i < n_nodes; ++i) {
                    int node_id = elem.get_node_index(i);
                    double v = fracture.get_v()(node_id);
                    v_val += v * N(i);
                    grad_v(0) += v * grad_N(0, i);
                    grad_v(1) += v * grad_N(1, i);
                }
                double surf = ((1.0 - v_val)*(1.0 - v_val)) / (4.0 * kappa) + kappa * grad_v.dot(grad_v);
                surf_sum += Gc * surf * dV;
            });
    }
    rec.surface_energy = surf_sum;
}

void Diagnostics::compute_nodal_extrema(const Polarization& polarization, const Mechanics& mechanics, 
                                        const Fracture& fracture, const Electrostatics& electrostatics, 
                                        EnergyRecord& rec) const {
    int n_nodes = mesh.get_num_nodes();
    rec.v_min = 1.0; rec.v_max = 0.0;
    rec.phi_min = 1e10; rec.phi_max = -1e10;

    for (int i = 0; i < n_nodes; ++i) {
        double v = fracture.get_v()(i);
        if (v < rec.v_min) rec.v_min = v;
        if (v > rec.v_max) rec.v_max = v;

        double phi = electrostatics.get_phi()(i);
        if (phi < rec.phi_min) rec.phi_min = phi;
        if (phi > rec.phi_max) rec.phi_max = phi;
    }
}

void Diagnostics::write_csv(const std::string& path) const {
    std::ofstream file(path);
    file << "LoadStep,Time,U,W,Chi,Elec,Bulk,Surface,Total,Vmin,Vmax\n";
    for (const auto& r : history) {
        file << r.load_step << "," << r.time << "," << r.U_total << "," << r.W_total << "," 
             << r.chi_total << "," << r.electric_total << "," << r.bulk_enthalpy << "," 
             << r.surface_energy << "," << r.total_energy << "," << r.v_min << "," << r.v_max << "\n";
    }
}

void Diagnostics::append_csv(const std::string& path) {
    std::ofstream file(path, std::ios::app);
    if (!append_header_written) {
        file << "LoadStep,Time,U,W,Chi,Elec,Bulk,Surface,Total,Vmin,Vmax\n";
        append_header_written = true;
    }
    const auto& r = history.back();
    file << r.load_step << "," << r.time << "," << r.U_total << "," << r.W_total << "," 
         << r.chi_total << "," << r.electric_total << "," << r.bulk_enthalpy << "," 
         << r.surface_energy << "," << r.total_energy << "," << r.v_min << "," << r.v_max << "\n";
}

void Diagnostics::compute_nodal_energies(const Polarization& polarization, const Mechanics& mechanics, 
                                         const Fracture& fracture, const Electrostatics& electrostatics, const Math& math) {
    // Note: Implémentation simplifiée par projection nodale (L2) si nécessaire
    // Cette partie est souvent spécifique au type d'export VTK choisi.
}