#pragma once

#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <vector>

#include "IO/Datafile.h"
#include "Core/Mesh.h"

class Polarization;
class Mechanics;
class Fracture;
class Electrostatics;
class Math;

struct EnergyRecord {
    double load_step   = 0.0;
    double time         = 0.0;
    double U_total      = 0.0;
    double W_total       = 0.0;
    double chi_total     = 0.0;
    double electric_total = 0.0;
    double bulk_enthalpy  = 0.0;
    double surface_energy = 0.0;
    double total_energy = 0.0;
    double v_min = 1.0;
    double v_max = 1.0;
    double crack_length_proxy = 0.0;
    double phi_min = 0.0, phi_max = 0.0;
    double Px_min = 0.0, Px_max = 0.0, Py_min = 0.0, Py_max = 0.0;
    double ux_max_abs = 0.0, uy_max_abs = 0.0;
};

class Diagnostics {
public:
    Diagnostics(const Datafile& config, const Mesh& mesh);

    EnergyRecord record(double load_step, double time,
                         const Polarization& polarization,
                         const Mechanics& mechanics,
                         const Fracture& fracture,
                         const Electrostatics& electrostatics,
                         const Math& math);

    void write_csv(const std::string& path) const;
    void append_csv(const std::string& path);
    const std::vector<EnergyRecord>& get_history() const { return history; }

    // =====================================================================
    // NOUVELLES FONCTIONS PUBLIQUES POUR L'EXPORT VTK
    // =====================================================================
    void compute_nodal_energies(const Polarization& polarization,
                                const Mechanics& mechanics,
                                const Fracture& fracture,
                                const Electrostatics& electrostatics,
                                const Math& math);

    const Eigen::VectorXd& get_U_nodal() const { return U_nodal; }
    const Eigen::VectorXd& get_W_nodal() const { return W_nodal; }
    const Eigen::VectorXd& get_chi_nodal() const { return chi_nodal; }
    const Eigen::VectorXd& get_elec_nodal() const { return elec_nodal; }
    const Eigen::VectorXd& get_surf_nodal() const { return surf_nodal; }

private:
    const Datafile& config;
    const Mesh& mesh;
    std::vector<EnergyRecord> history;
    bool append_header_written = false;

    // =====================================================================
    // VECTEURS DE STOCKAGE DES ÉNERGIES (Privés)
    // =====================================================================
    Eigen::VectorXd U_nodal;
    Eigen::VectorXd W_nodal;
    Eigen::VectorXd chi_nodal;
    Eigen::VectorXd elec_nodal;
    Eigen::VectorXd surf_nodal;

    void integrate_bulk_terms(const Polarization& polarization,
                               const Mechanics& mechanics,
                               const Fracture& fracture,
                               const Electrostatics& electrostatics,
                               const Math& math,
                               EnergyRecord& rec) const;

    void integrate_surface_term(const Fracture& fracture, EnergyRecord& rec) const;

    void compute_nodal_extrema(const Polarization& polarization,
                                const Mechanics& mechanics,
                                const Fracture& fracture,
                                const Electrostatics& electrostatics,
                                EnergyRecord& rec) const;
};