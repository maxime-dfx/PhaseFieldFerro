#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "Utils/toml.hpp"
#include "Utils/Logger.h"
#include "Core/Types.h"

// ------------------------------------------------------------
//  ENUMÉRATIONS 
// ------------------------------------------------------------

struct BoundaryRuleConfig {
    std::string field, bc_type, shape, edge_name, profile;
    double cx = 0.0, cy = 0.0, radius = 0.0;    
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;  
    double px = 0.0, py = 0.0, val = 0.0;       
    double val_start = 0.0, val_end = 0.0, t_end = 1.0;
    double amplitude = 0.0, frequency = 0.0, offset = 0.0;
    double max_val = 0.0, y_center = 0.0, width = 0.0;
    double P0 = 0.0, x0 = 0.0, epsilon = 1e-3, k = 0.0, omega = 0.0;
};

// ------------------------------------------------------------
//  STRUCTURES DE CONFIGURATION (Nouveau !)
// ------------------------------------------------------------
struct SimulationConfig {
    std::string output_dir;
    std::string mesh_create_file;
    int total_time;
    double dt;
    int save_frequency;
    double tol_ferro;
    double tol_vfield;
    bool debug_enabled;
};

struct MeshConfig {
    bool calcul_mesh;
    std::string get_mesh_file;
    double Lx, Ly;
    int nx, ny;
    double dx, dy;
    ElementType element_type;
};

struct ChronoConfig {
    bool mesh, polarization, mechanics, electrostatics, fracture, run;
};

struct PhysicsToggle {
    bool polarization, electrostatics, mechanics, fracture;
};

struct PhysicsInitConfig {
    InitializationType type;
    double val_x_0; // Sert pour Px_0, ux_0, Ex_0
    double val_y_0; // Sert pour Py_0, uy_0, Ey_0
};

struct MaterialConfig {
    double a0, b1, b2, b3, c1, c2, c3, eps0, t, P0, mu_p, mu_v;
    double xi, c0, alpha_1, alpha_11, alpha_111, alpha_1111;
    double alpha_12, alpha_112, alpha_1112, alpha_1122;
    double eta_k, Gc, kappa;
};

// ------------------------------------------------------------
//  CLASSE DATAFILE REFACTORISÉE
// ------------------------------------------------------------
class Datafile {
public:
    // Les sous-structures deviennent accessibles directement
    SimulationConfig simulation;
    MeshConfig mesh;
    ChronoConfig chrono;
    PhysicsToggle physics_toggle;
    PhysicsInitConfig polarization;
    PhysicsInitConfig mechanics;
    PhysicsInitConfig electrostatics;
    CrackBCType fracture_mode;
    MaterialConfig material;
    std::vector<BoundaryRuleConfig> boundary_rules;

    Datafile(const std::string& filename);

private:
    // Méthodes privées pour ranger la logique de parsing
    void parse_simulation(const toml::table& config);
    void parse_mesh(const toml::table& config);
    void parse_chrono(const toml::table& config);
    void parse_physics(const toml::table& config);
    void parse_material(const toml::table& config);
    void parse_boundaries(const toml::table& config);

    // Utilitaire interne
    template <typename T>
    T get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value);
};