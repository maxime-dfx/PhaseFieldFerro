#pragma once
#include <string>
#include "Utils/Types.h"

struct BoundaryRuleConfig {
    std::string field, bc_type, shape, edge_name, profile;
    double cx = 0.0, cy = 0.0, radius = 0.0;    
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;  
    double px = 0.0, py = 0.0, val = 0.0;       
    double val_start = 0.0, val_end = 0.0, t_start = 0.0, t_end = 1.0;
    double amplitude = 0.0, frequency = 0.0, offset = 0.0;
    double max_val = 0.0, y_center = 0.0, width = 0.0;
    double P0 = 0.0, x0 = 0.0, epsilon = 1e-3, k = 0.0, omega = 0.0;
};

struct SimulationConfig {
    std::string output_dir;
    std::string mesh_create_file;
    int total_time;
    double dt;
    double dt_relax;
    int save_frequency;
    double tol_ferro;
    double tol_vfield;
    int min_iter;
    int max_iter;
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

struct PrecrackConfig {
    PrecrackShape shape = PrecrackShape::NONE;
    double x0 = 0.0;
    double y0 = 50.0;
    double length = 5.0;
    double half_width = 0.5;
    double xmin = 0.0, xmax = 5.0, ymin = 48.0, ymax = 52.0;
    bool smooth = true;
    double smoothing_length = 1.0;
    bool growth_enable = false;
    double growth_length_start = 5.0;
    double growth_length_end   = 5.0;
    double growth_t_end        = 3.0;
};

struct FractureConfig {
    CrackBCType mode = CrackBCType::PERMEABLE;
    bool enable_precrack = false;
    PrecrackConfig precrack;
};