#include "Physics/Fracture.h"
#include "Physics/FractureAssembler.h"
#include "Solvers/LinearSolver.h"
#include "Utils/Logger.h"
#include <omp.h>
#include <algorithm>

Fracture::Fracture(const Datafile& config, const Mesh& mesh) : config(config), mesh(mesh) {
    size_t n_nodes = mesh.get_num_nodes();
    
    // Bonne pratique : allouer la mémoire proprement avant l'initialisation
    v_prev_iter.setZero(n_nodes);
    v_current.setZero(n_nodes); 
    v_n.setZero(n_nodes);
    
    appliquer_conditions_initiales();
}

void Fracture::appliquer_conditions_initiales() {
    size_t n_nodes = mesh.get_num_nodes();
    v_prev_iter.setOnes(n_nodes);
    v_current.setOnes(n_nodes);
    v_n.setOnes(n_nodes);

    if (!config.fracture.enable_precrack) {
        return;
    }

    const auto& precrack = config.fracture.precrack;

    for (size_t i = 0; i < n_nodes; ++i) {
        auto coord = mesh.get_node_coords(i); 
        double x = coord[0], y = coord[1];

        double d = calculer_distance_signee_precrack(x, y, precrack);

        double v_val;
        if (precrack.smooth) {

            v_val = 0.5 * (1.0 + std::tanh(d / precrack.smoothing_length));
        } else {
            v_val = (d < 0.0) ? 0.0 : 1.0;
        }

        v_current(i) = v_val;
        v_n(i)       = v_val;
        v_prev_iter(i) = v_val;
    }
}

void Fracture::update_precrack_geometry(double time) {
    if (!config.fracture.enable_precrack) return;

    PrecrackConfig pc = config.fracture.precrack;   // copie modifiable
    // La longueur croit lineairement avec le temps, comme vos BC mecaniques (time_ramp)
    double t_end = config.fracture.precrack.growth_t_end; // a definir dans le TOML
    double length_start = config.fracture.precrack.growth_length_start;
    double length_end = config.fracture.precrack.growth_length_end;

    double frac = std::min(1.0, std::max(0.0, time / t_end));
    pc.length = length_start + frac * (length_end - length_start);

    size_t n_nodes = mesh.get_num_nodes();
    for (size_t i = 0; i < n_nodes; ++i) {
        auto coord = mesh.get_node_coords(i);
        double x = coord[0], y = coord[1];
        double d = calculer_distance_signee_precrack(x, y, pc);

        double v_val;
        if (pc.smooth) {
            v_val = 0.5 * (1.0 + std::tanh(d / pc.smoothing_length));
        } else {
            v_val = (d < 0.0) ? 0.0 : 1.0;
        }

        v_current(i) = v_val;
        v_n(i) = v_val;
    }
}

double Fracture::calculer_distance_signee_precrack(double x, double y, const PrecrackConfig& pc) const {
    if (pc.shape == PrecrackShape::SEGMENT) {
        double dx = x - pc.x0;
        double dy = y - pc.y0;

        if (dx < 0.0) {
            return std::sqrt(dx * dx + dy * dy) - pc.half_width;
        } else if (dx > pc.length) {
            double dxt = dx - pc.length;
            return std::sqrt(dxt * dxt + dy * dy) - pc.half_width;
        } else {
            return std::abs(dy) - pc.half_width;
        }
    }
    else if (pc.shape == PrecrackShape::RECT) {
        double dx = std::max({pc.xmin - x, x - pc.xmax, 0.0});
        double dy = std::max({pc.ymin - y, y - pc.ymax, 0.0});
        if (dx == 0.0 && dy == 0.0) {
            double inside = std::min({x - pc.xmin, pc.xmax - x, y - pc.ymin, pc.ymax - y});
            return -inside;
        }
        return std::sqrt(dx * dx + dy * dy);
    }

    return 1e6;
}

void Fracture::update_v(double dt_relax, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
    size_t n_nodes = mesh.get_num_nodes();
    Eigen::SparseMatrix<double> K_global(n_nodes, n_nodes);
    Eigen::VectorXd F_global = Eigen::VectorXd::Zero(n_nodes);

    FractureAssembler::assemble_system(dt_relax, v_current, mesh, polarization, mechanics, electrostatics, math, config, K_global, F_global);

    Eigen::VectorXd v_new = LinearSolver::solve_with_guess(K_global, F_global, v_current, config.simulation.debug_enabled);

    if (v_new.allFinite()) {
        enforce_physical_bounds(v_new);
    } else {
        Logger::error("[Fracture] Echec resolution ou solution non finie");
    }
}

void Fracture::enforce_physical_bounds(const Eigen::VectorXd& v_new) {
    size_t n_nodes = mesh.get_num_nodes();
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_nodes; ++i) {
        // v doit rester entre 0 (cassé) et v_n (condition de non-guérison locale)
        v_current(i) = std::max(0.0, std::min(v_new(i), v_n(i))); 
    }
}

double Fracture::calculate_error() const { return (v_current - v_prev_iter).cwiseAbs().maxCoeff(); }
void Fracture::save_previous_iteration() { v_prev_iter = v_current; }
void Fracture::save_previous_state() { v_backup = v_current; }
void Fracture::restore_previous_state() { v_current = v_backup; }
void Fracture::update_history() { v_n = v_current; }