#include "BC/BoundaryManager.h"
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Profiles {
    inline BCProfile Constant(double val, double t_start) {
        return [=](double, double, double t) {
            if (t <= t_start) return 0.0;
            return val;
        };
    }
    
    inline BCProfile TimeRamp(double val_start, double val_end, double t_start, double t_end) {
        return [=](double, double, double t) {
            if (t <= t_start) return val_start;
            if (t >= t_end) return val_end;
            return val_start + (val_end - val_start) * ((t - t_start) / (t_end - t_start));
        };
    }
    
    inline BCProfile TimeSine(double amplitude, double frequency, double offset) {
        return [=](double, double, double t) {
            return offset + amplitude * std::sin(2.0 * M_PI * frequency * t);
        };
    }
    
    inline BCProfile SpatialTanhX(double P0, double x0, double epsilon) {
        return [=](double x, double, double) {
            return P0 * std::tanh((x - x0) / epsilon);
        };
    }
    
    inline BCProfile SpatialParabolaY(double max_val, double y_center, double width) {
        return [=](double, double y, double) {
            double dist = std::abs(y - y_center);
            if (dist >= width) return 0.0;
            return max_val * (1.0 - (dist * dist) / (width * width));
        };
    }
    
    inline BCProfile TravelingWave(double amplitude, double k, double omega) {
        return [=](double x, double, double t) {
            return amplitude * std::sin(k * x - omega * t);
        };
    }
}

BoundaryManager::BoundaryManager(const Mesh& mesh, const std::vector<BoundaryRuleConfig>& boundary_rules) 
    : m_mesh(mesh), m_boundary_rules(boundary_rules) 
{
    int n_nodes = m_mesh.get_num_nodes();
    m_bc_phi.resize(n_nodes); m_bc_ux.resize(n_nodes); m_bc_uy.resize(n_nodes);
    m_bc_px.resize(n_nodes); m_bc_py.resize(n_nodes);
}

// =========================================================================
// NOUVEAU : Création de la règle et pré-calcul (cache) des noeuds cibles
// =========================================================================
BCRule BoundaryManager::create_cached_rule(std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) {
    BCRule rule{shape, type, profile, {}};
    
    const auto& all_nodes = m_mesh.get_nodes();
    int n_nodes = static_cast<int>(all_nodes.size());
    
    // On n'évalue la géométrie qu'une seule fois à l'initialisation !
    for (int idx = 0; idx < n_nodes; ++idx) {
        if (shape->contains(all_nodes[idx].x, all_nodes[idx].y)) {
            rule.target_nodes.push_back(idx);
        }
    }
    return rule;
}

// --- Mise à jour des fonctions d'ajout pour utiliser le cache ---
void BoundaryManager::add_rule_phi(std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { 
    m_rules_phi.push_back(create_cached_rule(shape, type, profile)); 
}
void BoundaryManager::add_rule_ux (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { 
    m_rules_ux.push_back(create_cached_rule(shape, type, profile)); 
}
void BoundaryManager::add_rule_uy (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { 
    m_rules_uy.push_back(create_cached_rule(shape, type, profile)); 
}
void BoundaryManager::add_rule_px (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { 
    m_rules_px.push_back(create_cached_rule(shape, type, profile)); 
}
void BoundaryManager::add_rule_py (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { 
    m_rules_py.push_back(create_cached_rule(shape, type, profile)); 
}

void BoundaryManager::add_rule_phi_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start) { 
    add_rule_phi(shape, type, Profiles::Constant(value, t_start)); 
}
void BoundaryManager::add_rule_ux_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start)  { 
    add_rule_ux(shape, type, Profiles::Constant(value, t_start)); 
}
void BoundaryManager::add_rule_uy_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start)  { 
    add_rule_uy(shape, type, Profiles::Constant(value, t_start)); 
}
void BoundaryManager::add_rule_px_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start)  { 
    add_rule_px(shape, type, Profiles::Constant(value, t_start)); 
}
void BoundaryManager::add_rule_py_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start)  { 
    add_rule_py(shape, type, Profiles::Constant(value, t_start)); 
}

// =========================================================================
// OPTIMISATION EXTRÊME DE LA BOUCLE TEMPORELLE
// =========================================================================
void BoundaryManager::update_time(double t) {
    const auto& all_nodes = m_mesh.get_nodes();
    
    NodeBC def = {BCType::NEUMANN, 0.0};
    std::fill(m_bc_phi.begin(), m_bc_phi.end(), def);
    std::fill(m_bc_ux.begin(), m_bc_ux.end(), def);
    std::fill(m_bc_uy.begin(), m_bc_uy.end(), def);
    std::fill(m_bc_px.begin(), m_bc_px.end(), def);
    std::fill(m_bc_py.begin(), m_bc_py.end(), def);

    auto apply_rules = [&](const std::vector<BCRule>& rules, std::vector<NodeBC>& nodal_bcs) {
        if (rules.empty()) return;
        for (const auto& rule : rules) {
            // BOUCLE CACHÉE : On ne parcourt QUE les noeuds concernés !
            for (int idx : rule.target_nodes) {
                double x = all_nodes[idx].x;
                double y = all_nodes[idx].y;
                nodal_bcs[idx] = {rule.type, rule.profile(x, y, t)};
            }
        }
    };

    apply_rules(m_rules_phi, m_bc_phi);
    apply_rules(m_rules_ux, m_bc_ux);
    apply_rules(m_rules_uy, m_bc_uy);
    apply_rules(m_rules_px, m_bc_px);
    apply_rules(m_rules_py, m_bc_py);
    
    // int count_phi = 0, count_px = 0, count_ux = 0;
    // for (const auto& bc : m_bc_phi) { if (bc.type == BCType::DIRICHLET) count_phi++; }
    // for (const auto& bc : m_bc_px) { if (bc.type == BCType::DIRICHLET) count_px++; }
    // for (const auto& bc : m_bc_ux) { if (bc.type == BCType::DIRICHLET) count_ux++; }
    
    // Logger::debug("[DEBUG CL] Application a t=%f -> Noeuds Dirichlet : phi=%d, px=%d, ux=%d", t, count_phi, count_px, count_ux, m_config.simulation.debug_enabled);
}

void BoundaryManager::clear_all_rules() {
    m_rules_phi.clear(); m_rules_ux.clear(); m_rules_uy.clear(); m_rules_px.clear(); m_rules_py.clear();
}

void BoundaryManager::initialize_all_boundaries() {
    auto rules = m_boundary_rules;
    double Lx = m_mesh.get_Lx();
    double Ly = m_mesh.get_Ly();

    for (const auto& rule : m_boundary_rules) {
        std::shared_ptr<BoundaryShape> shape;
        if (rule.shape == "edge")        shape = std::make_shared<EdgeShape>(rule.edge_name, Lx, Ly);
        else if (rule.shape == "point")  shape = std::make_shared<PointShape>(rule.px, rule.py);
        else if (rule.shape == "circle") shape = std::make_shared<CircleShape>(rule.cx, rule.cy, rule.radius);
        else if (rule.shape == "rect")   shape = std::make_shared<RectShape>(rule.xmin, rule.xmax, rule.ymin, rule.ymax);
        else continue;

        BCProfile profile;
        if (rule.profile == "constant")              profile = Profiles::Constant(rule.val, rule.t_start);
        else if (rule.profile == "time_ramp")        profile = Profiles::TimeRamp(rule.val_start, rule.val_end, rule.t_start, rule.t_end);
        else if (rule.profile == "time_sine")        profile = Profiles::TimeSine(rule.amplitude, rule.frequency, rule.offset);
        else if (rule.profile == "spatial_tanh")     profile = Profiles::SpatialTanhX(rule.P0, rule.x0, rule.epsilon);
        else if (rule.profile == "spatial_parabola") profile = Profiles::SpatialParabolaY(rule.max_val, rule.y_center, rule.width);
        else if (rule.profile == "traveling_wave")   profile = Profiles::TravelingWave(rule.amplitude, rule.k, rule.omega);
        else profile = Profiles::Constant(rule.val, rule.t_start);

        BCType type = (rule.bc_type == "NEUMANN") ? BCType::NEUMANN : BCType::DIRICHLET;

        // Logger::debug("[DEBUG CL] Lecture TOML - field: '%s' | shape: %s | bc_type: %s", rule.field.c_str(), rule.shape.c_str(), rule.bc_type.c_str(), m_config.simulation.debug_enabled);

        if (rule.field == "phi")      add_rule_phi(shape, type, profile);
        else if (rule.field == "ux")  add_rule_ux(shape, type, profile);
        else if (rule.field == "uy")  add_rule_uy(shape, type, profile);
        else if (rule.field == "px")  add_rule_px(shape, type, profile);
        else if (rule.field == "py")  add_rule_py(shape, type, profile);
    }
}