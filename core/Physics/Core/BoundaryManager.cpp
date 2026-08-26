#include "Physics/Core/BoundaryManager.h"
#include "Utils/include/Profiling.h"
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
    // Vecteur par défaut (Neumann 0.0) retourné pour tout champ non-enregistré.
    m_default_bcs.assign(n_nodes, NodeBC{BCType::NEUMANN, 0.0});

    // Pré-déclare les champs présents dans la config pour qu'ils soient
    // correctement dimensionnés dès l'appel à update_time().
    for (const auto& rule : m_boundary_rules) {
        if (m_nodal_bcs.find(rule.field) == m_nodal_bcs.end()) {
            m_nodal_bcs[rule.field] = std::vector<NodeBC>(n_nodes);
        }
    }
}

BCRule BoundaryManager::create_cached_rule(const std::string& field, std::shared_ptr<BoundaryShape> shape,
                                            BCType type, BCProfile profile) {
    BCRule rule{field, shape, type, profile, {}};

    const auto& all_nodes = m_mesh.get_nodes();
    int n_nodes = static_cast<int>(all_nodes.size());

    for (int idx = 0; idx < n_nodes; ++idx) {
        if (shape->contains(all_nodes[idx].x, all_nodes[idx].y)) {
            rule.target_nodes.push_back(idx);
        }
    }
    return rule;
}

void BoundaryManager::add_rule(const std::string& field, std::shared_ptr<BoundaryShape> shape,
                                BCType type, BCProfile profile) {
    if (m_nodal_bcs.find(field) == m_nodal_bcs.end()) {
        m_nodal_bcs[field] = std::vector<NodeBC>(m_mesh.get_num_nodes());
    }
    m_rules[field].push_back(create_cached_rule(field, shape, type, profile));
}

void BoundaryManager::add_rule_constant(const std::string& field, std::shared_ptr<BoundaryShape> shape,
                                         BCType type, double value, double t_start) {
    add_rule(field, shape, type, Profiles::Constant(value, t_start));
}

const std::vector<NodeBC>& BoundaryManager::get_bcs(const std::string& field_name) const {
    auto it = m_nodal_bcs.find(field_name);
    if (it != m_nodal_bcs.end()) {
        return it->second;
    }
    // Champ jamais enregistré : comportement par défaut Neumann homogène.
    return m_default_bcs;
}

void BoundaryManager::update_time(double t) {
    PROFILE_ZONE_NC("BoundaryManager::update_time", PROFILE_COLOR_SEQUENTIAL);
    const auto& all_nodes = m_mesh.get_nodes();

    const NodeBC def = {BCType::NEUMANN, 0.0};
    for (auto& [field, bcs] : m_nodal_bcs) {
        std::fill(bcs.begin(), bcs.end(), def);
    }

    for (auto& [field, rules] : m_rules) {
        auto& nodal_bcs = m_nodal_bcs[field];
        for (const auto& rule : rules) {
            for (int idx : rule.target_nodes) {
                double x = all_nodes[idx].x;
                double y = all_nodes[idx].y;
                nodal_bcs[idx] = {rule.type, rule.profile(x, y, t)};
            }
        }
    }
}

void BoundaryManager::clear_all_rules() {
    m_rules.clear();
}

void BoundaryManager::initialize_all_boundaries() {
    PROFILE_ZONE_NC("BoundaryManager::initialize_all_boundaries", PROFILE_COLOR_SEQUENTIAL);
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

        // Data-driven : n'importe quel nom de champ est accepté sans modifier cette classe (OCP).
        add_rule(rule.field, shape, type, profile);
    }
    update_time(0.0);
}
