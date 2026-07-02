#include "Core/BoundaryManager.h"
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Profiles {
    inline BCProfile Constant(double val) {
        return [val](double, double, double) { return val; };
    }

    inline BCProfile TimeRamp(double val_start, double val_end, double t_end) {
        return [=](double, double, double t) {
            if (t >= t_end) return val_end;
            return val_start + (val_end - val_start) * (t / t_end);
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

BoundaryManager::BoundaryManager(const Mesh& mesh, const Datafile& config)
    : m_mesh(mesh), m_config(config) 
{
    int n_nodes = m_mesh.get_nx() * m_mesh.get_ny();
    m_bc_phi.resize(n_nodes); m_bc_ux.resize(n_nodes); m_bc_uy.resize(n_nodes);
    m_bc_px.resize(n_nodes); m_bc_py.resize(n_nodes);
    Logger::debug("[DEBUG][BCManager] phi_bcs.size()=" + std::to_string(get_phi_bcs().size()) + "\n", m_config.debug_enabled());
    size_t n_d = 0;
    for (auto& bc : get_phi_bcs()) if (bc.type == BCType::DIRICHLET) n_d++;
    Logger::debug("[DEBUG][BCManager] phi dirichlet count=" + std::to_string(n_d) + "\n", m_config.debug_enabled());
}

void BoundaryManager::add_rule_phi(std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { m_rules_phi.push_back({shape, type, profile}); }
void BoundaryManager::add_rule_ux (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { m_rules_ux.push_back({shape, type, profile}); }
void BoundaryManager::add_rule_uy (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { m_rules_uy.push_back({shape, type, profile}); }
void BoundaryManager::add_rule_px (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { m_rules_px.push_back({shape, type, profile}); }
void BoundaryManager::add_rule_py (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile) { m_rules_py.push_back({shape, type, profile}); }

void BoundaryManager::add_rule_phi_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value) { add_rule_phi(shape, type, Profiles::Constant(value)); }
void BoundaryManager::add_rule_ux_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value)  { add_rule_ux(shape, type, Profiles::Constant(value)); }
void BoundaryManager::add_rule_uy_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value)  { add_rule_uy(shape, type, Profiles::Constant(value)); }
void BoundaryManager::add_rule_px_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value)  { add_rule_px(shape, type, Profiles::Constant(value)); }
void BoundaryManager::add_rule_py_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value)  { add_rule_py(shape, type, Profiles::Constant(value)); }

void BoundaryManager::update_time(double t) {
    int nx = m_mesh.get_nx(), ny = m_mesh.get_ny();
    double dx = m_mesh.get_dx(), dy = m_mesh.get_dy();
    Logger::debug("[BCManager] nx=" + std::to_string(nx) + " ny=" + std::to_string(ny) + " dx=" + std::to_string(dx) + " dy=" + std::to_string(dy) + "\n", m_config.debug_enabled());

    NodeBC def = {BCType::NEUMANN, 0.0};
    std::fill(m_bc_phi.begin(), m_bc_phi.end(), def);
    std::fill(m_bc_ux.begin(), m_bc_ux.end(), def);
    std::fill(m_bc_uy.begin(), m_bc_uy.end(), def);
    std::fill(m_bc_px.begin(), m_bc_px.end(), def);
    std::fill(m_bc_py.begin(), m_bc_py.end(), def);

    auto apply_rules = [&](const std::vector<BCRule>& rules, std::vector<NodeBC>& nodal_bcs) {
        if (rules.empty()) return;
        for (const auto& rule : rules) {
            for (int j = 0; j < ny; ++j) {
                double y = j * dy;
                for (int i = 0; i < nx; ++i) {
                    double x = i * dx;
                    if (rule.shape->contains(x, y)) {
                        nodal_bcs[j * nx + i] = {rule.type, rule.profile(x, y, t)};
                    }
                }
            }
        }
    };

    apply_rules(m_rules_phi, m_bc_phi); 
    apply_rules(m_rules_ux, m_bc_ux); apply_rules(m_rules_uy, m_bc_uy); 
    apply_rules(m_rules_px, m_bc_px); apply_rules(m_rules_py, m_bc_py);
}

void BoundaryManager::clear_all_rules() {
    m_rules_phi.clear(); m_rules_ux.clear(); m_rules_uy.clear(); m_rules_px.clear(); m_rules_py.clear();
}

void BoundaryManager::initialize_all_boundaries() {
    auto rules = m_config.get_boundary_rules();
    double Lx = m_mesh.get_Lx();
    double Ly = m_mesh.get_Ly();

    for (const auto& r : rules) {
        std::shared_ptr<BoundaryShape> shape;
        if (r.shape == "edge")        shape = std::make_shared<EdgeShape>(r.edge_name, Lx, Ly);
        else if (r.shape == "point")  shape = std::make_shared<PointShape>(r.px, r.py);
        else if (r.shape == "circle") shape = std::make_shared<CircleShape>(r.cx, r.cy, r.radius);
        else if (r.shape == "rect")   shape = std::make_shared<RectShape>(r.xmin, r.xmax, r.ymin, r.ymax);
        else continue; 

        BCProfile profile;
        if (r.profile == "constant")              profile = Profiles::Constant(r.val);
        else if (r.profile == "time_ramp")        profile = Profiles::TimeRamp(r.val_start, r.val_end, r.t_end);
        else if (r.profile == "time_sine")        profile = Profiles::TimeSine(r.amplitude, r.frequency, r.offset);
        else if (r.profile == "spatial_tanh")     profile = Profiles::SpatialTanhX(r.P0, r.x0, r.epsilon);
        else if (r.profile == "spatial_parabola") profile = Profiles::SpatialParabolaY(r.max_val, r.y_center, r.width);
        else if (r.profile == "traveling_wave")   profile = Profiles::TravelingWave(r.amplitude, r.k, r.omega);
        else profile = Profiles::Constant(r.val); 

        BCType type = (r.bc_type == "NEUMANN") ? BCType::NEUMANN : BCType::DIRICHLET;
        
        if (r.field == "phi")      add_rule_phi(shape, type, profile);
        else if (r.field == "ux")  add_rule_ux(shape, type, profile);
        else if (r.field == "uy")  add_rule_uy(shape, type, profile);
        else if (r.field == "px")  add_rule_px(shape, type, profile);
        else if (r.field == "py")  add_rule_py(shape, type, profile);
    }
}