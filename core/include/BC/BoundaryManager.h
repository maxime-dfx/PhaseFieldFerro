#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <cmath>
#include "Mesh/Mesh.h"
#include "IO/Datafile.h"

enum class BCType { NEUMANN, DIRICHLET };

struct NodeBC {
    BCType type = BCType::NEUMANN;
    double value = 0.0;
};

using BCProfile = std::function<double(double, double, double)>;

struct BoundaryShape {
    virtual ~BoundaryShape() = default;
    virtual bool contains(double x, double y) const = 0;
};

struct EdgeShape : public BoundaryShape {
    std::string edge;
    double Lx, Ly;
    double tol = 1e-6;

    EdgeShape(const std::string& e, double l_x, double l_y) : edge(e), Lx(l_x), Ly(l_y) {}

    bool contains(double x, double y) const override {
        if (edge == "left")   return x <= tol;
        if (edge == "right")  return x >= Lx - tol;
        if (edge == "bottom") return y <= tol;
        if (edge == "top")    return y >= Ly - tol;
        return false;
    }
};

struct PointShape : public BoundaryShape {
    double px, py;
    double tol = 1e-6;

    PointShape(double x, double y) : px(x), py(y) {}

    bool contains(double x, double y) const override {
        return std::abs(x - px) <= tol && std::abs(y - py) <= tol;
    }
};

struct CircleShape : public BoundaryShape {
    double cx, cy, radius;

    CircleShape(double x, double y, double r) : cx(x), cy(y), radius(r) {}

    bool contains(double x, double y) const override {
        return (x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius * radius;
    }
};

struct RectShape : public BoundaryShape {
    double xmin, xmax, ymin, ymax;

    RectShape(double x0, double x1, double y0, double y1) : xmin(x0), xmax(x1), ymin(y0), ymax(y1) {}

    bool contains(double x, double y) const override {
        return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
    }
};

struct BCRule {
    std::shared_ptr<BoundaryShape> shape;
    BCType type;
    BCProfile profile;
    
    // NOUVEAU : Cache des indices des noeuds impactés par cette règle
    std::vector<int> target_nodes; 
};

class BoundaryManager {
public:
    BoundaryManager(const Mesh& mesh, const std::vector<BoundaryRuleConfig>& boundary_rules);

    const Mesh& get_mesh() const { return m_mesh; }

    void initialize_all_boundaries();
    void update_time(double t);
    void clear_all_rules();

    const std::vector<NodeBC>& get_phi_bcs() const { return m_bc_phi; }
    const std::vector<NodeBC>& get_ux_bcs()  const { return m_bc_ux; }
    const std::vector<NodeBC>& get_uy_bcs()  const { return m_bc_uy; }
    const std::vector<NodeBC>& get_px_bcs()  const { return m_bc_px; }
    const std::vector<NodeBC>& get_py_bcs()  const { return m_bc_py; }

    void add_rule_phi(std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
    void add_rule_ux (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
    void add_rule_uy (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
    void add_rule_px (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
    void add_rule_py (std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);

    void add_rule_phi_constant(std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start);
    void add_rule_ux_constant (std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start);
    void add_rule_uy_constant (std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start);
    void add_rule_px_constant (std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start);
    void add_rule_py_constant (std::shared_ptr<BoundaryShape> shape, BCType type, double value, double t_start);

private:
    const Mesh& m_mesh;
    const std::vector<BoundaryRuleConfig>& m_boundary_rules;

    std::vector<NodeBC> m_bc_phi, m_bc_ux, m_bc_uy, m_bc_px, m_bc_py;
    std::vector<BCRule> m_rules_phi, m_rules_ux, m_rules_uy, m_rules_px, m_rules_py;

    // NOUVEAU : Méthode utilitaire pour créer et cacher la règle
    BCRule create_cached_rule(std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
};