#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <unordered_map>
#include <cmath>
#include "Mesh/include/Mesh.h"
#include "IO/include/Datafile.h"

enum class BCType { NEUMANN, DIRICHLET };

struct NodeBC {
    BCType type = BCType::NEUMANN;
    double value = 0.0;
};

using BCProfile = std::function<double(double, double, double)>;

// ============================= Shapes (inchangés) =============================

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
        return (x - cx) * (x - cx) + (y - cy) * (y - cy) <= radius * radius;
    }
};

struct RectShape : public BoundaryShape {
    double xmin, xmax, ymin, ymax;

    RectShape(double x0, double x1, double y0, double y1) : xmin(x0), xmax(x1), ymin(y0), ymax(y1) {}

    bool contains(double x, double y) const override {
        return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
    }
};

// ============================= BCRule (data-driven) =============================
// NOUVEAU : porte désormais le nom du champ physique concerné ("phi", "ux", ...).
// Cela permet de supprimer toute connaissance codée en dur des champs dans le manager.
struct BCRule {
    std::string field;
    std::shared_ptr<BoundaryShape> shape;
    BCType type;
    BCProfile profile;

    // Cache des indices des noeuds impactés par cette règle
    std::vector<int> target_nodes;
};

// ============================= BoundaryManager =============================
// Remplace BoundaryManager. Ne connaît plus aucun champ physique à l'avance (OCP) :
// toutes les BC nodales sont stockées dans une map générique field-name -> BCs.
class BoundaryManager {
public:
    BoundaryManager(const Mesh& mesh, const std::vector<BoundaryRuleConfig>& boundary_rules);

    const Mesh& get_mesh() const { return m_mesh; }

    void initialize_all_boundaries();
    void update_time(double t);
    void clear_all_rules();

    // Ajoute une règle pour un champ arbitraire, identifié par son nom.
    void add_rule(const std::string& field, std::shared_ptr<BoundaryShape> shape, BCType type, BCProfile profile);
    void add_rule_constant(const std::string& field, std::shared_ptr<BoundaryShape> shape, BCType type,
                            double value, double t_start);

    // Retourne les BC nodales d'un champ. Si le champ n'a jamais été enregistré,
    // retourne un vecteur par défaut (Neumann 0.0, taille = nombre de noeuds).
    const std::vector<NodeBC>& get_bcs(const std::string& field_name) const;

private:
    const Mesh& m_mesh;
    const std::vector<BoundaryRuleConfig>& m_boundary_rules;

    // Clé = nom du champ ("phi", "ux", "uy", "px", "py", "v", ...)
    std::unordered_map<std::string, std::vector<NodeBC>> m_nodal_bcs;
    std::unordered_map<std::string, std::vector<BCRule>> m_rules;

    // Vecteur par défaut retourné pour un champ inconnu (Neumann homogène).
    std::vector<NodeBC> m_default_bcs;

    BCRule create_cached_rule(const std::string& field, std::shared_ptr<BoundaryShape> shape,
                               BCType type, BCProfile profile);
};
