#pragma once

#include <vector>
#include <array>
#include <cassert>
#include "Core/Types.h"

struct Node {
    double x, y, z;
    int ref_tag;
};

struct Element {
    // Optimisation HPC : Capacité maximale fixée à 8 (pour le Q8) pour éviter l'allocation dynamique
    std::array<int, 8> node_indices{}; 
    int num_nodes = 0; 
    int ref_tag = 0;

    int get_num_nodes() const { return num_nodes; }
    
    // Accès direct ultra-rapide
    int get_node_index(int i) const {
        assert(i < num_nodes && "Index hors limites pour cet element !");
        return node_indices[i];
    }

    // Itérateurs légers permettant de garder la compatibilité avec les boucles "for-range"
    // Exemple dans l'assembleur : for(int node_id : elem) { ... }
    const int* begin() const { return node_indices.data(); }
    const int* end() const { return node_indices.data() + num_nodes; }
};

class Mesh {
private:
    double m_Lx, m_Ly;
    int m_nx, m_ny;
    ElementType m_element_type;
    
    std::vector<Node> m_nodes;
    std::vector<Element> m_elements;

public:
    Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
         std::vector<Node> nodes, std::vector<Element> elements);

    const std::vector<Node>& get_nodes() const { return m_nodes; }
    const std::vector<Element>& get_elements() const { return m_elements; }
    
    ElementType get_element_type() const { return m_element_type; }
    
    double get_Lx() const { return m_Lx; }
    double get_Ly() const { return m_Ly; }
    
    int get_nx() const { return m_nx; }
    int get_ny() const { return m_ny; }
    
    double get_dx() const { return (m_nx > 1) ? m_Lx / (m_nx - 1) : 0.0; }
    double get_dy() const { return (m_ny > 1) ? m_Ly / (m_ny - 1) : 0.0; }
    
    int get_num_elements() const { return static_cast<int>(m_elements.size()); }
    int get_num_nodes() const { return static_cast<int>(m_nodes.size()); }

    std::vector<int> get_boundary_nodes(int edge_id) const;
    std::vector<std::array<double, 2>> get_element_coords(int elem_index) const;
    std::array<double, 2> get_node_coords(size_t node_index) const;
};