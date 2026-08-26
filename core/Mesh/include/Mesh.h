#pragma once
#include <vector>
#include <array>
#include <cassert>
#include <cstddef>
#include "Utils/include/Types.h"

struct Node {
    double x, y; 
    int ref_tag;
};

// Element ultra-leger : ne stocke que son offset dans le tableau de
// connectivite plat du Mesh (Data-Oriented Design), plus num_nodes et ref_tag.
struct Element {
    int offset = 0;
    int num_nodes = 0;
    int ref_tag = 0;

    int get_num_nodes() const { return num_nodes; }
};

class Mesh {
private:
    double m_Lx, m_Ly;
    int m_nx, m_ny;
    ElementType m_element_type;
    
    std::vector<Node> m_nodes;
    std::vector<Element> m_elements;

    // Tableau plat (Flat Array) contenant la connectivite globale de tous
    // les elements bout a bout : [e0_n0, e0_n1, ..., e1_n0, e1_n1, ...]
    std::vector<int> m_connectivity;
    std::vector<std::vector<int>> m_color_groups;

public:
    Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
         std::vector<Node> nodes, std::vector<Element> elements,
         std::vector<int> connectivity);

    const std::vector<Node>& get_nodes() const { return m_nodes; }
    const std::vector<Element>& get_elements() const { return m_elements; }
    const std::vector<int>& get_connectivity() const { return m_connectivity; }
    
    ElementType get_element_type() const { return m_element_type; }
    
    double get_Lx() const { return m_Lx; }
    double get_Ly() const { return m_Ly; }
    
    int get_num_elements() const { return static_cast<int>(m_elements.size()); }
    int get_num_nodes() const { return static_cast<int>(m_nodes.size()); }
    
    // Acces direct au noeud global d'un element (via le tableau plat)
    int get_node_index(int elem_index, int local_node_index) const {
        const Element& el = m_elements[elem_index];
        assert(local_node_index < el.num_nodes && "Index de noeud local hors limites !");
        return m_connectivity[el.offset + local_node_index];
    }

    // Pointeur brut vers le debut des noeuds de l'element (iteration rapide)
    const int* get_element_nodes_ptr(int elem_index) const {
        return &m_connectivity[m_elements[elem_index].offset];
    }

    std::array<std::array<double, 2>, 8> get_element_coords(int elem_index) const;
    std::array<double, 2> get_node_coords(size_t node_index) const;
    
    void apply_rcm();
    const std::vector<std::vector<int>>& get_color_groups() const { return m_color_groups; }
    void compute_coloring();
};
