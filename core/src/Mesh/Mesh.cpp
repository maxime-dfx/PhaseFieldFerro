#include "Mesh/Mesh.h"

Mesh::Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
           std::vector<Node> nodes, std::vector<Element> elements)
    : m_Lx(Lx), m_Ly(Ly), m_nx(nx), m_ny(ny), m_element_type(type),
      m_nodes(std::move(nodes)), m_elements(std::move(elements))
{
    // Pré-calcul et mise en cache des noeuds de bord à l'initialisation
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (m_nodes[i].ref_tag != 0) { // 0 représente généralement les noeuds internes
            m_boundary_nodes_cache[m_nodes[i].ref_tag].push_back(static_cast<int>(i));
        }
    }
}

const std::vector<int>& Mesh::get_boundary_nodes(int edge_id) const {
    auto it = m_boundary_nodes_cache.find(edge_id);
    if (it != m_boundary_nodes_cache.end()) {
        return it->second;
    }
    // Fallback static si le tag n'est pas trouvé (évite l'allocation d'un vecteur vide)
    static const std::vector<int> empty_vec;
    return empty_vec;
}

std::array<std::array<double, 2>, 8> Mesh::get_element_coords(int elem_index) const {
    assert(elem_index >= 0 && elem_index < m_elements.size() && "Index de l'element invalide");
    
    const Element& e = m_elements[elem_index];
    std::array<std::array<double, 2>, 8> coords{}; // Rempli de 0.0 par défaut
    
    int num_nodes = e.get_num_nodes();
    for (int i = 0; i < num_nodes; ++i) {
        int node_idx = e.get_node_index(i);
        coords[i] = {m_nodes[node_idx].x, m_nodes[node_idx].y};
    }
    
    return coords;
}

std::array<double, 2> Mesh::get_node_coords(size_t node_index) const {
    assert(node_index < m_nodes.size() && "Index du noeud invalide");
    return {m_nodes[node_index].x, m_nodes[node_index].y};
}