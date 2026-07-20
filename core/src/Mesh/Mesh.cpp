#include "Mesh/Mesh.h"
#include <cmath>
#include <cstddef>

Mesh::Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
           std::vector<Node> nodes, std::vector<Element> elements)
    : m_Lx(Lx), m_Ly(Ly), m_nx(nx), m_ny(ny), m_element_type(type),
      m_nodes(std::move(nodes)), m_elements(std::move(elements)) 
{
}

std::vector<int> Mesh::get_boundary_nodes(int edge_id) const {
    std::vector<int> boundary_nodes;
    switch (edge_id) {
        case 1: // bas
            for (int i = 0; i < m_nx; ++i) boundary_nodes.push_back(i);
            break;
        case 2: // droite
            for (int j = 0; j < m_ny; ++j) boundary_nodes.push_back(j * m_nx + (m_nx - 1));
            break;
        case 3: // haut
            for (int i = 0; i < m_nx; ++i) boundary_nodes.push_back((m_ny - 1) * m_nx + i);
            break;
        case 4: // gauche
            for (int j = 0; j < m_ny; ++j) boundary_nodes.push_back(j * m_nx);
            break;
        default:
            break;
    }
    return boundary_nodes;
}

std::vector<std::array<double, 2>> Mesh::get_element_coords(int elem_index) const {
    const Element& e = m_elements[elem_index];
    std::vector<std::array<double, 2>> coords;
    coords.reserve(e.node_indices.size());
    for (int node_idx : e.node_indices) {
        coords.push_back({m_nodes[node_idx].x, m_nodes[node_idx].y});
    }
    return coords;
}

std::array<double, 2> Mesh::get_node_coords(size_t node_index) const {
    return { m_nodes[node_index].x, m_nodes[node_index].y };
}