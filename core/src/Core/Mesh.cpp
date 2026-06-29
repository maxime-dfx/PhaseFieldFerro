#include "Core/Mesh.h"
#include "IO/Datafile.h"

using namespace std;

Mesh::Mesh(const Datafile& config) 
    : L_x(config.get_L_x()), 
      L_y(config.get_L_y()), 
      n_x(config.get_n_x()), 
      n_y(config.get_n_y()),
      m_element_type(config.get_element_type()) 
{
    m_nodes.reserve(get_num_nodes());

    double dx = get_dx();
    double dy = get_dy();

    for (int j = 0; j < n_y; ++j) {
        for (int i = 0; i < n_x; ++i) {
            Node node;
            node.x = i * dx;
            node.y = j * dy;
            node.z = 0.0;
            node.ref_tag = 0; 
            m_nodes.push_back(node);
        }
    }

    if (m_element_type == ElementType::QUAD) {
        for (int j = 0; j < n_y - 1; ++j) {
            for (int i = 0; i < n_x - 1; ++i) {
                int n1 = j * n_x + i;
                int n2 = j * n_x + (i + 1);
                int n3 = (j + 1) * n_x + (i + 1);
                int n4 = (j + 1) * n_x + i;
                
                int tag = 0;
                if (i == 0) tag = 1; else if (i == n_x - 2) tag = 2;
                else if (j == 0) tag = 3; else if (j == n_y - 2) tag = 4;

                Element e;
                e.node_indices = {n1, n2, n3, n4};
                e.ref_tag = tag;
                m_elements.push_back(e);
            }
        }
    }
    if (m_element_type == ElementType::TRIANGLE) {
        for (int j = 0; j < n_y - 1; ++j) {
            for (int i = 0; i < n_x - 1; ++i) {
                int n1 = j * n_x + i;
                int n2 = j * n_x + (i + 1);
                int n3 = (j + 1) * n_x + (i + 1);
                int n4 = (j + 1) * n_x + i;

                int tag = 0;
                if (i == 0) tag = 1; else if (i == n_x - 2) tag = 2;
                else if (j == 0) tag = 3; else if (j == n_y - 2) tag = 4;

                Element t1;
                t1.node_indices = {n1, n2, n3};
                t1.ref_tag = tag;
                m_elements.push_back(t1);

                Element t2;
                t2.node_indices = {n1, n3, n4};
                t2.ref_tag = tag;
                m_elements.push_back(t2);
            }
        }
    }
}

std::vector<int> Mesh::get_boundary_nodes(int edge_id) const {
    std::vector<int> nodes;
    int nx = get_nx(), ny = get_ny();
    switch (edge_id) {
        case 1: // bas
            for (int i = 0; i < nx; ++i) nodes.push_back(i);
            break;
        case 2: // droite
            for (int j = 0; j < ny; ++j) nodes.push_back(j * nx + (nx - 1));
            break;
        case 3: // haut
            for (int i = 0; i < nx; ++i) nodes.push_back((ny - 1) * nx + i);
            break;
        case 4: // gauche
            for (int j = 0; j < ny; ++j) nodes.push_back(j * nx);
            break;
        default:
            break;
    }
    return nodes;
}

std::vector<std::array<double, 2>> Mesh::get_element_coords(int elem_index) const {
    const Element& e = m_elements[elem_index];
    std::vector<std::array<double, 2>> coords;
    for (int node_idx : e.node_indices) {
        coords.push_back({m_nodes[node_idx].x, m_nodes[node_idx].y});
    }
    return coords;
}
