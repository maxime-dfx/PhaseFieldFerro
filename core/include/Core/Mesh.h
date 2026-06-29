#pragma once
#include <vector>
#include <array>
#include "IO/Datafile.h"

struct Node {
    double x, y, z;
    int ref_tag; 
};

struct Quad {
    int n1, n2, n3, n4; 
    int ref_tag;
};

struct Triangle {
    int n1, n2, n3; 
    int ref_tag;
};

struct Element {
    std::vector<int> node_indices; 
    int ref_tag; 
    int get_num_nodes() const { return static_cast<int>(node_indices.size()); }
    const std::vector<int>& get_node_indices() const { return node_indices; }
};

class Mesh
{
private:
    const double L_x, L_y;
    const int n_x, n_y;
    ElementType m_element_type;

    std::vector<Node> m_nodes;
    std::vector<Element> m_elements;
public:
    Mesh(const Datafile& config);
    const std::vector<Node>& get_nodes() const { return m_nodes; }

    ElementType get_element_type() const { return m_element_type; }
    
    const std::vector<Element>& get_elements() const { return m_elements; }
    
    double get_Lx() const { return L_x; }
    double get_Ly() const { return L_y; }
    int get_nx() const { return n_x; }
    int get_ny() const { return n_y; }
    double get_dx() const { return L_x / (n_x - 1); }
    double get_dy() const { return L_y / (n_y - 1); }
    
    int get_num_elements() const { return static_cast<int>(m_elements.size()); }
    int get_num_nodes() const { return n_x * n_y; }

    std::vector<int> get_boundary_nodes(int edge_id) const;
    std::vector<std::array<double, 2>> get_element_coords(int elem_index) const;
};