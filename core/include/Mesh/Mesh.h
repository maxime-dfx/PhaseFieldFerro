// #pragma once
// #include <vector>
// #include <array>
// #include "Core/Types.h" // On inclut juste les types de base, pas Datafile !

// struct Node {
//     double x, y, z;
//     int ref_tag; 
// };

// struct Element {
//     std::vector<int> node_indices; 
//     int ref_tag; 
//     int get_num_nodes() const { return static_cast<int>(node_indices.size()); }
//     const std::vector<int>& get_node_indices() const { return node_indices; }
// };

// class Mesh {
// private:
//     double m_Lx, m_Ly;
//     int m_nx, m_ny;
//     ElementType m_element_type;
//     std::vector<Node> m_nodes;
//     std::vector<Element> m_elements;

// public:
//     Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
//          std::vector<Node> nodes, std::vector<Element> elements);

//     const std::vector<Node>& get_nodes() const { return m_nodes; }
//     const std::vector<Element>& get_elements() const { return m_elements; }
    
//     ElementType get_element_type() const { return m_element_type; }
//     double get_Lx() const { return m_Lx; }
//     double get_Ly() const { return m_Ly; }

        
//     int get_nx() const { return m_nx; }
//     int get_ny() const { return m_ny; }
    

//     double get_dx() const { return (m_nx > 1) ? m_Lx / (m_nx - 1) : 0.0; }
//     double get_dy() const { return (m_ny > 1) ? m_Ly / (m_ny - 1) : 0.0; }
    
//     int get_num_elements() const { return static_cast<int>(m_elements.size()); }


//     int get_num_nodes() const { return static_cast<int>(m_nodes.size()); }

//     std::vector<int> get_boundary_nodes(int edge_id) const;
//     std::vector<std::array<double, 2>> get_element_coords(int elem_index) const;
//     std::array<double, 2> get_node_coords(size_t node_index) const;
// };