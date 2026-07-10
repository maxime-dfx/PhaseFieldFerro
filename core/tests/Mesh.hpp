// #pragma once
// #include <vector>
// #include "Element.hpp" // Doit inclure votre définition de Node et Element

// class Mesh {
// private:
//     double Lx_;
//     double Ly_;
//     int nx_;
//     int ny_;
//     int element_type_; // Non utilisé directement si on se fie aux vecteurs, mais gardé pour la compatibilité
    
//     std::vector<Node> nodes_;
//     std::vector<Element> elements_;

// public:
//     // Constructeur complet
//     // Note : Le type de retour ElementType a été remplacé par un int pour éviter l'erreur de type introuvable
//     Mesh(double Lx, double Ly, int nx, int ny, int element_type, 
//          const std::vector<Node>& nodes, const std::vector<Element>& elements)
//         : Lx_(Lx), Ly_(Ly), nx_(nx), ny_(ny), element_type_(element_type), 
//           nodes_(nodes), elements_(elements) {}

//     // Getters pour les dimensions (utiles pour BoundaryManager)
//     double get_Lx() const { return Lx_; }
//     double get_Ly() const { return Ly_; }
    
//     // Getters pour la géométrie et la connectivité
//     int get_num_nodes() const { return static_cast<int>(nodes_.size()); }
//     int get_num_elements() const { return static_cast<int>(elements_.size()); }
    
//     const std::vector<Node>& get_nodes() const { return nodes_; }
//     const std::vector<Element>& get_elements() const { return elements_; }
// };