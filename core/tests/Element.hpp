// #pragma once
// #include <vector>
// #include <Eigen/Dense>

// using namespace Eigen;

// // Représente un nœud du maillage
// struct Node {
//     double x, y;
// };

// // Représente un Élément Fini 2D (Triangle Linéaire T3)
// class Element {
// public:
//     std::vector<int> node_indices; // 3 indices globaux pour un triangle
//     std::vector<Node> coords;      // Coordonnées locales des 3 nœuds

//     // Points de Gauss et Poids (Règle à 3 points pour le triangle)
//     // Coordonnées barycentriques (milieux des arêtes)
//     const double gauss_pts[3][2] = {
//         {1.0 / 6.0, 1.0 / 6.0},
//         {2.0 / 3.0, 1.0 / 6.0},
//         {1.0 / 6.0, 2.0 / 3.0}
//     };
    
//     // La somme des poids d'un triangle de référence vaut son aire (1/2)
//     const double gauss_weights[3] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

//     int num_gauss_points() const { return 3; }

//     // Fonctions de forme N_i(xi, eta) pour un T3
//     VectorXd N(int q) const {
//         double xi = gauss_pts[q][0];
//         double eta = gauss_pts[q][1];
//         VectorXd n(3);
//         n(0) = 1.0 - xi - eta; // N1
//         n(1) = xi;             // N2
//         n(2) = eta;            // N3
//         return n;
//     }

//     // Dérivées des fonctions de forme par rapport à xi et eta (Constantes pour un T3)
//     MatrixXd dN_dxi(int /*q*/) const {
//         MatrixXd dn(2, 3);
//         // dN / dxi
//         dn(0, 0) = -1.0; 
//         dn(0, 1) =  1.0; 
//         dn(0, 2) =  0.0;
        
//         // dN / deta
//         dn(1, 0) = -1.0;  
//         dn(1, 1) =  0.0;  
//         dn(1, 2) =  1.0;
//         return dn;
//     }

//     // Calcul de la matrice Jacobienne J au point de Gauss q
//     Matrix2d jacobian(int q) const {
//         MatrixXd dn = dN_dxi(q);
//         Matrix2d J = Matrix2d::Zero();
//         for (int i = 0; i < 3; ++i) {
//             J(0, 0) += dn(0, i) * coords[i].x; J(0, 1) += dn(0, i) * coords[i].y;
//             J(1, 0) += dn(1, i) * coords[i].x; J(1, 1) += dn(1, i) * coords[i].y;
//         }
//         return J;
//     }

//     // Poids d'intégration physique : det(J) * Poids_Gauss
//     double weight(int q) const {
//         return jacobian(q).determinant() * gauss_weights[q];
//     }

//     // Gradient spatial réel : dN/dx, dN/dy = J^{-1} * dN/dxi
//     MatrixXd grad_N(int q) const {
//         Matrix2d J_inv = jacobian(q).inverse();
//         MatrixXd dn = dN_dxi(q);
//         return J_inv * dn; // Taille (2 lignes, 3 colonnes)
//     }

//     // Matrice B mécanique : relie le vecteur déplacement aux déformations (Taille 3x6)
//     MatrixXd B_mech(int q) const {
//         MatrixXd grad = grad_N(q);
//         MatrixXd B = MatrixXd::Zero(3, 6);
//         for (int i = 0; i < 3; ++i) {
//             B(0, 2 * i)     = grad(0, i); // dNi/dx
//             B(1, 2 * i + 1) = grad(1, i); // dNi/dy
//             B(2, 2 * i)     = grad(1, i); // dNi/dy
//             B(2, 2 * i + 1) = grad(0, i); // dNi/dx
//         }
//         return B;
//     }
// };