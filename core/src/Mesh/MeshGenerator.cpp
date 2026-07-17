// #include "Mesh/MeshGenerator.h"

// Mesh MeshGenerator::generate_structured_mesh(const MeshConfig& config) {
//     std::vector<Node> nodes;
//     std::vector<Element> elements;
    
//     int num_nodes = config.nx * config.ny;
//     nodes.reserve(num_nodes);

//     // 1. Génération des Nœuds
//     for (int j = 0; j < config.ny; ++j) {
//         for (int i = 0; i < config.nx; ++i) {
//             Node node;
//             node.x = i * config.dx;
//             node.y = j * config.dy;
//             node.z = 0.0;
//             node.ref_tag = 0; 
//             nodes.push_back(node);
//         }
//     }

//     // 2. Génération des Éléments
//     if (config.element_type == ElementType::QUAD) {
//         for (int j = 0; j < config.ny - 1; ++j) {
//             for (int i = 0; i < config.nx - 1; ++i) {
//                 int n1 = j * config.nx + i;
//                 int n2 = j * config.nx + (i + 1);
//                 int n3 = (j + 1) * config.nx + (i + 1);
//                 int n4 = (j + 1) * config.nx + i;
                
//                 int tag = 0;
//                 if (i == 0) tag = 1; else if (i == config.nx - 2) tag = 2;
//                 else if (j == 0) tag = 3; else if (j == config.ny - 2) tag = 4;

//                 Element e;
//                 e.node_indices = {n1, n2, n3, n4};
//                 e.ref_tag = tag;
//                 elements.push_back(e);
//             }
//         }
//     } else if (config.element_type == ElementType::TRIANGLE) {
//         for (int j = 0; j < config.ny - 1; ++j) {
//             for (int i = 0; i < config.nx - 1; ++i) {
//                 int n1 = j * config.nx + i;
//                 int n2 = j * config.nx + (i + 1);
//                 int n3 = (j + 1) * config.nx + (i + 1);
//                 int n4 = (j + 1) * config.nx + i;

//                 int tag = 0;
//                 if (i == 0) tag = 1; else if (i == config.nx - 2) tag = 2;
//                 else if (j == 0) tag = 3; else if (j == config.ny - 2) tag = 4;

//                 Element t1, t2;

//                 if ((i + j) % 2 == 0) {
//                     // Diagonale n1-n3 (comme avant)
//                     t1.node_indices = {n1, n2, n3};
//                     t2.node_indices = {n1, n3, n4};
//                 } else {
//                     // Diagonale n2-n4 (alternee)
//                     t1.node_indices = {n1, n2, n4};
//                     t2.node_indices = {n2, n3, n4};
//                 }

//                 t1.ref_tag = tag;
//                 t2.ref_tag = tag;
                
//                 elements.push_back(t1);
//                 elements.push_back(t2);
//             }
//         }
//     }

//     // 3. Retourne l'objet Mesh finalisé
//     return Mesh(config.Lx, config.Ly, config.nx, config.ny, 
//                 config.element_type, std::move(nodes), std::move(elements));
// }