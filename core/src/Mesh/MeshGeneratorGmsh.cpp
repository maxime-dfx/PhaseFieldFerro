// #include "Mesh/MeshGeneratorGmsh.h"
// #include <gmsh.h> 
// #include <stdexcept>
// #include <unordered_map>
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <stdexcept>
// #include <unordered_map>
// #include <vector>

// namespace {

// struct RawElement {
//     int elm_type;
//     int physical_tag;
//     std::vector<int> node_ids;
// };

// } 
// // Dans MeshGeneratorGmsh.cpp



// Mesh MeshGeneratorGmsh::generate_with_api(const MeshConfig& config) {
//     // 1. Initialiser le contexte Gmsh
//     gmsh::initialize();
//     gmsh::option::setNumber("General.Terminal", 0); // Désactive les logs verbeux de Gmsh

//     gmsh::model::add("api_rectangle");

//     double Lx = config.Lx;
//     double Ly = config.Ly;
//     double lc = config.mesh_size;

//     // 2. Définition géométrique (CAO) du rectangle
//     int p1 = gmsh::model::geo::addPoint(0.0, 0.0, 0.0, lc);
//     int p2 = gmsh::model::geo::addPoint(Lx,  0.0, 0.0, lc);
//     int p3 = gmsh::model::geo::addPoint(Lx,  Ly,  0.0, lc);
//     int p4 = gmsh::model::geo::addPoint(0.0, Ly,  0.0, lc);

//     int l1 = gmsh::model::geo::addLine(p1, p2); // bas
//     int l2 = gmsh::model::geo::addLine(p2, p3); // droite
//     int l3 = gmsh::model::geo::addLine(p3, p4); // haut
//     int l4 = gmsh::model::geo::addLine(p4, p1); // gauche

//     int cl = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
//     int surf = gmsh::model::geo::addPlaneSurface({cl});

//     gmsh::model::geo::synchronize();

//     // 3. Groupes physiques (Alignés avec Mesh.cpp : 1=bas, 2=droite, 3=haut, 4=gauche)
//     gmsh::model::addPhysicalGroup(1, {l1}, 1); gmsh::model::setPhysicalName(1, 1, "bottom");
//     gmsh::model::addPhysicalGroup(1, {l2}, 2); gmsh::model::setPhysicalName(1, 2, "right");
//     gmsh::model::addPhysicalGroup(1, {l3}, 3); gmsh::model::setPhysicalName(1, 3, "top");
//     gmsh::model::addPhysicalGroup(1, {l4}, 4); gmsh::model::setPhysicalName(1, 4, "left");
//     gmsh::model::addPhysicalGroup(2, {surf}, 0); gmsh::model::setPhysicalName(2, 0, "domain");

//     // Recombinaison automatique en QUAD si demandé par ton config.toml
//     if (config.element_type == ElementType::QUAD) {
//         gmsh::option::setNumber("Mesh.RecombineAll", 1);
//         gmsh::option::setNumber("Mesh.Algorithm", 8); // Algorithme Frontal-DelQuad
//     }

//     // 4. Générer le maillage 2D en mémoire
//     gmsh::model::mesh::generate(2);

//     // 5. Extraction des Nœuds depuis l'API de Gmsh
//     std::vector<size_t> nodeTags;
//     std::vector<double> coord;
//     std::vector<double> parametricCoord;
//     gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, 2, surf, true);

//     std::vector<Node> nodes(nodeTags.size());
//     std::unordered_map<size_t, int> gmsh_id_to_local;

//     for (size_t i = 0; i < nodeTags.size(); ++i) {
//         nodes[i].x = coord[3 * i];
//         nodes[i].y = coord[3 * i + 1];
//         nodes[i].z = 0.0;
//         nodes[i].ref_tag = 0; // Par défaut, nœud intérieur
//         gmsh_id_to_local[nodeTags[i]] = static_cast<int>(i);
//     }

//     // 6. Marquage des ref_tag des nœuds situés sur les frontières (Corrigé)
//     std::vector<int> boundary_lines = {l1, l2, l3, l4};
//     std::unordered_map<int, int> line_to_ref_tag = {{l1, 1}, {l2, 2}, {l3, 3}, {l4, 4}};

//     for (int line_id : boundary_lines) {
//         std::vector<size_t> bNodeTags;
//         std::vector<double> bCoord;
//         gmsh::model::mesh::getNodes(bNodeTags, bCoord, parametricCoord, 1, line_id, true);
        
//         int ref_tag = line_to_ref_tag[line_id];
//         for (size_t tag : bNodeTags) {
//             if (gmsh_id_to_local.find(tag) != gmsh_id_to_local.end()) {
//                 nodes[gmsh_id_to_local[tag]].ref_tag = ref_tag;
//             }
//         }
//     }

//     // 7. Extraction des Éléments (Triangles ou Quads) depuis l'API Gmsh
//     std::vector<int> elementTypes;
//     std::vector<std::vector<size_t>> elementTags;
//     std::vector<std::vector<size_t>> nodeTagsPerElement;
//     gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, 2, surf);

//     std::vector<Element> elements;
//     for (size_t t = 0; t < elementTypes.size(); ++t) {
//         int type = elementTypes[t];
//         if (type != 2 && type != 3) continue; // type 2 = triangle, type 3 = quadrangle

//         size_t num_elements_of_type = elementTags[t].size();
//         size_t num_nodes_per_elem = (type == 2) ? 3 : 4;

//         for (size_t e = 0; e < num_elements_of_type; ++e) {
//             Element elem;
//             elem.node_indices.resize(num_nodes_per_elem);
//             for (size_t n = 0; n < num_nodes_per_elem; ++n) {
//                 size_t gmsh_node_id = nodeTagsPerElement[t][e * num_nodes_per_elem + n];
//                 elem.node_indices[n] = gmsh_id_to_local.at(gmsh_node_id);
//             }
//             elem.ref_tag = 0; 
//             elements.push_back(elem);
//         }
//     }

//     // 8. Nettoyer et fermer la session Gmsh
//     gmsh::finalize();

//     std::cout << "[MeshGeneratorGmsh-API] Maillage genere en memoire via l'API. ("
//               << nodes.size() << " noeuds, " << elements.size() << " elements)" << std::endl;

//     // nx=0 et ny=0 car c'est un maillage non-structuré
//     return Mesh(Lx, Ly, 0, 0, config.element_type, std::move(nodes), std::move(elements));
// }

// Mesh MeshGeneratorGmsh::load_from_msh(const std::string& msh_path, double Lx, double Ly) {
//     std::ifstream file(msh_path);
//     if (!file.is_open()) {
//         throw std::runtime_error("[MeshGeneratorGmsh] Impossible d'ouvrir le fichier : " + msh_path);
//     }

//     // --- 1. Table des noms physiques : tag -> nom ("left","right","top","bottom","domain") ---
//     std::unordered_map<int, std::string> physical_names;

//     // --- 2. Noeuds : id Gmsh (1-based) -> index local (0-based) ---
//     std::unordered_map<int, int> gmsh_id_to_local;
//     std::vector<Node> nodes;

//     // --- 3. Elements bruts ---
//     std::vector<RawElement> raw_elements;

//     std::string line;
//     while (std::getline(file, line)) {
//         if (line.find("$PhysicalNames") != std::string::npos) {
//             int n_phys;
//             file >> n_phys;
//             std::getline(file, line); 
//             for (int k = 0; k < n_phys; ++k) {
//                 int dim, tag;
//                 std::string name;
//                 file >> dim >> tag >> name;
//                 if (!name.empty() && name.front() == '"' && name.back() == '"') {
//                     name = name.substr(1, name.size() - 2);
//                 }
//                 physical_names[tag] = name;
//             }
//             std::getline(file, line);
//         } else if (line.find("$Nodes") != std::string::npos) {
//             int n_nodes;
//             file >> n_nodes;
//             nodes.reserve(n_nodes);
//             for (int k = 0; k < n_nodes; ++k) {
//                 int id;
//                 double x, y, z;
//                 file >> id >> x >> y >> z;

//                 Node node;
//                 node.x = x;
//                 node.y = y;
//                 node.z = 0.0;
//                 node.ref_tag = 0; 
//                 gmsh_id_to_local[id] = static_cast<int>(nodes.size());
//                 nodes.push_back(node);
//             }
//             std::getline(file, line);
//         } else if (line.find("$Elements") != std::string::npos) {
//             int n_elems;
//             file >> n_elems;
//             raw_elements.reserve(n_elems);
//             for (int k = 0; k < n_elems; ++k) {
//                 int id, elm_type, n_tags;
//                 file >> id >> elm_type >> n_tags;

//                 std::vector<int> tags(n_tags);
//                 for (int t = 0; t < n_tags; ++t) file >> tags[t];
//                 int physical_tag = (n_tags > 0) ? tags[0] : 0;

//                 int n_nodes_elem = 0;
//                 if (elm_type == 1) n_nodes_elem = 2;      // ligne (bord)
//                 else if (elm_type == 2) n_nodes_elem = 3; // triangle
//                 else {
//                     throw std::runtime_error(
//                         "[MeshGeneratorGmsh] Type d'element Gmsh non supporte : " + std::to_string(elm_type) +
//                         ". Seuls les lignes (1) et triangles (2) sont geres. "
//                         "Verifie que Mesh.ElementOrder=1 et qu'aucune recombinaison en quads n'est active.");
//                 }

//                 RawElement re;
//                 re.elm_type = elm_type;
//                 re.physical_tag = physical_tag;
//                 re.node_ids.resize(n_nodes_elem);
//                 for (int nnode = 0; nnode < n_nodes_elem; ++nnode) {
//                     file >> re.node_ids[nnode];
//                 }
//                 raw_elements.push_back(re);
//             }
//             std::getline(file, line);
//         }
//     }
//     file.close();

//     if (nodes.empty()) {
//         throw std::runtime_error("[MeshGeneratorGmsh] Aucun noeud lu dans " + msh_path +
//                                   " - verifie le format (MSH v2.2 ASCII attendu, "
//                                   "genere avec 'gmsh -2 rectangle.geo -format msh2 -o rectangle.msh').");
//     }

//     // --- 4. Tag des noeuds de bord, en s'appuyant sur les Physical Line ---

//     std::unordered_map<std::string, int> name_to_ref_tag = {
//         {"left", 1}, {"right", 2}, {"bottom", 3}, {"top", 4}
//     };

//     for (const auto& re : raw_elements) {
//         if (re.elm_type != 1) continue; 
//         auto it_name = physical_names.find(re.physical_tag);
//         if (it_name == physical_names.end()) continue;
//         auto it_tag = name_to_ref_tag.find(it_name->second);
//         if (it_tag == name_to_ref_tag.end()) continue;

//         int ref_tag = it_tag->second;
//         for (int gmsh_id : re.node_ids) {
//             int local_idx = gmsh_id_to_local.at(gmsh_id);
//             nodes[local_idx].ref_tag = ref_tag;
//         }
//     }

//     // --- 5. Construction des elements triangulaires 
//     std::vector<Element> elements;
//     elements.reserve(raw_elements.size());

//     for (const auto& re : raw_elements) {
//         if (re.elm_type != 2) continue; 

//         Element e;
//         e.node_indices.resize(3);
//         for (int k = 0; k < 3; ++k) {
//             e.node_indices[k] = gmsh_id_to_local.at(re.node_ids[k]);
//         }
//         e.ref_tag = 0;
//         elements.push_back(e);
//     }

//     if (elements.empty()) {
//         throw std::runtime_error("[MeshGeneratorGmsh] Aucun triangle lu dans " + msh_path);
//     }

//     return Mesh(Lx, Ly, /*nx=*/0, /*ny=*/0,
//                 ElementType::TRIANGLE, std::move(nodes), std::move(elements));
// }