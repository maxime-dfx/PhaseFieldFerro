#include "Mesh/include/MeshGenerators.h"
#include "Utils/include/Profiling.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace {
    struct RawElement {
        int elm_type;
        int physical_tag;
        std::vector<int> node_ids;
    };
}

namespace MeshGenerators {

Mesh generate_structured(const MeshConfig& config) {
    PROFILE_ZONE_NC("MeshGenerators::generate_structured", PROFILE_COLOR_SEQUENTIAL);

    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> connectivity;

    int num_nodes_x = config.nx + 1;
    int num_nodes_y = config.ny + 1;
    int num_nodes = num_nodes_x * num_nodes_y;
    nodes.reserve(num_nodes);

    // 1. Génération des Noeuds (et affectation des tags de bord géométriques)
    for (int j = 0; j < num_nodes_y; ++j) {
        for (int i = 0; i < num_nodes_x; ++i) {
            Node node;
            node.x = i * config.dx;
            node.y = j * config.dy;

            int tag = 0; // Noeud interne (aucun bord)
            if (i == 0) tag = 1;                       // Bord gauche
            else if (i == num_nodes_x - 1) tag = 2;     // Bord droit
            else if (j == 0) tag = 3;                   // Bord bas
            else if (j == num_nodes_y - 1) tag = 4;     // Bord haut

            node.ref_tag = tag;
            nodes.push_back(node);
        }
    }

    // 2. Génération des Éléments
    if (config.element_type == ElementType::QUAD4) {
        for (int j = 0; j < config.ny; ++j) {
            for (int i = 0; i < config.nx; ++i) {
                int n1 = j * num_nodes_x + i;
                int n2 = j * num_nodes_x + (i + 1);
                int n3 = (j + 1) * num_nodes_x + (i + 1);
                int n4 = (j + 1) * num_nodes_x + i;

                Element e;
                e.num_nodes = 4;
                e.ref_tag = 0; // Matériau de base du domaine géométrique
                e.offset = static_cast<int>(connectivity.size());
                connectivity.insert(connectivity.end(), {n1, n2, n3, n4});
                elements.push_back(e);
            }
        }
    } else if (config.element_type == ElementType::TRIANGLE3) {
        for (int j = 0; j < config.ny; ++j) {
            for (int i = 0; i < config.nx; ++i) {
                int n1 = j * num_nodes_x + i;
                int n2 = j * num_nodes_x + (i + 1);
                int n3 = (j + 1) * num_nodes_x + (i + 1);
                int n4 = (j + 1) * num_nodes_x + i;

                Element t1, t2;
                t1.num_nodes = 3;
                t2.num_nodes = 3;
                t1.ref_tag = 0;
                t2.ref_tag = 0;

                t1.offset = static_cast<int>(connectivity.size());
                if ((i + j) % 2 == 0) {
                    connectivity.insert(connectivity.end(), {n1, n2, n3});
                } else {
                    connectivity.insert(connectivity.end(), {n1, n2, n4});
                }

                t2.offset = static_cast<int>(connectivity.size());
                if ((i + j) % 2 == 0) {
                    connectivity.insert(connectivity.end(), {n1, n3, n4});
                } else {
                    connectivity.insert(connectivity.end(), {n2, n3, n4});
                }

                elements.push_back(t1);
                elements.push_back(t2);
            }
        }
    }

    return Mesh(config.Lx, config.Ly, config.nx, config.ny, config.element_type,
                std::move(nodes), std::move(elements), std::move(connectivity));
}

Mesh load_from_gmsh(const std::string& msh_path, double Lx, double Ly) {
    PROFILE_ZONE_NC("MeshGenerators::load_from_gmsh", PROFILE_COLOR_SEQUENTIAL);

    std::ifstream file(msh_path);
    if (!file.is_open()) {
        throw std::runtime_error("[MeshGenerators::load_from_gmsh] Impossible d'ouvrir le fichier : " + msh_path);
    }

    std::unordered_map<int, std::string> physical_names;
    std::unordered_map<int, int> gmsh_id_to_local;
    std::vector<Node> nodes;
    std::vector<RawElement> raw_elements;
    std::string line;

    while (std::getline(file, line)) {
        // --- 0. Vérification stricte de la version de Gmsh ---
        if (line.find("$MeshFormat") != std::string::npos) {
            std::getline(file, line);
            std::stringstream ss(line);
            double version;
            ss >> version;
            if (version >= 3.0) {
                throw std::runtime_error(
                    "ERREUR FATALE : Le maillage est au format Gmsh v" + std::to_string(version) +
                    ". Ce code requiert le format v2.2.\n" +
                    "-> Re-generez le maillage en ajoutant l'option '-format msh2' a la commande Gmsh."
                );
            }
        }
        else if (line.find("$PhysicalNames") != std::string::npos) {
            int n_phys;
            if (!(file >> n_phys)) break;
            std::getline(file, line);
            for (int k = 0; k < n_phys; ++k) {
                int dim, tag;
                std::string name;
                file >> dim >> tag >> name;
                if (!name.empty() && name.front() == '"' && name.back() == '"') {
                    name = name.substr(1, name.size() - 2);
                }
                physical_names[tag] = name;
            }
            std::getline(file, line);
        } else if (line.find("$Nodes") != std::string::npos) {
            int n_nodes;
            if (!(file >> n_nodes)) break;
            nodes.reserve(n_nodes);
            for (int k = 0; k < n_nodes; ++k) {
                int id;
                double x, y, z;
                file >> id >> x >> y >> z;

                Node node;
                node.x = x;
                node.y = y;
                node.ref_tag = 0; // Valeur par défaut

                gmsh_id_to_local[id] = static_cast<int>(nodes.size());
                nodes.push_back(node);
            }
            std::getline(file, line);
        } else if (line.find("$Elements") != std::string::npos) {
            int n_elems;
            if (!(file >> n_elems)) break;
            raw_elements.reserve(n_elems);
            for (int k = 0; k < n_elems; ++k) {
                int id, elm_type, n_tags;
                file >> id >> elm_type >> n_tags;

                std::vector<int> tags(n_tags);
                for (int t = 0; t < n_tags; ++t) file >> tags[t];
                int physical_tag = (n_tags > 0) ? tags[0] : 0;

                int n_nodes_elem = 0;
                if (elm_type == 15) n_nodes_elem = 1;     // Point
                else if (elm_type == 1) n_nodes_elem = 2; // Segment
                else if (elm_type == 2) n_nodes_elem = 3; // Triangle
                else if (elm_type == 3) n_nodes_elem = 4; // Quadrangle
                else {
                    // Type inconnu : on lit et on ignore proprement
                    std::string dummy;
                    std::getline(file, dummy);
                    continue;
                }

                RawElement re;
                re.elm_type = elm_type;
                re.physical_tag = physical_tag;
                re.node_ids.resize(n_nodes_elem);
                for (int nnode = 0; nnode < n_nodes_elem; ++nnode) {
                    file >> re.node_ids[nnode];
                }
                raw_elements.push_back(re);
            }
            std::getline(file, line);
        }
    }
    file.close();

    // --- Tag des noeuds de bord ---
    for (const auto& re : raw_elements) {
        if (re.elm_type != 1) continue;
        for (int gmsh_id : re.node_ids) {
            auto it = gmsh_id_to_local.find(gmsh_id);
            if (it != gmsh_id_to_local.end()) {
                nodes[it->second].ref_tag = re.physical_tag;
            }
        }
    }

    // --- Construction des éléments et de la connectivité ---
    std::vector<Element> elements;
    std::vector<int> connectivity;
    elements.reserve(raw_elements.size());
    connectivity.reserve(raw_elements.size() * 3);

    ElementType mesh_elem_type = ElementType::TRIANGLE3;

    for (const auto& re : raw_elements) {
        // On ne garde que les Triangles (2) ou Quads (3)
        if (re.elm_type != 2 && re.elm_type != 3) continue;

        Element e;
        e.ref_tag = re.physical_tag;
        e.num_nodes = (re.elm_type == 2) ? 3 : 4;
        e.offset = static_cast<int>(connectivity.size());

        bool valid = true;
        for (int k = 0; k < e.num_nodes; ++k) {
            auto it = gmsh_id_to_local.find(re.node_ids[k]);
            if (it != gmsh_id_to_local.end()) {
                connectivity.push_back(it->second);
            } else {
                valid = false;
                break; // Noeud corrompu, on abandonne l'élément
            }
        }

        if (valid) {
            elements.push_back(e);
            if (re.elm_type == 3) mesh_elem_type = ElementType::QUAD4;
        } else {
            connectivity.resize(e.offset);
        }
    }

    return Mesh(Lx, Ly, 0, 0, mesh_elem_type, std::move(nodes), std::move(elements), std::move(connectivity));
}

} // namespace MeshGenerators
