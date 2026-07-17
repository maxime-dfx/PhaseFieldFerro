#include "Core/MeshGeneratorGmsh.h"
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

Mesh MeshGeneratorGmsh::load_from_msh(const std::string& msh_path, double Lx, double Ly) {
    std::ifstream file(msh_path);
    if (!file.is_open()) {
        throw std::runtime_error("[MeshGeneratorGmsh] Impossible d'ouvrir le fichier : " + msh_path);
    }

    // --- 1. Table des noms physiques : tag -> nom ("left","right","top","bottom","domain") ---
    std::unordered_map<int, std::string> physical_names;

    // --- 2. Noeuds : id Gmsh (1-based) -> index local (0-based) ---
    std::unordered_map<int, int> gmsh_id_to_local;
    std::vector<Node> nodes;

    // --- 3. Elements bruts ---
    std::vector<RawElement> raw_elements;

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("$PhysicalNames") != std::string::npos) {
            int n_phys;
            file >> n_phys;
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
            file >> n_nodes;
            nodes.reserve(n_nodes);
            for (int k = 0; k < n_nodes; ++k) {
                int id;
                double x, y, z;
                file >> id >> x >> y >> z;

                Node node;
                node.x = x;
                node.y = y;
                node.z = 0.0;
                node.ref_tag = 0; 
                gmsh_id_to_local[id] = static_cast<int>(nodes.size());
                nodes.push_back(node);
            }
            std::getline(file, line);
        } else if (line.find("$Elements") != std::string::npos) {
            int n_elems;
            file >> n_elems;
            raw_elements.reserve(n_elems);
            for (int k = 0; k < n_elems; ++k) {
                int id, elm_type, n_tags;
                file >> id >> elm_type >> n_tags;

                std::vector<int> tags(n_tags);
                for (int t = 0; t < n_tags; ++t) file >> tags[t];
                int physical_tag = (n_tags > 0) ? tags[0] : 0;

                int n_nodes_elem = 0;
                if (elm_type == 1) n_nodes_elem = 2;      // ligne (bord)
                else if (elm_type == 2) n_nodes_elem = 3; // triangle
                else {
                    throw std::runtime_error(
                        "[MeshGeneratorGmsh] Type d'element Gmsh non supporte : " + std::to_string(elm_type) +
                        ". Seuls les lignes (1) et triangles (2) sont geres. "
                        "Verifie que Mesh.ElementOrder=1 et qu'aucune recombinaison en quads n'est active.");
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

    if (nodes.empty()) {
        throw std::runtime_error("[MeshGeneratorGmsh] Aucun noeud lu dans " + msh_path +
                                  " - verifie le format (MSH v2.2 ASCII attendu, "
                                  "genere avec 'gmsh -2 rectangle.geo -format msh2 -o rectangle.msh').");
    }

    // --- 4. Tag des noeuds de bord, en s'appuyant sur les Physical Line ---

    std::unordered_map<std::string, int> name_to_ref_tag = {
        {"left", 1}, {"right", 2}, {"bottom", 3}, {"top", 4}
    };

    for (const auto& re : raw_elements) {
        if (re.elm_type != 1) continue; 
        auto it_name = physical_names.find(re.physical_tag);
        if (it_name == physical_names.end()) continue;
        auto it_tag = name_to_ref_tag.find(it_name->second);
        if (it_tag == name_to_ref_tag.end()) continue;

        int ref_tag = it_tag->second;
        for (int gmsh_id : re.node_ids) {
            int local_idx = gmsh_id_to_local.at(gmsh_id);
            nodes[local_idx].ref_tag = ref_tag;
        }
    }

    // --- 5. Construction des elements triangulaires 
    std::vector<Element> elements;
    elements.reserve(raw_elements.size());

    for (const auto& re : raw_elements) {
        if (re.elm_type != 2) continue; 

        Element e;
        e.ref_tag = re.physical_tag; 
        e.num_nodes = 3; 

        for (int k = 0; k < 3; ++k) {
            e.node_indices[k] = gmsh_id_to_local.at(re.node_ids[k]);
        }

        elements.push_back(e);
    }

    if (elements.empty()) {
        throw std::runtime_error("[MeshGeneratorGmsh] Aucun triangle lu dans " + msh_path);
    }

    return Mesh(Lx, Ly, /*nx=*/0, /*ny=*/0,
                ElementType::TRIANGLE3, std::move(nodes), std::move(elements));
}