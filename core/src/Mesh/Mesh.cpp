#include "Mesh/Mesh.h"
#include <queue>
#include <algorithm>
#include <numeric>
#include <stdexcept>

Mesh::Mesh(double Lx, double Ly, int nx, int ny, ElementType type, 
           std::vector<Node> nodes, std::vector<Element> elements)
    : m_Lx(Lx), m_Ly(Ly), m_nx(nx), m_ny(ny), m_element_type(type),
      m_nodes(std::move(nodes)), m_elements(std::move(elements))
{
    // Pré-calcul et mise en cache des noeuds de bord à l'initialisation
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (m_nodes[i].ref_tag != 0) { // 0 représente généralement les noeuds internes
            m_boundary_nodes_cache[m_nodes[i].ref_tag].push_back(static_cast<int>(i));
        }
    }
}

const std::vector<int>& Mesh::get_boundary_nodes(int edge_id) const {
    auto it = m_boundary_nodes_cache.find(edge_id);
    if (it != m_boundary_nodes_cache.end()) {
        return it->second;
    }
    // Fallback static si le tag n'est pas trouvé (évite l'allocation d'un vecteur vide)
    static const std::vector<int> empty_vec;
    return empty_vec;
}

std::array<std::array<double, 2>, 8> Mesh::get_element_coords(int elem_index) const {
    assert(elem_index >= 0 && elem_index < m_elements.size() && "Index de l'element invalide");
    
    const Element& e = m_elements[elem_index];
    std::array<std::array<double, 2>, 8> coords{}; // Rempli de 0.0 par défaut
    
    int num_nodes = e.get_num_nodes();
    for (int i = 0; i < num_nodes; ++i) {
        int node_idx = e.get_node_index(i);
        coords[i] = {m_nodes[node_idx].x, m_nodes[node_idx].y};
    }
    
    return coords;
}

std::array<double, 2> Mesh::get_node_coords(size_t node_index) const {
    assert(node_index < m_nodes.size() && "Index du noeud invalide");
    return {m_nodes[node_index].x, m_nodes[node_index].y};
}



// Algorithme de renumérotation Reverse Cuthill-McKee (RCM)

void Mesh::apply_rcm() {
    int num_nodes = static_cast<int>(m_nodes.size());
    if (num_nodes <= 2) return;

    // 1. Construction du graphe d'adjacence
    std::vector<std::vector<int>> adj(num_nodes);
    for (size_t e_idx = 0; e_idx < m_elements.size(); ++e_idx) {
        const auto& elem = m_elements[e_idx];
        int n_n = elem.get_num_nodes();
        for (int i = 0; i < n_n; ++i) {
            int u = elem.get_node_index(i);
            if (u < 0 || u >= num_nodes) {
                throw std::runtime_error("RCM Erreur : Indice de nœud u hors limites.");
            }
            for (int j = 0; j < n_n; ++j) {
                if (i == j) continue;
                int v = elem.get_node_index(j);
                if (v < 0 || v >= num_nodes) {
                    throw std::runtime_error("RCM Erreur : Indice de nœud voisin v hors limites.");
                }
                adj[u].push_back(v);
            }
        }
    }

    // Dédoublonnage des listes d'adjacence
    for (int i = 0; i < num_nodes; ++i) {
        auto& neighbors = adj[i];
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    // Calcul des degrés de chaque nœud
    std::vector<int> degree(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        degree[i] = static_cast<int>(adj[i].size());
    }

    // 2. Recherche d'un nœud de départ (heuristique du degré minimum)
    std::vector<bool> visited(num_nodes, false);
    std::vector<int> rcm_order;
    rcm_order.reserve(num_nodes);

    auto find_pseudo_peripheral = [&](int start) {
        int current = start;
        bool found_better = true;
        while (found_better) {
            found_better = false;
            std::queue<int> q;
            std::vector<int> dist(num_nodes, -1);
            q.push(current);
            dist[current] = 0;
            int farthest = current;

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v : adj[u]) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        q.push(v);
                        if (dist[v] > dist[farthest]) {
                            farthest = v;
                        }
                    }
                }
            }

            std::vector<int> candidates;
            int max_d = dist[farthest];
            for (int i = 0; i < num_nodes; ++i) {
                if (dist[i] == max_d) candidates.push_back(i);
            }

            int min_deg_node = candidates[0];
            for (int c : candidates) {
                if (degree[c] < degree[min_deg_node]) {
                    min_deg_node = c;
                }
            }

            if (degree[min_deg_node] < degree[current]) {
                current = min_deg_node;
                found_better = true;
            }
        }
        return current;
    };

    // 3. Parcours BFS avec tri des voisins par degré croissant
    for (int i = 0; i < num_nodes; ++i) {
        if (visited[i]) continue;

        int start_node = i;
        for (int j = 0; j < num_nodes; ++j) {
            if (!visited[j] && degree[j] < degree[start_node]) {
                start_node = j;
            }
        }

        start_node = find_pseudo_peripheral(start_node);

        std::queue<int> bfs_q;
        bfs_q.push(start_node);
        visited[start_node] = true;

        while (!bfs_q.empty()) {
            int u = bfs_q.front();
            bfs_q.pop();
            rcm_order.push_back(u);

            std::vector<int> unvisited_neighbors;
            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    unvisited_neighbors.push_back(v);
                }
            }

            std::sort(unvisited_neighbors.begin(), unvisited_neighbors.end(), [&](int a, int b) {
                return degree[a] < degree[b];
            });

            for (int v : unvisited_neighbors) {
                bfs_q.push(v);
            }
        }
    }

    // 4. Inversion de l'ordre (Reverse Cuthill-McKee)
    std::reverse(rcm_order.begin(), rcm_order.end());

    std::vector<int> old_to_new(num_nodes);
    for (int new_idx = 0; new_idx < num_nodes; ++new_idx) {
        old_to_new[rcm_order[new_idx]] = new_idx;
    }

    // 5. Application de la permutation (CORRECTION : Redimensionnement obligatoire !)
    std::vector<Node> old_nodes = std::move(m_nodes);
    m_nodes.resize(num_nodes); // <-- C'ICI QUE LE CODE PLANQUAIT LE SEGFAULT
    for (int old_idx = 0; old_idx < num_nodes; ++old_idx) {
        int new_idx = old_to_new[old_idx];
        m_nodes[new_idx] = old_nodes[old_idx];
    }

    for (auto& elem : m_elements) {
        int n_n = elem.get_num_nodes();
        for (int k = 0; k < n_n; ++k) {
            elem.node_indices[k] = old_to_new[elem.node_indices[k]];
        }
    }

    // Regénération du cache des bords
    m_boundary_nodes_cache.clear();
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (m_nodes[i].ref_tag != 0) {
            m_boundary_nodes_cache[m_nodes[i].ref_tag].push_back(static_cast<int>(i));
        }
    }
}