#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "Element.hpp" 

class GmshReader {
public:
    static void read_msh(const std::string& filename, std::vector<Node>& out_nodes, std::vector<Element>& out_elements) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Erreur: Impossible d'ouvrir le fichier maillage " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            
            // --- 1. Lecture des Noeuds ---
            if (line.find("$Nodes") != std::string::npos) {
                int num_nodes;
                file >> num_nodes;
                out_nodes.resize(num_nodes);
                
                int id;
                double x, y, z;
                for (int i = 0; i < num_nodes; ++i) {
                    file >> id >> x >> y >> z;
                    // L'identifiant Gmsh commence a 1, nos vecteurs commencent a 0
                    out_nodes[id - 1] = {x, y}; 
                }
            }
            
            // --- 2. Lecture des Elements ---
            else if (line.find("$Elements") != std::string::npos) {
                int num_elements;
                file >> num_elements;
                
                int id, type, num_tags;
                for (int i = 0; i < num_elements; ++i) {
                    file >> id >> type >> num_tags;
                    
                    // Passer les métadonnées (tags) propres à GMSH
                    int tag;
                    for (int j = 0; j < num_tags; ++j) file >> tag;

                    // Type 2 = Triangle T3 (3 noeuds)
                    if (type == 2) {
                        Element el;
                        el.node_indices.resize(3);
                        el.coords.resize(3);
                        
                        for (int k = 0; k < 3; ++k) {
                            int node_id;
                            file >> node_id;
                            el.node_indices[k] = node_id - 1; 
                            el.coords[k] = out_nodes[node_id - 1]; // Geometrie locale
                        }
                        out_elements.push_back(el);
                    } 
                    // Type 1 (Ligne) ou Type 15 (Point) : On les ignore pour l'assemblage volumique
                    else {
                        std::string dummy;
                        std::getline(file, dummy); 
                    }
                }
            }
        }
        std::cout << "-> Maillage TRI charge : " << out_nodes.size() << " noeuds, " 
                  << out_elements.size() << " elements T3." << std::endl;
    }
};