#pragma once

#include <string>
#include "Mesh/include/Mesh.h"
#include "IO/include/Datafile.h"

// Remplace IMeshProvider / MeshManager / StructuredMeshGenerator / GmshReader.
// Un maillage est un simple objet de données (Mesh) produit par une fonction pure ;
// aucune hiérarchie objet n'est nécessaire pour ça.
namespace MeshGenerators {

    // Génère un maillage structuré régulier (quad ou triangle) à partir d'une MeshConfig.
    Mesh generate_structured(const MeshConfig& config);

    // Charge un maillage depuis un fichier Gmsh (.msh, format v2.2).
    Mesh load_from_gmsh(const std::string& msh_path, double Lx, double Ly);

} // namespace MeshGenerators
