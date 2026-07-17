
// #pragma once
// #include "Mesh/Mesh.h"
// #include "IO/Datafile.h" 
// #include "Mesh/MeshGenerator.h"
// #include "Mesh/MeshGeneratorGmsh.h"
// #include <stdexcept>

// class MeshFactory {
// public:
//     static Mesh build(const MeshConfig& config) {
//         if (config.source == MeshSource::INTERNAL_STRUCTURED) {
//             return MeshGenerator::generate_structured_mesh(config);
            
//         } else if (config.source == MeshSource::GMSH_FILE) {
//             return MeshGeneratorGmsh::load_from_msh(config.get_mesh_file, config.Lx, config.Ly);
            
//         } else if (config.source == MeshSource::GMSH_API) {
//             return MeshGeneratorGmsh::generate_with_api(config); // <--- Nouveau !
            
//         } else {
//             throw std::runtime_error("[MeshFactory] Source de maillage inconnue !");
//         }
//     }
// };