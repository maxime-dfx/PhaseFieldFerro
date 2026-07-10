// #include <iostream>
// #include <vector>
// #include <string>
// #include <exception>

// #include "Datafile.hpp"
// #include "Element.hpp"
// #include "GmshReader.hpp"
// #include "BoundaryManager.hpp"
// #include "FerroelectricMaterial.hpp"
// #include "MultiphysicsAssembler.hpp"
// #include "StaggeredSolver.hpp"
// #include "Mesh.hpp"

// int main(int argc, char** argv) {
//     std::string config_file = (argc > 1) ? argv[1] : "config.toml";

//     try {
//         std::cout << "1. Lecture de la configuration : " << config_file << std::endl;
//         Datafile config(config_file);

//         std::cout << "2. Initialisation du Maillage Externe..." << std::endl;
//         std::vector<Node> nodes;
//         std::vector<Element> elements;
        
//         std::string mesh_filename = config.simulation.mesh_create_file; 
//         if (mesh_filename.empty() || mesh_filename == "mesh") {
//             mesh_filename = "../input/rectangle.msh"; 
//         }
        
//         GmshReader::read_msh(mesh_filename, nodes, elements);

//         Mesh mesh(config.mesh.Lx, config.mesh.Ly, config.mesh.nx, config.mesh.ny, 
//                   config.mesh.element_type, nodes, elements);
                  
//         int total_nodes = mesh.get_num_nodes();

//         std::cout << "3. Traitement des Conditions aux Limites..." << std::endl;
//         BoundaryManager bc_manager(mesh, config);
//         bc_manager.initialize_all_boundaries();

//         std::cout << "4. Instanciation de la Physique et du Materiau..." << std::endl;
//         FerroelectricMaterial material(
//             config.material.Gc, 
//             config.material.kappa, 
//             config.material.eta_k, 
//             config.material.eps0, 
//             config.material.a0,
//             config.material.mu_p,
//             config.material.mu_v
//         );

//         MultiphysicsAssembler assembler(material, config.simulation.dt);

//         std::cout << "5. Lancement du Solveur Staggered..." << std::endl;
//         StaggeredSolver solver(assembler, bc_manager, config.simulation.dt, total_nodes);

//         double t_max = config.simulation.total_time * config.simulation.dt;
        
//         solver.run_simulation(t_max, mesh, config.simulation.output_dir);

//         std::cout << "=== Simulation terminee avec succes ! ===" << std::endl;

//     } catch (const std::exception& e) {
//         std::cerr << "\nERREUR FATALE : " << e.what() << std::endl;
//         return EXIT_FAILURE;
//     }

//     return EXIT_SUCCESS;
// }