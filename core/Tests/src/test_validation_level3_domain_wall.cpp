#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Physics/include/Core/PhysicsManager.h"
#include "Mesh/include/MeshGenerators.h"
#include "Physics/include/Modules/Polarization.h"
#include <cmath>
#include <fstream>
#include <cstdio>

TEST_CASE("Validation Niv3 : Profil analytique de paroi de domaine a 180 degres") {
    // 1. Initialisation de la configuration avec les constantes du papier
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    mat_cfg.alpha_1 = -0.0023; // Force la valeur pour le calcul analytique
    mat_cfg.a0 = 0.1;
    
    // Calcul de la solution analytique attendue
    double P0 = std::sqrt(std::abs(mat_cfg.alpha_1) / mat_cfg.alpha_11); // P_spontané simplifié
    double delta = std::sqrt(mat_cfg.a0 / (2.0 * std::abs(mat_cfg.alpha_1))); // Epaisseur paroi
    double x0 = 10.0; // Centre du maillage
    
    // 2. Configuration d'un maillage 1D-like (longueur 20, largeur 1)
    MeshConfig mesh_cfg;
    mesh_cfg.calcul_mesh = true;
    mesh_cfg.Lx = 20.0; mesh_cfg.Ly = 1.0;
    mesh_cfg.nx = 100;  mesh_cfg.ny = 1;
    mesh_cfg.element_type = ElementType::QUAD4;
    Mesh mesh = MeshGenerators::generate_structured(mesh_cfg);
    
    // 3. Setup des managers : Création du fichier TOML temporaire
    std::string path = "test_dummy.toml";
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n";
    out.close();
    
    Datafile config(path); 
    std::remove(path.c_str()); // Nettoyage immédiat
    
    config.simulation.dt_relax = 0.05;
    config.simulation.max_iter = 500;
    config.simulation.tol_ferro = 1e-4;
    config.physics_toggle.polarization = true;
    config.physics_toggle.mechanics = false;
    config.physics_toggle.electrostatics = false;
    config.physics_toggle.fracture = false;
    
    MaterialManager materials;
    materials.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    
    std::vector<BoundaryRuleConfig> empty_rules;
    BoundaryManager bc(mesh, empty_rules);
    Polycrystal polycrystal; // Vide (monocristal)
    
    PhysicsManager pm(config, mesh, bc, materials, polycrystal);
    
    // 4. Initialisation artificielle "marche d'escalier" pour forcer la paroi
    const PhysicsState& state = pm.get_current_state();
    Eigen::VectorXd* Py = const_cast<Eigen::VectorXd*>(state.get_field("Py"));
    for(int i = 0; i < mesh.get_num_nodes(); ++i) {
        auto coords = mesh.get_node_coords(i);
        (*Py)(i) = (coords[0] < x0) ? P0 : -P0; // +P0 à gauche, -P0 à droite
    }
    
    // 5. Relaxation du système jusqu'à l'équilibre
    int iters = pm.compute_one_step_physics(0.0, 0.1);
    CHECK_TRUE(iters > 0 && iters < config.simulation.max_iter); // Vérifie la convergence
    
    // 6. Validation contre la solution analytique
    double max_error = 0.0;
    for(int i = 0; i < mesh.get_num_nodes(); ++i) {
        auto coords = mesh.get_node_coords(i);
        double x = coords[0];
        
        // On évite les bords stricts pour les effets de bord de Neumann homogène
        if(x > 2.0 && x < 18.0) {
            double py_num = (*Py)(i);
            double py_ana = -P0 * std::tanh((x - x0) / delta); // Analytique
            double error = std::abs(py_num - py_ana);
            max_error = std::max(max_error, error);
        }
    }
    
    // Si l'erreur est petite, c'est que l'assemblage de a0 (gradient) et alpha (Landau) est PARFAIT
    CHECK_NEAR(max_error, 0.0, 5e-2); 
}