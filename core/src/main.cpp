#include "Simulation.h"
#include "IO/Datafile.h"
#include "Mesh/Mesh.h"
#include "Mesh/MeshGenerator.h" 
#include "Mesh/MeshGeneratorGmsh.h"
#include "IO/ResultsExporter.h"
#include "Utils/Logger.h"
#include "Physics/FerroelectricMaterial.h"
#include "Physics/MaterialModel.h"
#include <string>
#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>

// Macro pour l'affichage des tests
#define CHECK_TOLERANCE(nom, fd, analytique, tol) \
    std::cout << std::left << std::setw(30) << nom \
              << " | FD: " << std::setw(12) << fd \
              << " | Analytique: " << std::setw(12) << analytique \
              << " | Err: " << std::abs(fd - analytique) \
              << (std::abs(fd - analytique) < tol ? " [OK]" : " [FAIL]") << "\n";

#define CHECK_PARAM(nom, valeur_code, valeur_papier) \
    do { \
        double diff = std::abs((valeur_code) - (valeur_papier)); \
        bool ok = (diff < 1e-5); \
        std::cout << std::left << std::setw(30) << nom \
                  << " | Code: " << std::setw(10) << (valeur_code) \
                  << " | Papier: " << std::setw(10) << (valeur_papier) \
                  << (ok ? " [OK]" : " [ERREUR]") << "\n"; \
    } while(0)

// =========================================================================
// 1. VÉRIFICATION DES PARAMÈTRES
// =========================================================================
void test_paper_configuration(const Datafile& config, const Mesh& mesh) {
    std::cout << "========================================================\n";
    std::cout << "       VERIFICATION DES PARAMETRES (Arias 2011)         \n";
    std::cout << "========================================================\n";

    std::cout << "\n--- 1. Geometrie et Maillage (Sec. 3.1) ---\n";
    CHECK_PARAM("Dimension Lx", config.mesh.Lx, 200.0);
    CHECK_PARAM("Dimension Ly", config.mesh.Ly, 200.0);
    
    size_t num_elements = mesh.get_num_elements();
    std::cout << std::left << std::setw(30) << "Nombre d'elements" 
              << " | Code: " << std::setw(10) << num_elements 
              << " | Papier: ~80000   " 
              << ((num_elements > 75000 && num_elements < 85000) ? " [OK]" : " [ATTENTION]") << "\n";

    std::cout << "\n--- 2. Parametres de Fracture et Phase-Field ---\n";
    CHECK_PARAM("Gc normalise", config.material.Gc, 4.0);
    CHECK_PARAM("Regularisation kappa", config.material.kappa, 2.0);
    CHECK_PARAM("Raideur residuelle eta_k", config.material.eta_k, 1e-6);
    CHECK_PARAM("Largeur de paroi a0", config.material.a0, 0.1);
    CHECK_PARAM("Inverse mobilite P (mu_p)", config.material.mu_p, 1.0); 
    CHECK_PARAM("Inverse mobilite v (mu_v)", config.material.mu_v, 15.0);

    std::cout << "\n--- 3. Parametres Numeriques (Algo 1) ---\n";
    CHECK_PARAM("Pas de temps global (dt)", config.simulation.dt, 0.03);
    CHECK_PARAM("Tolerance ferro", config.simulation.tol_ferro, 1e-3);
    CHECK_PARAM("Tolerance vfield", config.simulation.tol_vfield, 1e-3);
    
    std::cout << "\n--- 4. Conditions Initiales ---\n";
    CHECK_PARAM("Polarisation initiale Px", config.polarization.val_x_0, 1.0);
    CHECK_PARAM("Polarisation initiale Py", config.polarization.val_y_0, 0.0);

    std::cout << "\n--- 5. Coefficients de Landau (Tableau 1) ---\n";
    CHECK_PARAM("c1", config.material.c1, 185.0);
    CHECK_PARAM("c2", config.material.c2, 111.0);
    CHECK_PARAM("c3", config.material.c3, 74.0);
    CHECK_PARAM("b1", config.material.b1, 1.4282);
    CHECK_PARAM("b2", config.material.b2, -0.185);
    CHECK_PARAM("b3", config.material.b3, 0.8066);
    CHECK_PARAM("alpha_1", config.material.alpha_1, -0.0023);
    CHECK_PARAM("alpha_11", config.material.alpha_11, -0.0029);
    CHECK_PARAM("alpha_12", config.material.alpha_12, -0.0011);
    CHECK_PARAM("alpha_111", config.material.alpha_111, 0.003);
    CHECK_PARAM("alpha_112", config.material.alpha_112, -0.00068);
    CHECK_PARAM("alpha_1111", config.material.alpha_1111, 0.001);
    CHECK_PARAM("alpha_1112", config.material.alpha_1112, 0.0093);
    CHECK_PARAM("alpha_1122", config.material.alpha_1122, 1.24);
    CHECK_PARAM("Permittivite eps0", config.material.eps0, 0.131);

    std::cout << "========================================================\n\n";
}


// =========================================================================
// 2. VÉRIFICATION DES DÉRIVÉES (APPEL DIRECT AUX MÉTHODES DE MATH.H)
// =========================================================================
void run_all_physics_tests(const MaterialModel& material, const Datafile& config) {
    std::cout << "========================================================\n";
    std::cout << "       VERIFICATION MULTIPHYSIQUE (Arias 2011)          \n";
    std::cout << "========================================================\n";

    // Etat arbitraire de test (non-nul pour éviter les annulations triviales)
    Eigen::Vector2d P(0.8, 0.3);
    Eigen::Matrix2d strain;
    strain << 0.01, -0.005,
             -0.005, 0.02;
    Eigen::Vector2d E(0.001, -0.002);
    
    double v = 0.5; // Endommagement partiel
    double penalite = (v * v) + material.get_eta_k();
    bool is_impermeable = true; // Mode testé
    
    double delta = 1e-6;
    double tol = 1e-4; 

    // ---------------------------------------------------------
    // 1. FORCES DE POLARISATION (Dérivées premières de l'Enthalpie)
    // ---------------------------------------------------------
    std::cout << "\n--- 1. Forces de Ginzburg-Landau (F_px, F_py) ---\n";
    
    // Enthalpie totale H = penalite * W + chi - E.P (Modèle de ton Math.h)
    auto compute_H_tot = [&](const Eigen::Vector2d& p) {
        return penalite * material.W_energy(p, strain) + material.chi_energy(p) - p.dot(E);
    };

    double H_p1_plus  = compute_H_tot(P + Eigen::Vector2d(delta, 0));
    double H_p1_minus = compute_H_tot(P - Eigen::Vector2d(delta, 0));
    double dH_dp1_FD  = (H_p1_plus - H_p1_minus) / (2.0 * delta);
    
    double H_p2_plus  = compute_H_tot(P + Eigen::Vector2d(0, delta));
    double H_p2_minus = compute_H_tot(P - Eigen::Vector2d(0, delta));
    double dH_dp2_FD  = (H_p2_plus - H_p2_minus) / (2.0 * delta);

    // On récupère les forces calculées par ta méthode :
    GinzburgLandauTerms GL = material.compute_GL_terms(P, strain, E, penalite, is_impermeable);

    CHECK_TOLERANCE("Force Px (dH/dP1)", dH_dp1_FD, GL.force_px, tol);
    CHECK_TOLERANCE("Force Py (dH/dP2)", dH_dp2_FD, GL.force_py, tol);

    // ---------------------------------------------------------
    // 2. JACOBIEN (Dérivées secondes de l'Enthalpie)
    // ---------------------------------------------------------
    std::cout << "\n--- 2. Jacobien de Ginzburg-Landau (J_11, J_22, J_12) ---\n";
    
    // Le Jacobien est la dérivée des forces. On rappelle math.compute_GL_terms avec FD.
    auto GL_p1_plus  = material.compute_GL_terms(P + Eigen::Vector2d(delta, 0), strain, E, penalite, is_impermeable);
    auto GL_p1_minus = material.compute_GL_terms(P - Eigen::Vector2d(delta, 0), strain, E, penalite, is_impermeable);
    double J11_FD = (GL_p1_plus.force_px - GL_p1_minus.force_px) / (2.0 * delta);
    double J12_FD = (GL_p1_plus.force_py - GL_p1_minus.force_py) / (2.0 * delta);

    auto GL_p2_plus  = material.compute_GL_terms(P + Eigen::Vector2d(0, delta), strain, E, penalite, is_impermeable);
    auto GL_p2_minus = material.compute_GL_terms(P - Eigen::Vector2d(0, delta), strain, E, penalite, is_impermeable);
    double J22_FD = (GL_p2_plus.force_py - GL_p2_minus.force_py) / (2.0 * delta);

    // On vérifie que ta matrice analytique correspond bien :
    CHECK_TOLERANCE("J_11 (d2H/dP1^2)", J11_FD, GL.J_11, tol);
    CHECK_TOLERANCE("J_22 (d2H/dP2^2)", J22_FD, GL.J_22, tol);
    CHECK_TOLERANCE("J_12 (d2H/dP1dP2)", J12_FD, GL.J_12, tol);

    // ---------------------------------------------------------
    // 3. CONTRAINTE ÉLASTIQUE (Dérivée de W par rapport à strain)
    // ---------------------------------------------------------
    std::cout << "\n--- 3. Contrainte Mecanique (Sigma) ---\n";
    
    // FD par rapport à eps_11
    Eigen::Matrix2d strain_eps11_plus = strain; strain_eps11_plus(0,0) += delta;
    Eigen::Matrix2d strain_eps11_minus = strain; strain_eps11_minus(0,0) -= delta;
    double dW_deps11_FD = (material.W_energy(P, strain_eps11_plus) - material.W_energy(P, strain_eps11_minus)) / (2.0 * delta);
    
    // Dans Math.h, compute_sigma_0 ne renvoie que la partie piezo (couplage).
    // La contrainte totale est sigma_0 + C * eps.
    Eigen::Vector3d sigma_0 = material.compute_sigma_0(P);
    Eigen::Matrix3d C = material.get_elastic_matrix();
    double sigma_11_total_analytic = sigma_0(0) + C(0,0) * strain(0,0) + C(0,1) * strain(1,1);
    
    CHECK_TOLERANCE("Sigma 11 (dW/deps11)", dW_deps11_FD, sigma_11_total_analytic, tol);

    std::cout << "========================================================\n\n";
}


int main(int argc, char* argv[]) {
    try {
        std::string config_file = (argc > 1) ? argv[1] : "../config.toml";
        
        Datafile config(config_file);
        FerroelectricMaterial material(config.material);
        Logger::set_level(config.simulation.debug_enabled ? LogLevel::DEBUG : LogLevel::INFO);
        Logger::info("Configuration loaded from: " + config_file);

        Mesh mesh = [&]() {
            if (config.mesh.calcul_mesh) {
                Logger::info("Mesh generated successfully.");
                return MeshGenerator::generate_structured_mesh(config.mesh);
            } else {
                Logger::info("Mesh loaded from file: " + config.mesh.get_mesh_file);
                return MeshGeneratorGmsh::load_from_msh(config.mesh.get_mesh_file, config.mesh.Lx, config.mesh.Ly);
            }
        }(); 

        ResultsExporter exporter(mesh);
        
        // -----------------------------------------------------------
        // TESTS DE ROBUSTESSE ET DE VALIDATION PHYSIQUE
        // -----------------------------------------------------------
        test_paper_configuration(config, mesh);
        // On passe directement l'instance de FerroelectricMaterial, qui hérite de MaterialModel
        run_all_physics_tests(material, config);
        // -----------------------------------------------------------
        
        // Retrait de 'material' des arguments car Simulation le génère en interne ou ne prend que 3 paramètres
        Simulation sim(config, mesh, exporter, material);
        sim.initialize_physics();
        sim.run();
        
    } catch (const std::exception& e) {
        std::cerr << "[ERREUR FATALE] " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}