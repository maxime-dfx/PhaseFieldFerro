#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>
#include "IO/Datafile.h"
#include "Physics/MaterialModel.h"

class Mesh;
class Polarization;
class Mechanics;
class Fracture;
class Electrostatics;

// Résultat d'un test individuel
struct ValidationResult {
    std::string name;
    bool passed;
    double error;       // erreur relative ou absolue mesurée
    double tolerance;   // tolérance utilisée pour ce test
};

// Classe regroupant tous les tests de cohérence physique/numérique du modèle.
// Ne dépend d'aucun maillage ni assembleur : teste uniquement Math.h
// via comparaison analytique <-> différences finies.
class Validation {
public:
    explicit Validation(const MaterialModel& material);

    // Lance tous les tests et retourne les résultats
    std::vector<ValidationResult> run_all();

    // Affiche un résumé formaté dans la console
    void print_report(const std::vector<ValidationResult>& results) const;

    // Retourne true si tous les tests sont passés
    static bool all_passed(const std::vector<ValidationResult>& results);

    ValidationResult test_domain_wall_profile(const Datafile& config, const Mesh& mesh, Polarization& pol, Fracture& frac, Mechanics& mec, Electrostatics& elec);
    ValidationResult test_phase_field_crack_profile(const Datafile& config, const Mesh& mesh, Fracture& frac, Polarization& pol, Mechanics& mec, Electrostatics& elec);
    ValidationResult test_electrostatics_patch(const Datafile& config, const Mesh& mesh, Electrostatics& elec, Polarization& pol, Fracture& frac, Mechanics& mec);

private:
    const MaterialModel& material;

    static constexpr double FD_EPSILON = 1e-6;   // pas de différences finies
    static constexpr double FD_EPSILON_2ND = 1e-4;   // epsilon plus grand pour dérivées secondes (stabilité numérique)
    static constexpr double DEFAULT_TOL = 1e-4;  // tolérance relative par défaut

    // --- Tests individuels ---
    ValidationResult test_dchi_dp1();
    ValidationResult test_dchi_dp2();
    ValidationResult test_d2chi_dp1dp2_symmetry();
    ValidationResult test_dW_dp1();
    ValidationResult test_dW_dp2();
    ValidationResult test_d2W_dp1dp2_symmetry();
    ValidationResult test_sigma0_conjugate_to_W();
    ValidationResult test_H_drive_positivity();
    ValidationResult test_H_drive_zero_at_equilibrium();

    // --- Utilitaires génériques de différences finies ---

    // Dérivée par rapport à P(0) ou P(1) d'une fonction scalaire de P
    double finite_diff_dP(const std::function<double(const Eigen::Vector2d&)>& f,
                           const Eigen::Vector2d& P, int component) const;

    // Dérivée seconde croisée par rapport à P(0) et P(1)
    double finite_diff_d2P(const std::function<double(const Eigen::Vector2d&)>& f,
                            const Eigen::Vector2d& P) const;

    ValidationResult make_result(const std::string& name, double analytic, double numeric, double tol) const;

    // --- Nouveaux tests Niveau 1 : dynamique locale 0D (sans maillage) ---
    ValidationResult test_relaxation_to_minimum();
    ValidationResult test_stationary_equilibrium();
    ValidationResult test_switching_under_stress();

    // --- Niveau 1bis : patch test mécanique (1 élément QUAD) ---
    ValidationResult test_mechanics_patch_energy();

    // --- Utilitaire : Descente de gradient adaptative pour relaxation locale (0D) ---
    Eigen::Vector2d integrate_polarization_point(Eigen::Vector2d P, const Eigen::Matrix2d& strain,
                                                const Eigen::Vector2d& E, int n_steps) const;
                                                
    // --- Niveau 2 : Couplage de la fracture (Champ de phase v) ---
    ValidationResult test_fracture_penalty_GL_terms();
    ValidationResult test_crack_permeability_D_field();
    ValidationResult test_dh_dv_driving_force();
};