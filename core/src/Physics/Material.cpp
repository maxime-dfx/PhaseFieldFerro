#include "Physics/Material.h"
#include "IO/Datafile.h"

Material::Material(const Datafile& config) {

    // Initialisation des coefficients du modèle de Landau-Devonshire
    a0 = config.get_a0(); // Coefficient linéaire
    b1 = config.get_b1(); b2 = config.get_b2(); b3 = config.get_b3(); // Coefficients quadratiques
    c1 = config.get_c1(); c2 = config.get_c2(); c3 = config.get_c3(); // Coefficients cubiques

    eps0 = config.get_eps0(); // Permittivité du vide
    t = config.get_t(); // Température
    P0 = config.get_P0(); // Polarisation de saturation
    mu_p = config.get_mu_p(); // Coefficient de viscosité
    mu_v = config.get_mu_v(); // Coefficient de viscosité
    xi = config.get_xi(); // Coefficient de couplage électromécanique
    c0 = config.get_c0(); // Coefficient de rigidité
    eta_k = config.get_eta_k(); // Coefficient de pénalisation pour la fracture

    // Coefficients d'ordre supérieur pour le modèle de Landau-Devonshire
    alpha_1 = config.get_alpha_1();
    alpha_11 = config.get_alpha_11();
    alpha_111 = config.get_alpha_111();
    alpha_1111 = config.get_alpha_1111();
    alpha_12 = config.get_alpha_12();
    alpha_112 = config.get_alpha_112();
    alpha_1112 = config.get_alpha_1112();
    alpha_1122 = config.get_alpha_1122();
}