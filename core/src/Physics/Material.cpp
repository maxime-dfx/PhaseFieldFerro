#include "Physics/Material.h"
#include "IO/Datafile.h"

Material::Material(const Datafile& config) {
    // Initialisation des coefficients du modèle de Landau-Devonshire
    a0 = config.material.a0; 
    b1 = config.material.b1; b2 = config.material.b2; b3 = config.material.b3; 
    c1 = config.material.c1; c2 = config.material.c2; c3 = config.material.c3; 

    eps0 = config.material.eps0; 
    t = config.material.t; 
    P0 = config.material.P0; 
    mu_p = config.material.mu_p; 
    mu_v = config.material.mu_v; 
    xi = config.material.xi; 
    c0 = config.material.c0; 
    eta_k = config.material.eta_k; 

    // Coefficients d'ordre supérieur pour le modèle de Landau-Devonshire
    alpha_1 = config.material.alpha_1;
    alpha_11 = config.material.alpha_11;
    alpha_111 = config.material.alpha_111;
    alpha_1111 = config.material.alpha_1111;
    alpha_12 = config.material.alpha_12;
    alpha_112 = config.material.alpha_112;
    alpha_1112 = config.material.alpha_1112;
    alpha_1122 = config.material.alpha_1122;
}