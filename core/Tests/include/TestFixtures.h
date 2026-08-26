#pragma once
#include "IO/include/ConfigTypes.h"

// TestFixtures
// -----------------------------------------------------------------------
// Donnees partagees entre les fichiers de test. Centralisees ici pour
// eviter que chaque fichier de test reconstruise sa propre copie
// (potentiellement divergente) des parametres materiau.
namespace TestFixtures {

// Jeu de coefficients de la Table 1 du papier (Abdollahi & Arias, Acta
// Materialia 59 (2011) 4733-4746) - identique a celui utilise dans les
// fichiers de configuration TOML du projet. MaterialConfig n'a pas
// d'initialisation par defaut (agregat simple), donc CHAQUE champ doit etre
// renseigne explicitement ici pour eviter un comportement indefini.
inline MaterialConfig make_paper_material_config() {
    MaterialConfig m{};
    m.alpha_1    = -0.0023;
    m.alpha_11   = -0.0029;
    m.alpha_12   = -0.0011;
    m.alpha_111  = 0.003;
    m.alpha_112  = -0.00068;
    m.alpha_1111 = 0.001;
    m.alpha_1112 = 0.0093;
    m.alpha_1122 = 1.24;
    m.b1 = 1.4282;
    m.b2 = -0.185;
    m.b3 = 0.8066;
    m.c1 = 185.0;
    m.c2 = 111.0;
    m.c3 = 74.0;
    m.eps0 = 0.131;
    m.t    = 300.0;
    m.P0   = 1.0;
    m.c0   = 1.0;
    m.xi   = 1.0;
    m.a0    = 0.1;
    m.Gc    = 4.0;
    m.kappa = 2.0;
    m.mu_p  = 1.0;
    m.mu_v  = 15.0;
    m.eta_k = 1e-6;
    return m;
}

} // namespace TestFixtures
