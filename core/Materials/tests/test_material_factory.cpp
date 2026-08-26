#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Core/MaterialManager.h"

// =========================================================================
// Tests MaterialManager
// -------------------------------------------------------------------------
// Verifie que chaque MaterialType est bien route vers la classe concrete
// attendue (get_type() renvoie le type demande), et que les materiaux
// passifs (PureElastic, PolymerElastic, GrainBoundary) sont bien
// non-ferroelectriques (P=0 impose par PolarizationDofMapper), tandis que
// Ferroelectric reste actif. C'est le seul test qui exercerait une
// regression de dispatch dans MaterialManager::create sans avoir besoin de
// monter une simulation complete.
// =========================================================================

TEST_CASE("MaterialManager : type Ferroelectric instancie un materiau actif") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.type = MaterialType::Ferroelectric;
    
    MaterialManager manager;
    auto material = manager.create(cfg);

    CHECK_TRUE(material->get_type() == MaterialType::Ferroelectric);
    CHECK_TRUE(material->is_ferroelectric());
}

TEST_CASE("MaterialManager : type PureElastic instancie un materiau passif") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.type = MaterialType::PureElastic;
    
    MaterialManager manager;
    auto material = manager.create(cfg);

    CHECK_TRUE(material->get_type() == MaterialType::PureElastic);
    CHECK_TRUE(!material->is_ferroelectric());
    CHECK_TRUE(material->get_a0() == 0.0);
    CHECK_TRUE(material->get_mu_p() == 0.0);
}

TEST_CASE("MaterialManager : type PolymerElastic instancie un materiau passif avec permittivite") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.type = MaterialType::PolymerElastic;
    
    MaterialManager manager;
    auto material = manager.create(cfg);

    CHECK_TRUE(material->get_type() == MaterialType::PolymerElastic);
    CHECK_TRUE(!material->is_ferroelectric());
    // Contrairement a PureElastic, la permittivite reste non nulle : la
    // matrice reste un dielectrique valide pour l'equation de Poisson.
    CHECK_NEAR(material->compute_effective_permittivity(1.0, cfg.eta_k, false), cfg.eps0, 1e-12);
}

TEST_CASE("MaterialManager : type GrainBoundary instancie par defaut un materiau passif avec ses propres Gc/kappa/mu_v") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.type = MaterialType::GrainBoundary;
    // is_ferroelectric n'est PAS deduit de `type` pour GrainBoundary (seul
    // type dans ce cas, cf. Datafile::parse_materials_list) : il faut le
    // positionner explicitement, comme le ferait le TOML.
    cfg.is_ferroelectric = false;
    cfg.Gc = 0.5;    // volontairement different du materiau de volume
    cfg.kappa = 3.0;
    cfg.mu_v = 7.0;
    
    MaterialManager manager;
    auto material = manager.create(cfg);

    CHECK_TRUE(material->get_type() == MaterialType::GrainBoundary);
    CHECK_TRUE(!material->is_ferroelectric());
    CHECK_NEAR(material->get_Gc(), 0.5, 1e-12);
    CHECK_NEAR(material->get_kappa(), 3.0, 1e-12);
    CHECK_NEAR(material->get_mu_v(), 7.0, 1e-12);
    // Passif : mu_p/a0 neutralises, aucun terme de Landau/couplage.
    CHECK_TRUE(material->get_mu_p() == 0.0);
    CHECK_TRUE(material->get_a0() == 0.0);
}

TEST_CASE("MaterialManager : type GrainBoundary avec is_ferroelectric=true devient actif (Option A)") {
    // Un joint de grain "ferroelectrique affaibli" : memes Gc/kappa/mu_v
    // dedies, mais avec sa propre physique de Landau (a0/mu_p/b1/b2/b3/
    // alpha_*) plutot qu'un simple verrou P=0.
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.type = MaterialType::GrainBoundary;
    cfg.is_ferroelectric = true;
    cfg.Gc = 0.5;
    cfg.kappa = 3.0;
    cfg.mu_v = 7.0;
    cfg.mu_p = 2.0;
    cfg.a0 = 0.05;
    
    MaterialManager manager;
    auto material = manager.create(cfg);

    CHECK_TRUE(material->get_type() == MaterialType::GrainBoundary);
    CHECK_TRUE(material->is_ferroelectric());
    CHECK_NEAR(material->get_mu_p(), 2.0, 1e-12);
    CHECK_NEAR(material->get_a0(), 0.05, 1e-12);

    // Contrairement au mode passif, la polarisation effective n'est plus
    // ecrasee a zero.
    Eigen::Vector2d P(0.3, 0.1);
    Eigen::Vector2d P_eff = material->compute_effective_polarization(P, 1.0, cfg.eta_k, false);
    CHECK_TRUE(P_eff.norm() > 1e-12);
}

TEST_CASE("MaterialManager : type derive de is_ferroelectric legacy si absent du TOML (via Datafile)") {
    // Ce test documente le contrat de compatibilite ascendante porte par
    // Datafile::parse_materials_list, pas MaterialManager lui-meme : un
    // MaterialConfig construit a la main avec is_ferroelectric=false mais
    // type par defaut (Ferroelectric) N'EST PAS auto-corrige par
    // MaterialManager - c'est Datafile qui fait la derivation au parsing.
    // MaterialManager::create fait toujours confiance a `type` tel quel.
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    cfg.is_ferroelectric = false;
    cfg.type = MaterialType::PolymerElastic; // derivation deja faite en amont
    
    MaterialManager manager;
    auto material = manager.create(cfg);
    
    CHECK_TRUE(!material->is_ferroelectric());
}