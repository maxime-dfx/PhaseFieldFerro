#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Materials/Models/Ferroelectric.h"
#include "Materials/Microstructure/Crystal.h"
#include <cmath>

// =========================================================================
// Tests OrientedMaterial
// -------------------------------------------------------------------------
// C'est la partie la plus a risque du module polycristal : rotation des
// champs (P, strain, E, grad_P) vers le repere local du grain, appel des
// formules SingleCrystalMaterial inchangees, puis rotation retour des
// resultats (forces, contraintes, tenseur de rigidite). Deux familles de
// tests, chacune independante de "faire confiance a l'implementation" :
//
//  1) theta=0 doit etre l'IDENTITE stricte par rapport a SingleCrystalMaterial
//     utilise directement (a theta=0, R(theta)=Identite, donc to_local/
//     to_global sont des no-op). C'est le test de non-regression le plus
//     important : si OrientedMaterial(theta=0) ne matche pas exactement
//     SingleCrystalMaterial, il y a un bug dans la logique de rotation
//     elle-meme, independamment de toute question de convention.
//
//  2) Periodicite en pi : la physique de ce modele (Abdollahi & Arias) est
//     construite exclusivement a partir de puissances PAIRES de P (Landau
//     centrosymetrique - cf Table 1, aucun terme alpha_* n'est de degre
//     impair) et de couplages bilineaires pairs (P⊗P pour sigma_0, P.E pour
//     l'electrostatique). Consequence physique : un grain d'orientation
//     theta et un grain d'orientation theta+pi sont IDENTIQUES pour ce
//     modele - toutes les grandeurs publiques d'OrientedMaterial doivent
//     donc etre rigoureusement pi-periodiques en theta. C'est un invariant
//     fort, derivable independamment du code (cf. notes de derivation),
//     qui detecterait la plupart des erreurs de signe dans les rotations.
// =========================================================================

namespace {

const double kTol = 1e-9;

} // namespace

TEST_CASE("OrientedMaterial : theta=0 est l'identite stricte") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial base(cfg);
    OrientedMaterial oriented(base, 0.0);

    CHECK_NEAR(oriented.get_eta_k(), base.get_eta_k(), kTol);
    CHECK_NEAR(oriented.get_eps0(), base.get_eps0(), kTol);

    // Plusieurs points d'essai varies (dont un P asymetrique et un P nul).
    const std::vector<Eigen::Vector2d> P_samples = {
        {0.3, -0.2}, {1.0, 0.0}, {0.0, 0.0}, {-0.5, 0.7}
    };
    const std::vector<Eigen::Matrix2d> strain_samples = {
        (Eigen::Matrix2d() << 0.01, 0.002, 0.002, -0.015).finished(),
        Eigen::Matrix2d::Zero()
    };
    const Eigen::Vector2d E(0.1, -0.05);
    const Eigen::Matrix2d grad_P = (Eigen::Matrix2d() << 0.02, -0.01, 0.03, 0.015).finished();
    const double v = 0.8;
    const double penalite = v * v + cfg.eta_k;

    for (const auto& P : P_samples) {
        for (const auto& strain : strain_samples) {
            CHECK_NEAR(oriented.W_energy(P, strain), base.W_energy(P, strain), kTol);
            CHECK_NEAR(oriented.chi_energy(P), base.chi_energy(P), kTol);

            for (bool is_imp : {true, false}) {
                CHECK_NEAR(oriented.compute_H_drive(grad_P, P, strain, E, is_imp),
                           base.compute_H_drive(grad_P, P, strain, E, is_imp), kTol);

                auto GL_o = oriented.compute_GL_terms(P, strain, E, penalite, is_imp);
                auto GL_b = base.compute_GL_terms(P, strain, E, penalite, is_imp);
                CHECK_NEAR(GL_o.force_px, GL_b.force_px, kTol);
                CHECK_NEAR(GL_o.force_py, GL_b.force_py, kTol);
                CHECK_NEAR(GL_o.J_11, GL_b.J_11, kTol);
                CHECK_NEAR(GL_o.J_12, GL_b.J_12, kTol);
                CHECK_NEAR(GL_o.J_22, GL_b.J_22, kTol);
            }
        }
        CHECK_NEAR(oriented.U_energy(grad_P), base.U_energy(grad_P), kTol);

        Eigen::Vector3d sigma_o = oriented.compute_sigma_0(P);
        Eigen::Vector3d sigma_b = base.compute_sigma_0(P);
        CHECK_NEAR(sigma_o(0), sigma_b(0), kTol);
        CHECK_NEAR(sigma_o(1), sigma_b(1), kTol);
        CHECK_NEAR(sigma_o(2), sigma_b(2), kTol);

        for (bool is_imp : {true, false}) {
            Eigen::Vector2d Peff_o = oriented.compute_effective_polarization(P, v, cfg.eta_k, is_imp);
            Eigen::Vector2d Peff_b = base.compute_effective_polarization(P, v, cfg.eta_k, is_imp);
            CHECK_NEAR(Peff_o(0), Peff_b(0), kTol);
            CHECK_NEAR(Peff_o(1), Peff_b(1), kTol);
        }
    }

    for (bool is_imp : {true, false}) {
        CHECK_NEAR(oriented.compute_effective_permittivity(v, cfg.eta_k, is_imp),
                   base.compute_effective_permittivity(v, cfg.eta_k, is_imp), kTol);
    }

    Eigen::Matrix3d C_o = oriented.get_elastic_matrix();
    Eigen::Matrix3d C_b = base.get_elastic_matrix();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK_NEAR(C_o(i,j), C_b(i,j), kTol);
}

TEST_CASE("OrientedMaterial : pi-periodicite en theta (physique centrosymetrique)") {
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial base(cfg);

    const double theta = 0.4123; // valeur arbitraire, non triviale (pas 0, pi/2, ...)
    OrientedMaterial oriented_theta(base, theta);
    OrientedMaterial oriented_theta_plus_pi(base, theta + M_PI);

    const Eigen::Vector2d P(0.6, -0.35);
    const Eigen::Matrix2d strain = (Eigen::Matrix2d() << 0.01, 0.004, 0.004, -0.007).finished();
    const Eigen::Vector2d E(0.08, 0.02);
    const Eigen::Matrix2d grad_P = (Eigen::Matrix2d() << 0.015, -0.01, 0.02, 0.01).finished();
    const double v = 0.6;
    const double penalite = v * v + cfg.eta_k;

    CHECK_NEAR(oriented_theta.W_energy(P, strain), oriented_theta_plus_pi.W_energy(P, strain), kTol);
    CHECK_NEAR(oriented_theta.chi_energy(P), oriented_theta_plus_pi.chi_energy(P), kTol);
    CHECK_NEAR(oriented_theta.U_energy(grad_P), oriented_theta_plus_pi.U_energy(grad_P), kTol);

    Eigen::Vector3d sigma_t = oriented_theta.compute_sigma_0(P);
    Eigen::Vector3d sigma_p = oriented_theta_plus_pi.compute_sigma_0(P);
    CHECK_NEAR(sigma_t(0), sigma_p(0), kTol);
    CHECK_NEAR(sigma_t(1), sigma_p(1), kTol);
    CHECK_NEAR(sigma_t(2), sigma_p(2), kTol);

    Eigen::Matrix3d C_t = oriented_theta.get_elastic_matrix();
    Eigen::Matrix3d C_p = oriented_theta_plus_pi.get_elastic_matrix();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK_NEAR(C_t(i,j), C_p(i,j), kTol);

    for (bool is_imp : {true, false}) {
        CHECK_NEAR(oriented_theta.compute_H_drive(grad_P, P, strain, E, is_imp),
                   oriented_theta_plus_pi.compute_H_drive(grad_P, P, strain, E, is_imp), kTol);

        auto GL_t = oriented_theta.compute_GL_terms(P, strain, E, penalite, is_imp);
        auto GL_p = oriented_theta_plus_pi.compute_GL_terms(P, strain, E, penalite, is_imp);
        CHECK_NEAR(GL_t.force_px, GL_p.force_px, kTol);
        CHECK_NEAR(GL_t.force_py, GL_p.force_py, kTol);
        CHECK_NEAR(GL_t.J_11, GL_p.J_11, kTol);
        CHECK_NEAR(GL_t.J_12, GL_p.J_12, kTol);
        CHECK_NEAR(GL_t.J_22, GL_p.J_22, kTol);

        Eigen::Vector2d Peff_t = oriented_theta.compute_effective_polarization(P, v, cfg.eta_k, is_imp);
        Eigen::Vector2d Peff_p = oriented_theta_plus_pi.compute_effective_polarization(P, v, cfg.eta_k, is_imp);
        CHECK_NEAR(Peff_t(0), Peff_p(0), kTol);
        CHECK_NEAR(Peff_t(1), Peff_p(1), kTol);
    }
}

TEST_CASE("OrientedMaterial : la rotation change bien le resultat a theta != 0 (mod pi/2)") {
    // Contrepoint des deux tests precedents : s'assurer qu'on n'a pas ecrit
    // une rotation qui compile mais ne fait rien (ce qui ferait passer les
    // tests ci-dessus de facon triviale et inutile).
    //
    // IMPORTANT : W_energy (= W_mec + W_elas) N'EST PAS le bon candidat ici.
    // Avec les coefficients du papier (Table 1), b1 - b2 = 2*b3 EXACTEMENT
    // (1.4282 - (-0.185) = 1.6132 = 2*0.8066), qui est precisement la
    // condition d'isotropie 2D du couplage electrostrictif P<->eps. Le
    // terme W_energy est donc mathematiquement invariant par rotation,
    // quels que soient P, strain et theta - ce n'est pas un artefact,
    // c'est une propriete reelle des constantes de couplage BaTiO3 du
    // papier (l'anisotropie du cristal tetragonal vient du terme de
    // Landau-Devonshire chi, pas du couplage electrostrictif W_mec).
    //
    // On teste donc chi_energy (alpha_11 != alpha_12 dans le papier),
    // qui est genuinement anisotrope et doit varier avec theta.
    MaterialConfig cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial base(cfg);
    OrientedMaterial oriented_0(base, 0.0);
    OrientedMaterial oriented_30(base, M_PI / 6.0);

    const Eigen::Vector2d P(0.8, 0.3);

    double chi0 = oriented_0.chi_energy(P);
    double chi30 = oriented_30.chi_energy(P);
    CHECK_TRUE(std::abs(chi0 - chi30) > 1e-6);
}