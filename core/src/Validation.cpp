// #include "Validation.h"
// #include "Core/Mesh.h"
// #include "Physics/Polarization.h"
// #include "Physics/Fracture.h"
// #include "Physics/Electrostatics.h"
// #include "Physics/Mechanics.h"  
// #include <functional>
// #include <cmath>
// #include <iostream>
// #include <iomanip>

// Validation::Validation(const Math& math_) : math(math_) {}

// // ==========================================
// // Utilitaires de différences finies
// // ==========================================

// double Validation::finite_diff_dP(const std::function<double(const Eigen::Vector2d&)>& f,
//                                    const Eigen::Vector2d& P, int component) const {
//     Eigen::Vector2d P_plus = P, P_minus = P;
//     P_plus(component)  += FD_EPSILON;
//     P_minus(component) -= FD_EPSILON;
//     return (f(P_plus) - f(P_minus)) / (2.0 * FD_EPSILON);
// }

// double Validation::finite_diff_d2P(const std::function<double(const Eigen::Vector2d&)>& f,
//                                     const Eigen::Vector2d& P) const {
//     const double h = FD_EPSILON_2ND;
//     Eigen::Vector2d Pxy_pp = P, Pxy_pm = P, Pxy_mp = P, Pxy_mm = P;
//     Pxy_pp(0) += h; Pxy_pp(1) += h;
//     Pxy_pm(0) += h; Pxy_pm(1) -= h;
//     Pxy_mp(0) -= h; Pxy_mp(1) += h;
//     Pxy_mm(0) -= h; Pxy_mm(1) -= h;

//     return (f(Pxy_pp) - f(Pxy_pm) - f(Pxy_mp) + f(Pxy_mm)) / (4.0 * h * h);
// }

// ValidationResult Validation::make_result(const std::string& name, double analytic, double numeric, double tol) const {
//     double denom = std::max(1e-12, std::abs(analytic));
//     double err = std::abs(analytic - numeric) / denom;
//     return { name, err <= tol, err, tol };
// }

// // ==========================================
// // NIVEAU 0 : Tests Mathématiques et Jacobiens
// // ==========================================

// ValidationResult Validation::test_dchi_dp1() {
//     Eigen::Vector2d P(0.6, 0.3);
//     auto f = [&](const Eigen::Vector2d& p) { return math.chi_energy(p); };
//     double numeric = finite_diff_dP(f, P, 0);

//     Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);
    
//     return make_result("[Niv 0] dchi_dp1 (Landau-Devonshire, strain=0)", GL.force_px, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_dchi_dp2() {
//     Eigen::Vector2d P(0.6, 0.3);
//     auto f = [&](const Eigen::Vector2d& p) { return math.chi_energy(p); };
//     double numeric = finite_diff_dP(f, P, 1);

//     Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);
    
//     return make_result("[Niv 0] dchi_dp2 (Landau-Devonshire, strain=0)", GL.force_py, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_d2chi_dp1dp2_symmetry() {
//     Eigen::Vector2d P(0.5, 0.4);
//     auto f = [&](const Eigen::Vector2d& p) { return math.chi_energy(p); };
//     double numeric = finite_diff_d2P(f, P);

//     Eigen::Matrix2d zero_strain = Eigen::Matrix2d::Zero();
//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, zero_strain, zero_E, 1.0, false);
    
//     return make_result("[Niv 0] d2chi_dp1dp2 (Jacobien croisé, strain=0)", GL.J_12, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_dW_dp1() {
//     Eigen::Vector2d P(0.5, 0.3);
//     Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;

//     auto f = [&](const Eigen::Vector2d& p) { return math.W_energy(p, strain); };
//     double numeric = finite_diff_dP(f, P, 0);

//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, strain, zero_E, 1.0, false);
//     auto GL_chi_only = math.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
//     double analytic_dW_dp1 = GL.force_px - GL_chi_only.force_px;

//     return make_result("[Niv 0] dW_dp1 (Couplage électrostrictif P <-> Eps)", analytic_dW_dp1, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_dW_dp2() {
//     Eigen::Vector2d P(0.5, 0.3);
//     Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;

//     auto f = [&](const Eigen::Vector2d& p) { return math.W_energy(p, strain); };
//     double numeric = finite_diff_dP(f, P, 1);

//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, strain, zero_E, 1.0, false);
//     auto GL_chi_only = math.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
//     double analytic_dW_dp2 = GL.force_py - GL_chi_only.force_py;

//     return make_result("[Niv 0] dW_dp2 (Couplage électrostrictif P <-> Eps)", analytic_dW_dp2, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_d2W_dp1dp2_symmetry() {
//     Eigen::Vector2d P(0.4, 0.4);
//     Eigen::Matrix2d strain; strain << 0.01, 0.006, 0.006, -0.004;

//     auto f = [&](const Eigen::Vector2d& p) { return math.W_energy(p, strain); };
//     double numeric = finite_diff_d2P(f, P);

//     Eigen::Vector2d zero_E = Eigen::Vector2d::Zero();
//     auto GL = math.compute_GL_terms(P, strain, zero_E, 1.0, false);
//     auto GL_chi_only = math.compute_GL_terms(P, Eigen::Matrix2d::Zero(), zero_E, 1.0, false);
//     double analytic = GL.J_12 - GL_chi_only.J_12;

//     return make_result("[Niv 0] d2W_dp1dp2 (Jacobien croisé électrostrictif)", analytic, numeric, DEFAULT_TOL);
// }

// ValidationResult Validation::test_sigma0_conjugate_to_W() {
//     Eigen::Vector2d P(0.5, 0.3);
//     Eigen::Vector2d P_zero(0.0, 0.0);
//     double eps11 = 0.02, eps22 = -0.015, eps12 = 0.008;

//     auto W_mec_only = [&](double e11, double e22, double e12) {
//         Eigen::Matrix2d strain; strain << e11, e12, e12, e22;
//         return math.W_energy(P, strain) - math.W_energy(P_zero, strain);
//     };

//     double dW_deps11_num = (W_mec_only(eps11 + FD_EPSILON, eps22, eps12) - W_mec_only(eps11 - FD_EPSILON, eps22, eps12)) / (2.0 * FD_EPSILON);
//     double dW_deps22_num = (W_mec_only(eps11, eps22 + FD_EPSILON, eps12) - W_mec_only(eps11, eps22 - FD_EPSILON, eps12)) / (2.0 * FD_EPSILON);
//     double dW_deps12_num = (W_mec_only(eps11, eps22, eps12 + FD_EPSILON) - W_mec_only(eps11, eps22, eps12 - FD_EPSILON)) / (2.0 * FD_EPSILON);

//     Eigen::Vector3d sigma_0 = math.compute_sigma_0(P);
//     double max_err = std::max({std::abs(sigma_0(0) - dW_deps11_num), std::abs(sigma_0(1) - dW_deps22_num), std::abs(sigma_0(2) - 0.5 * dW_deps12_num)});
//     double scale = std::max({std::abs(dW_deps11_num), std::abs(dW_deps22_num), std::abs(dW_deps12_num), 1e-12});

//     return { "[Niv 0] sigma_0 conjugué à W_mec (Test thermodynamique)", (max_err / scale) <= DEFAULT_TOL, max_err / scale, DEFAULT_TOL };
// }

// ValidationResult Validation::test_H_drive_positivity() {
//     Eigen::Matrix2d grad_P; grad_P << 0.1, -0.05, 0.02, 0.08;
//     Eigen::Vector2d P(0.7, 0.2);
//     Eigen::Matrix2d strain; strain << 0.01, 0.002, 0.002, -0.005;
//     Eigen::Vector2d E(0.01, -0.02);

//     double H_perm = math.compute_H_drive(grad_P, P, strain, E, false);
//     double H_imperm = math.compute_H_drive(grad_P, P, strain, E, true);

//     bool ok = (H_perm >= 0.0) && (H_imperm >= 0.0);
//     return { "[Niv 0] H_drive strictement positif (Critère d'irréversibilité)", ok, ok ? 0.0 : 1.0, 0.0 };
// }

// ValidationResult Validation::test_H_drive_zero_at_equilibrium() {
//     Eigen::Vector2d P(1.0, 0.0);
//     double H = math.compute_H_drive(Eigen::Matrix2d::Zero(), P, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), false);
//     bool ok = std::abs(H) < 1e-10;
//     return { "[Niv 0] H_drive nul à l'équilibre homogène", ok, std::abs(H), 1e-10 };
// }

// // ==========================================
// // NIVEAU 1 : Dynamique Locale 0D (Polarisation)
// // ==========================================

// Eigen::Vector2d Validation::integrate_polarization_point(Eigen::Vector2d P, const Eigen::Matrix2d& strain,
//                                                            const Eigen::Vector2d& E, int n_steps) const {
//     auto enthalpy = [&](const Eigen::Vector2d& p_val) {
//         return math.W_energy(p_val, strain) + math.chi_energy(p_val) - p_val(0)*E(0) - p_val(1)*E(1);
//     };

//     double alpha = 0.05; // Pas de départ beaucoup plus petit

//     for (int step = 0; step < n_steps; ++step) {
//         auto GL = math.compute_GL_terms(P, strain, E, 1.0, false);
//         Eigen::Vector2d grad(GL.force_px, GL.force_py);

//         if (grad.norm() < 1e-8) break;

//         double current_h = enthalpy(P);
//         bool step_accepted = false;

//         for (int ls = 0; ls < 20; ++ls) {
//             Eigen::Vector2d step_vec = alpha * grad;
            
//             // SÉCURITÉ CRITIQUE : On bloque la taille du pas pour éviter l'explosion de P^8
//             if (step_vec.norm() > 0.1) step_vec = step_vec.normalized() * 0.1;
            
//             Eigen::Vector2d P_new = P - step_vec;
//             double new_h = enthalpy(P_new);
            
//             if (new_h < current_h) {
//                 P = P_new;
//                 step_accepted = true;
//                 break;
//             }
//             alpha *= 0.5;
//         }

//         if (!step_accepted) break;
//         alpha *= 1.2;
//         if (alpha > 1.0) alpha = 1.0;
//     }
//     return P;
// }

// ValidationResult Validation::test_relaxation_to_minimum() {
//     Eigen::Vector2d P0(0.01, 0.001); 
//     Eigen::Vector2d P_final = integrate_polarization_point(P0, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);

//     auto GL_final = math.compute_GL_terms(P_final, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 1.0, false);
//     double force_norm = std::sqrt(GL_final.force_px * GL_final.force_px + GL_final.force_py * GL_final.force_py);
//     double max_stiffness = std::max({std::abs(GL_final.J_11), std::abs(GL_final.J_22), 1e-6});

//     return make_result("[Niv 1] Dynamique 0D : Relaxation vers un équilibre stable", 0.0, force_norm / max_stiffness, 1e-4);
// }

// ValidationResult Validation::test_stationary_equilibrium() {
//     Eigen::Vector2d P_min = integrate_polarization_point(Eigen::Vector2d(0.01, 0.0), Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);
    
//     // Perturbation du point stable
//     Eigen::Vector2d P_perturbed = P_min;
//     P_perturbed(0) *= 0.8; 
//     P_perturbed(1) += 0.05;

//     // L'algorithme doit retrouver EXACTEMENT le même point
//     Eigen::Vector2d P_recovered = integrate_polarization_point(P_perturbed, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);
//     double drift = (P_recovered - P_min).norm() / std::max(1e-6, P_min.norm());
    
//     return make_result("[Niv 1] Dynamique 0D : Stabilité locale du puits de potentiel", 0.0, drift, 1e-4);
// }

// ValidationResult Validation::test_switching_under_stress() {
//     Eigen::Vector2d P = integrate_polarization_point(Eigen::Vector2d(0.01, 0.001), Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero(), 20000);
//     bool switched = false;

//     // Rampe de traction selon Y pour forcer le basculement ferroélastique
//     for (int i = 1; i <= 200; ++i) {
//         Eigen::Matrix2d strain; strain << 0.0, 0.0, 0.0, i * 0.005; 
//         P = integrate_polarization_point(P, strain, Eigen::Vector2d::Zero(), 2000);

//         if (std::abs(P(1)) > 1.5 * std::abs(P(0))) {
//             switched = true;
//             break;
//         }
//     }

//     return { "[Niv 1] Dynamique 0D : Ferroelastic Switching 90° sous traction", switched, switched ? 0.0 : 1.0, 0.0 };
// }

// ValidationResult Validation::test_mechanics_patch_energy() {
//     double a = 0.01, b = -0.005;
//     Eigen::Matrix3d C = math.get_elastic_matrix();
//     Eigen::Vector3d eps_voigt(a, b, 0.0);
//     double W_analytic = 0.5 * eps_voigt.transpose() * C * eps_voigt;

//     Eigen::Matrix2d strain; strain << a, 0.0, 0.0, b;
//     double W_numeric = math.W_energy(Eigen::Vector2d::Zero(), strain);

//     return make_result("[Niv 1] Mécanique : Patch test W_elas avec matrice C", W_analytic, W_numeric, 1e-8);
// }

// // ==========================================
// // NIVEAU 2 : Couplage de la Fracture (Champ de phase v)
// // ==========================================

// ValidationResult Validation::test_fracture_penalty_GL_terms() {
//     Eigen::Vector2d P(0.5, 0.3);
//     Eigen::Matrix2d strain; strain << 0.01, 0.005, 0.005, -0.008;
//     Eigen::Vector2d E = Eigen::Vector2d::Zero();
    
//     auto GL_intact = math.compute_GL_terms(P, strain, E, 1.0, false);
//     auto GL_broken = math.compute_GL_terms(P, strain, E, math.eta_k, false);
//     auto GL_chi_only = math.compute_GL_terms(P, Eigen::Matrix2d::Zero(), E, 0.0, false);
    
//     double J11_mec_intact = GL_intact.J_11 - GL_chi_only.J_11;
//     double J11_mec_broken = GL_broken.J_11 - GL_chi_only.J_11;
    
//     // Correction du std::abs() implémentée
//     double ratio = std::abs(J11_mec_broken) / std::max(1e-12, std::abs(J11_mec_intact));
    
//     return make_result("[Niv 2] Fracture : Contrainte libérée (Traction-Free) pour v=0", math.eta_k, ratio, 1e-10);
// }

// ValidationResult Validation::test_crack_permeability_D_field() {
//     Eigen::Vector2d P(0.0, 0.5);
//     Eigen::Vector2d E(0.0, 0.01);
    
//     Eigen::Vector2d D_perm = math.compute_effective_permittivity(0.0, math.eta_k, false) * E 
//                            + math.compute_effective_polarization(P, 0.0, math.eta_k, false);
                           
//     Eigen::Vector2d D_imperm = math.compute_effective_permittivity(0.0, math.eta_k, true) * E 
//                              + math.compute_effective_polarization(P, 0.0, math.eta_k, true);
    
//     double expected_perm = (math.eps0 * E + P).norm();
//     double expected_imperm = math.eta_k * expected_perm;
    
//     double err_perm = std::abs(D_perm.norm() - expected_perm);
//     double err_imperm = std::abs(D_imperm.norm() - expected_imperm);
    
//     bool passed = (err_perm < 1e-10) && (err_imperm < 1e-10);
//     return { "[Niv 2] Fracture : Champ D (Condition Perméable vs Imperméable)", passed, std::max(err_perm, err_imperm), 1e-8 };
// }

// ValidationResult Validation::test_dh_dv_driving_force() {
//     Eigen::Matrix2d grad_P; grad_P << 0.1, -0.05, 0.02, 0.08;
//     Eigen::Vector2d P(0.4, 0.4);
//     Eigen::Matrix2d strain; strain << 0.01, 0.0, 0.0, -0.01;
//     Eigen::Vector2d E(0.01, 0.01);
//     double v_test = 0.5;
    
//     auto enthalpy = [&](double v_val, bool imperm) {
//         double penalty = v_val * v_val + math.eta_k;
//         double U = math.U_energy(grad_P);
//         double W = math.W_energy(P, strain);
//         double chi = math.chi_energy(P);
//         double W_elec_base = -P.dot(E);
//         double W_elec_imperm = W_elec_base - 0.5 * math.eps0 * E.squaredNorm();
        
//         return imperm ? penalty * (U + W + W_elec_imperm) + chi 
//                       : penalty * (U + W) + chi + W_elec_base - 0.5 * math.eps0 * E.squaredNorm(); 
//     };
    
//     double dh_dv_perm_num = (enthalpy(v_test + FD_EPSILON, false) - enthalpy(v_test - FD_EPSILON, false)) / (2.0 * FD_EPSILON);
//     double dh_dv_imperm_num = (enthalpy(v_test + FD_EPSILON, true) - enthalpy(v_test - FD_EPSILON, true)) / (2.0 * FD_EPSILON);
    
//     double dh_dv_perm_analy = 2.0 * v_test * math.compute_H_drive(grad_P, P, strain, E, false);
//     double dh_dv_imperm_analy = 2.0 * v_test * math.compute_H_drive(grad_P, P, strain, E, true);
    
//     double err_perm = std::abs(dh_dv_perm_num - dh_dv_perm_analy) / std::max(1e-12, std::abs(dh_dv_perm_num));
//     double err_imperm = std::abs(dh_dv_imperm_num - dh_dv_imperm_analy) / std::max(1e-12, std::abs(dh_dv_imperm_num));
    
//     bool passed = (err_perm < DEFAULT_TOL) && (err_imperm < DEFAULT_TOL);
//     return { "[Niv 2] Fracture : Force thermodynamique dh/dv (H_drive)", passed, std::max(err_perm, err_imperm), DEFAULT_TOL };
// }

// // ==========================================
// // NIVEAU 3 : Validation Spatiale (Éléments Finis)
// // ==========================================

// ValidationResult Validation::test_domain_wall_profile(const Datafile& config, const Mesh& mesh, Polarization& pol, Fracture& frac, Mechanics& mec, Electrostatics& elec) {
//     int n_nodes = (int)mesh.get_nodes().size();
//     if (n_nodes == 0) return { "[Niv 3] Spatial : Profil analytique de paroi 180°", false, 0.0, 0.0 };

//     frac.set_v(Eigen::VectorXd::Ones(n_nodes));
//     elec.set_phi(Eigen::VectorXd::Zero(n_nodes));
//     mec.set_ux(Eigen::VectorXd::Zero(n_nodes)); mec.set_uy(Eigen::VectorXd::Zero(n_nodes));

//     // Utiliser la VRAIE dimension du maillage
//     double max_x = 1e-9;
//     for (int i = 0; i < n_nodes; ++i) max_x = std::max(max_x, mesh.get_nodes()[i].x);
//     double center_x = max_x / 2.0;
//     double delta_th = std::sqrt(config.material.a0 / std::abs(config.material.alpha_1)); 

//     Eigen::VectorXd Px_init(n_nodes);
//     for(int i = 0; i < n_nodes; ++i) {
//         Px_init(i) = -std::tanh((mesh.get_nodes()[i].x - center_x) / delta_th);
//     }
//     pol.set_Px(Px_init);
//     pol.set_Py(Eigen::VectorXd::Zero(n_nodes));

//     pol.save_previous_state();
//     pol.save_previous_iteration();
//     pol.update_P(0.0, frac, mec, elec, math);

//     double max_error = 0.0;
//     for (int i = 0; i < n_nodes; ++i) {
//         double Px_analytic = -std::tanh((mesh.get_nodes()[i].x - center_x) / delta_th);
//         max_error = std::max(max_error, std::abs(pol.get_Px()(i) - Px_analytic));
//     }
    
//     // Tolérance adaptée à la discrétisation d'un maillage
//     return make_result("[Niv 3] Spatial : Profil analytique de paroi 180°", 1.0, 1.0 + max_error, 1.5e-1);
// }

// ValidationResult Validation::test_phase_field_crack_profile(const Datafile& config, const Mesh& mesh, Fracture& frac, Polarization& pol, Mechanics& mec, Electrostatics& elec) {
//     int n_nodes = (int)mesh.get_nodes().size();
//     if (n_nodes == 0) return { "[Niv 3] Spatial : Profil fissure", false, 0.0, 0.0 };

//     pol.set_Px(Eigen::VectorXd::Zero(n_nodes)); pol.set_Py(Eigen::VectorXd::Zero(n_nodes));
//     elec.set_phi(Eigen::VectorXd::Zero(n_nodes));
//     mec.set_ux(Eigen::VectorXd::Zero(n_nodes)); mec.set_uy(Eigen::VectorXd::Zero(n_nodes));

//     double max_x = 1e-9;
//     for (int i = 0; i < n_nodes; ++i) max_x = std::max(max_x, mesh.get_nodes()[i].x);
//     double center_x = max_x / 2.0;
//     double kappa = config.material.kappa;
    
//     Eigen::VectorXd v_init(n_nodes);
//     for (int i = 0; i < n_nodes; ++i) {
//         v_init(i) = 1.0 - std::exp(-std::abs(mesh.get_nodes()[i].x - center_x) / (2.0 * kappa));
//     }
//     frac.set_v(v_init);

//     frac.save_previous_state();
//     frac.save_previous_iteration();
//     frac.update_v(config.simulation.dt, pol, mec, elec, math);

//     double max_error = 0.0;
//     for (int i = 0; i < n_nodes; ++i) {
//         double v_analytic = 1.0 - std::exp(-std::abs(mesh.get_nodes()[i].x - center_x) / (2.0 * kappa));
//         max_error = std::max(max_error, std::abs(frac.get_v()(i) - v_analytic));
//     }
    
//     return make_result("[Niv 3] Spatial : Profil fissure", 1.0, 1.0 + max_error, 5e-2);
// }

// ValidationResult Validation::test_electrostatics_patch(const Datafile& config, const Mesh& mesh, Electrostatics& elec, Polarization& pol, Fracture& frac, Mechanics& mec) {
//     int n_nodes = (int)mesh.get_nodes().size();
//     if (n_nodes == 0) return { "[Niv 3] Spatial : Patch test Électrostatique", false, 0.0, 0.0 };
    
//     pol.set_Px(Eigen::VectorXd::Zero(n_nodes)); pol.set_Py(Eigen::VectorXd::Zero(n_nodes));
//     frac.set_v(Eigen::VectorXd::Ones(n_nodes)); 
//     mec.set_ux(Eigen::VectorXd::Zero(n_nodes)); mec.set_uy(Eigen::VectorXd::Zero(n_nodes));

//     // Utilisation de la VRAIE dimension max_x pour calculer la pente exacte du Patch Test
//     double max_x = 1e-9;
//     for (int i = 0; i < n_nodes; ++i) max_x = std::max(max_x, mesh.get_nodes()[i].x);
    
//     Eigen::VectorXd phi_init(n_nodes);
//     for (int i = 0; i < n_nodes; ++i) {
//         phi_init(i) = (mesh.get_nodes()[i].x / max_x) * 100.0;
//     }
//     elec.set_phi(phi_init);

//     elec.save_previous_state();
//     elec.save_previous_iteration();
//     elec.update_phi(0.0, pol, frac, math);

//     double max_error = 0.0;
//     for (int i = 0; i < n_nodes; ++i) {
//         double phi_analytic = (mesh.get_nodes()[i].x / max_x) * 100.0;
//         max_error = std::max(max_error, std::abs(elec.get_phi()(i) - phi_analytic));
//     }
    
//     return make_result("[Niv 3] Spatial : Patch test Électrostatique", 1.0, 1.0 + max_error, 1e-4);
// }

// // ==========================================
// // Orchestration et Affichage
// // ==========================================

// std::vector<ValidationResult> Validation::run_all() {
//     return {
//         // Niveau 0
//         test_dchi_dp1(), test_dchi_dp2(), test_d2chi_dp1dp2_symmetry(),
//         test_dW_dp1(), test_dW_dp2(), test_d2W_dp1dp2_symmetry(),
//         test_sigma0_conjugate_to_W(), test_H_drive_positivity(), test_H_drive_zero_at_equilibrium(),
//         // Niveau 1
//         test_relaxation_to_minimum(), test_stationary_equilibrium(), test_switching_under_stress(),
//         test_mechanics_patch_energy(),
//         // Niveau 2
//         test_fracture_penalty_GL_terms(), test_crack_permeability_D_field(), test_dh_dv_driving_force()
//     };
// }

// bool Validation::all_passed(const std::vector<ValidationResult>& results) {
//     for (const auto& r : results) { if (!r.passed) return false; }
//     return true;
// }

// void Validation::print_report(const std::vector<ValidationResult>& results) const {
//     std::cout << "\n========================================================================\n";
//     std::cout << "           RAPPORT DE VALIDATION - Physique et Mathématiques\n";
//     std::cout << "========================================================================\n\n";

//     size_t n_pass = 0;
//     std::string current_level = "";

//     for (const auto& r : results) {
//         // Esthétique : Sépare les blocs de niveaux dans l'affichage
//         std::string level = r.name.substr(0, 7); // Extrait "[Niv X]"
//         if (level != current_level && r.name.rfind("[Niv", 0) == 0) {
//             std::cout << "------------------------------------------------------------------------\n";
//             current_level = level;
//         }

//         std::cout << (r.passed ? "\033[1;32m[PASS]\033[0m " : "\033[1;31m[FAIL]\033[0m ")
//                   << std::left << std::setw(65) << r.name
//                   << " err=" << std::scientific << std::setprecision(2) << r.error
//                   << " (tol=" << r.tolerance << ")\n";
        
//         if (r.passed) n_pass++;
//     }

//     std::cout << "\n========================================================================\n";
//     if (n_pass == results.size()) {
//          std::cout << " \033[1;32mSUCCÈS TOTAL : " << n_pass << " / " << results.size() << " tests réussis.\033[0m\n";
//     } else {
//          std::cout << " \033[1;31mÉCHEC : " << n_pass << " / " << results.size() << " tests réussis. Vérifiez les erreurs ci-dessus.\033[0m\n";
//     }
//     std::cout << "========================================================================\n\n";
// }