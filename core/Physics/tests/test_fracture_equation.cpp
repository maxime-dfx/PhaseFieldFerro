#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Physics/include/Modules/Fracture.h"
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Models/Ferroelectric.h"
#include "Mesh/include/Mesh.h"
#include "IO/include/Datafile.h"
#include <Eigen/Dense>
#include <fstream>
#include <cstdio>
#include <memory>

// =========================================================================
// Tests d'assemblage element pour FractureEquation (module Physics).
// Meme convention geometrique que les autres tests Physics/tests : element
// Q4 unite carre. K_local/F_local sont pour le champ scalaire v (endommagement).
// =========================================================================

namespace {

Mesh make_single_q4_mesh() {
    std::vector<Node> nodes = {
        {0.0, 0.0, 0}, {1.0, 0.0, 0}, {1.0, 1.0, 0}, {0.0, 1.0, 0}
    };
    Element e;
    e.num_nodes = 4;
    e.offset = 0;
    std::vector<int> connectivity = {0, 1, 2, 3};
    std::vector<Element> elems = {e};
    return Mesh(1.0, 1.0, 1, 1, ElementType::QUAD4, nodes, elems, connectivity);
}

std::array<std::array<double, 2>, 8> unit_q4_coords() {
    std::array<std::array<double, 2>, 8> coords{};
    coords[0] = {0.0, 0.0};
    coords[1] = {1.0, 0.0};
    coords[2] = {1.0, 1.0};
    coords[3] = {0.0, 1.0};
    return coords;
}

Datafile make_datafile(const std::string& path, const std::string& fracture_mode = "PERMEABLE") {
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n"
        << "[mesh]\nL_x=1.0\nL_y=1.0\nn_x=1\nn_y=1\nelement_type=\"Q4\"\n"
        << "[material]\nmu_v=15.0\nkappa=2.0\nGc=4.0\n"
        << "[crystal]\nnum_grains=1\n"
        << "[fracture]\nmode=\"" << fracture_mode << "\"\n";
    out.close();
    Datafile config(path);
    std::remove(path.c_str());
    return config;
}

} // namespace

TEST_CASE("FractureEquation : K_local est symetrique (champs nuls)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_fracture_symm.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;

    FractureEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    CHECK_TRUE((K - K.transpose()).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("FractureEquation : v_n=1 uniforme sans force motrice -> v=1 est bien la solution d'equilibre") {

    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_fracture_equilibrium.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;
    Eigen::VectorXd v_n = Eigen::VectorXd::Ones(4);
    state.fields["v_prev"] = &v_n;

    FractureEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    Eigen::VectorXd v_uniform = Eigen::VectorXd::Ones(4);
    Eigen::VectorXd residual = K * v_uniform - F;
    CHECK_TRUE(residual.cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("FractureEquation : mode PERMEABLE annule le champ electrique dans H_drive") {

    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config_permeable = make_datafile("test_fracture_permeable.toml", "PERMEABLE");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config_permeable;
    state.dt = 0.01;
    Eigen::VectorXd v_n = Eigen::VectorXd::Ones(4);
    // phi non uniforme -> E_gp non nul si pris en compte.
    Eigen::VectorXd phi(4);
    phi << 0.0, 1.0, 1.0, 0.0;
    state.fields["v_prev"] = &v_n;
    state.fields["phi_prev"] = &phi;

    FractureEquation eq;
    Eigen::MatrixXd K_with_phi = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F_with_phi = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K_with_phi, F_with_phi, global_dofs);

    state.fields.erase("phi_prev");
    Eigen::MatrixXd K_no_phi = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F_no_phi = Eigen::VectorXd::Zero(4);
    eq.compute_element_matrices(ctx, state, K_no_phi, F_no_phi, global_dofs);

    CHECK_TRUE((K_with_phi - K_no_phi).cwiseAbs().maxCoeff() < 1e-12);
    CHECK_TRUE((F_with_phi - F_no_phi).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("FractureEquation : garde-fou mass_coeff_floor empeche un coefficient diagonal negatif ou nul") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_fracture_floor.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;
    Eigen::VectorXd v_n = Eigen::VectorXd::Ones(4);
    Eigen::VectorXd px = Eigen::VectorXd::Constant(4, 5.0);
    Eigen::VectorXd py = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd ux(4), uy(4);

    ux << 0.0, -0.5, -0.5, 0.0;
    uy << 0.0, 0.0, 0.0, 0.0;
    state.fields["v_prev"] = &v_n;
    state.fields["Px_prev"] = &px;
    state.fields["Py_prev"] = &py;
    state.fields["ux_prev"] = &ux;
    state.fields["uy_prev"] = &uy;

    FractureEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    for (int i = 0; i < 4; ++i) {
        CHECK_TRUE(K(i, i) > 0.0);
    }
}

TEST_CASE("FractureEquation : DIAGNOSTIC - relaxation TDGL de v sur plusieurs iterations Picard (m)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_fracture_relax_ratio.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.1; // dt_relax de production (input/config.toml)

    Eigen::VectorXd v_prev = Eigen::VectorXd::Ones(4); // v^0 = v^(n-1) = 1
    state.fields["v_prev"] = &v_prev;

    // H_drive impose indirectement via une polarisation uniforme + un

    Eigen::VectorXd px = Eigen::VectorXd::Constant(4, 8.0);
    Eigen::VectorXd py = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd ux(4), uy(4);
    ux << 0.0, 0.5, 0.5, 0.0; // eps_11 = +0.5 (traction, pas de compression)
    uy << 0.0, 0.0, 0.0, 0.0;
    state.fields["Px_prev"] = &px;
    state.fields["Py_prev"] = &py;
    state.fields["ux_prev"] = &ux;
    state.fields["uy_prev"] = &uy;

    FractureEquation eq;
    const int n_picard_iters = 30;
    double v_new = 1.0;
    for (int m = 1; m <= n_picard_iters; ++m) {
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
        Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
        std::vector<int> global_dofs;
        eq.compute_element_matrices(ctx, state, K, F, global_dofs);
        Eigen::VectorXd v_solution = K.colPivHouseholderQr().solve(F);
        v_new = v_solution(0); // uniforme par symetrie du cas teste
        v_prev = v_solution;   // save_previous_iteration() : v^(m-1) <- v^m
    }

    double chute_pourcent = (1.0 - v_new) * 100.0;

    Logger::info("[DIAGNOSTIC relax_ratio] v^0=1.0 -> v^", n_picard_iters,
                   "=", v_new, " (chute cumulee de ",
                   chute_pourcent,
                   "% sur ", n_picard_iters,
                   " iterations Picard, dt_relax=0.1, H_drive eleve)");

    CHECK_TRUE(v_new >= 0.0 && v_new <= 1.0);
    CHECK_TRUE(chute_pourcent > 5.0);
    CHECK_TRUE(chute_pourcent > 5.0);
}