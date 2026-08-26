#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Physics/include/Modules/Electrostatics.h"
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
// Tests d'assemblage element pour ElectrostaticsEquation (module Physics).
// Meme convention que test_mechanics_equation.cpp : element Q4 unite carre.
// K_local ici est 4x4 (un seul dof scalaire phi par noeud).
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

// Ecrit un TOML minimal sur disque (meme pattern que test_datafile.cpp) et
// retourne un Datafile parse. mode "PERMEABLE" par defaut (comme le papier,
// Section 3.2.1).
Datafile make_datafile(const std::string& path, const std::string& fracture_mode = "PERMEABLE") {
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n"
        << "[mesh]\nL_x=1.0\nL_y=1.0\nn_x=1\nn_y=1\nelement_type=\"Q4\"\n"
        << "[material]\n[crystal]\nnum_grains=1\n"
        << "[fracture]\nmode=\"" << fracture_mode << "\"\n";
    out.close();
    Datafile config(path);
    std::remove(path.c_str());
    return config;
}

} // namespace

TEST_CASE("ElectrostaticsEquation : K_local est symetrique (element Q4 unite)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_electro_symm.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;

    ElectrostaticsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    CHECK_TRUE((K - K.transpose()).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("ElectrostaticsEquation : potentiel constant -> K*phi nul (patch test de Laplace)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_electro_patch.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;

    ElectrostaticsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    // phi uniforme -> gradient nul -> K*phi = 0, quel que soit eps_eff.
    Eigen::VectorXd phi_uniform = Eigen::VectorXd::Constant(4, 0.75);
    Eigen::VectorXd residual = K * phi_uniform;
    CHECK_TRUE(residual.cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("ElectrostaticsEquation : mode PERMEABLE ignore l'endommagement v dans la permittivite") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_electro_permeable.toml", "PERMEABLE");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    // v_prev = 0 partout : materiau totalement rompu. En PERMEABLE, la
    // fissure ne doit pas isoler electriquement le milieu (cf papier,
    // Section 2.3) -> K_local doit rester identique au cas v=1 (materiau
    // sain), a l'oppose du comportement attendu en IMPERMEABLE.
    Eigen::VectorXd v_broken = Eigen::VectorXd::Zero(4);
    state.fields["v_prev"] = &v_broken;

    ElectrostaticsEquation eq;
    Eigen::MatrixXd K_broken = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F_broken = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K_broken, F_broken, global_dofs);

    Eigen::VectorXd v_sound = Eigen::VectorXd::Ones(4);
    state.fields["v_prev"] = &v_sound;
    Eigen::MatrixXd K_sound = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F_sound = Eigen::VectorXd::Zero(4);
    eq.compute_element_matrices(ctx, state, K_sound, F_sound, global_dofs);

    CHECK_TRUE((K_broken - K_sound).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("ElectrostaticsEquation : mode IMPERMEABLE degrade bien la permittivite avec v") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_electro_impermeable.toml", "IMPERMEABLE");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    Eigen::VectorXd v_broken = Eigen::VectorXd::Zero(4);
    state.fields["v_prev"] = &v_broken;

    ElectrostaticsEquation eq;
    Eigen::MatrixXd K_broken = Eigen::MatrixXd::Zero(4, 4);
    Eigen::VectorXd F_broken = Eigen::VectorXd::Zero(4);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K_broken, F_broken, global_dofs);

    // En IMPERMEABLE, avec v=0 partout, eps_eff = eps0 * eta_k (quasi nul)
    // -> K_local doit etre quasi nulle (fissure isolante electriquement).
    CHECK_TRUE(K_broken.cwiseAbs().maxCoeff() < 1e-4);
}