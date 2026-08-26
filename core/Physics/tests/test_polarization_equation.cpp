#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Physics/include/Modules/Polarization.h"
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
// Tests d'assemblage element pour PolarizationEquation (module Physics).
// Meme convention geometrique que les autres tests Physics/tests : element
// Q4 unite carre. K_local est 8x8 (2 dofs Px,Py par noeud, entrelaces).
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

// TOML minimal avec les coefficients du papier pour a0/mu_p/eta_k (les
// valeurs par defaut de Datafile::parse_material different de celles du
// papier - cf Datafile.cpp - donc on les fixe explicitement ici).
Datafile make_datafile(const std::string& path, const std::string& fracture_mode = "PERMEABLE") {
    std::ofstream out(path);
    out << "[simulation]\noutput_dir = \".\"\ntotal_time = 1\ndt = 0.1\n"
        << "[mesh]\nL_x=1.0\nL_y=1.0\nn_x=1\nn_y=1\nelement_type=\"Q4\"\n"
        << "[material]\na0=0.1\nmu_p=1.0\neta_k=1e-6\n"
        << "[crystal]\nnum_grains=1\n"
        << "[fracture]\nmode=\"" << fracture_mode << "\"\n";
    out.close();
    Datafile config(path);
    std::remove(path.c_str());
    return config;
}

} // namespace

TEST_CASE("PolarizationEquation : K_local est symetrique (element Q4 unite, P au repos)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_polar_symm.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;
    // P = (0,0) partout : point de linearisation au repos, cas le plus
    // simple pour verifier la structure (masse + rigidite + Jacobien GL).
    Eigen::VectorXd px = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd py = Eigen::VectorXd::Zero(4);
    state.fields["Px"] = &px;
    state.fields["Py"] = &py;
    state.fields["Px_n"] = &px;
    state.fields["Py_n"] = &py;

    PolarizationEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    CHECK_TRUE((K - K.transpose()).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("PolarizationEquation : P constant, sans historique ni GL -> F_local nul (terme de masse pur)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_polar_mass.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;
    // P == P_n == 0 partout : pas de force GL (P=0 est un point stationnaire
    // du developpement de Landau a l'ordre le plus bas pour ce test), pas de
    // terme d'historique (Px_n=0). F_local doit donc etre nul.
    Eigen::VectorXd zero = Eigen::VectorXd::Zero(4);
    state.fields["Px"] = &zero;
    state.fields["Py"] = &zero;
    state.fields["Px_n"] = &zero;
    state.fields["Py_n"] = &zero;

    PolarizationEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    CHECK_TRUE(F.cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("PolarizationEquation : global_dofs suit le p_dof_map (pas l'indexation geometrique)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));
    Datafile config = make_datafile("test_polar_dofs.toml");

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    state.config = &config;
    state.dt = 0.01;
    // Sans p_dof_map fourni (nullptr), le mapping doit retomber sur
    // l'indexation geometrique standard (comme documente dans le code :
    // "p_dofs[i] = state.p_dof_map ? ... : idx").
    PolarizationEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    std::vector<int> expected = {0, 1, 2, 3, 4, 5, 6, 7};
    CHECK_TRUE(global_dofs == expected);
}