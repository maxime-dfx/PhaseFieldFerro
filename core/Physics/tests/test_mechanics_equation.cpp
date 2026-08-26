#include "Tests/include/TestFramework.h"
#include "Tests/include/TestFixtures.h"
#include "Physics/include/Modules/Mechanics.h"
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Materials/Core/MaterialManager.h"
#include "Materials/Models/Ferroelectric.h"
#include "Mesh/include/Mesh.h"
#include <Eigen/Dense>
#include <memory>

// =========================================================================
// Tests d'assemblage element pour MechanicsEquation (module Physics).
// Contrairement aux tests "Validation Niv0/1/2" (Materials/tests), qui
// verifient les formules materiau isolement, on teste ici l'assemblage FEM
// complet (K_local, F_local) tel qu'utilise reellement par SystemAssembler.
//
// Convention element : Q4 unite carre, noeuds (0,0)-(1,0)-(1,1)-(0,1),
// dofs entrelaces [ux0,uy0,ux1,uy1,ux2,uy2,ux3,uy3].
// =========================================================================

namespace {

// Construit un unique element Q4 carre unite + son Mesh associe.
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

} // namespace

TEST_CASE("MechanicsEquation : K_local est symetrique (element Q4 unite, materiau papier)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    SingleCrystalMaterial material(mat_cfg);
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    // Pas de champ Px/Py/v_prev fourni -> P=0, v=0 partout (degradation = eta_k).

    MechanicsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;

    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    CHECK_TRUE((K - K.transpose()).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("MechanicsEquation : mouvement de corps rigide -> force interne K*u nulle") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    // v_prev = 1 partout pour avoir une degradation = 1 (materiau sain, cas
    // le plus discriminant pour un patch test standard de rigidite).
    Eigen::VectorXd v_field = Eigen::VectorXd::Ones(4);
    state.fields["v_prev"] = &v_field;

    MechanicsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    // Translation rigide u = (0.37, -0.21) identique sur les 4 noeuds :
    // aucune deformation, donc K*u doit etre nul (a la precision numerique
    // pres) quelle que soit la matrice elastique C.
    Eigen::VectorXd u_rigid(8);
    for (int i = 0; i < 4; ++i) {
        u_rigid(2*i)   = 0.37;
        u_rigid(2*i+1) = -0.21;
    }
    Eigen::VectorXd f_internal = K * u_rigid;
    CHECK_TRUE(f_internal.cwiseAbs().maxCoeff() < 1e-9);
}

TEST_CASE("MechanicsEquation : polarisation nulle -> pas de contrainte spontanee (F_local == 0)") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;
    Eigen::VectorXd v_field = Eigen::VectorXd::Ones(4);
    Eigen::VectorXd px_field = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd py_field = Eigen::VectorXd::Zero(4);
    state.fields["v_prev"] = &v_field;
    state.fields["Px_prev"] = &px_field;
    state.fields["Py_prev"] = &py_field;

    MechanicsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    // Sans polarisation, sigma_0(P=0) doit etre nul (couplage purement
    // electrostrictif, cf Eq. 4 du papier) -> F_local = -B^T*sigma_0*dV = 0.
    CHECK_TRUE(F.cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("MechanicsEquation : global_dofs correspond bien a l'entrelacement [ux,uy] par noeud") {
    MaterialConfig mat_cfg = TestFixtures::make_paper_material_config();
    MaterialManager registry;
    registry.register_material(0, std::make_unique<SingleCrystalMaterial>(mat_cfg));

    Mesh mesh = make_single_q4_mesh();
    auto coords = unit_q4_coords();
    ElementContext ctx{mesh.get_elements()[0], coords, 0, 0, mesh};

    PhysicsState state;
    state.materials_manager = &registry;

    MechanicsEquation eq;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(8);
    std::vector<int> global_dofs;
    eq.compute_element_matrices(ctx, state, K, F, global_dofs);

    std::vector<int> expected = {0, 1, 2, 3, 4, 5, 6, 7};
    CHECK_TRUE(global_dofs == expected);
}   