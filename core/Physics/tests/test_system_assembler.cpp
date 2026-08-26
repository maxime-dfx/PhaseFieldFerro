#include "Tests/include/TestFramework.h"
#include "Physics/include/Assembly/SystemAssembler.h"
#include "Physics/include/Core/PhysicsConcepts.h"
#include "Mesh/include/Mesh.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

// =========================================================================
// Tests SystemAssembler : verifient l'assemblage global (scatter-add des
// K_local/F_local elementaires) independamment de toute physique reelle,
// via une IElementEquation bouchon dont la contribution est connue a
// l'avance et calculable a la main.
// Maillage : 2 elements Q4 partageant une arete (noeuds 1 et 2 communs).
//   3---2---5
//   |E0 |E1 |
//   0---1---4
// element0 = [0,1,2,3], element1 = [1,4,5,2]. Champ scalaire (1 dof/noeud,
// comme ElectrostaticsEquation/FractureEquation) -> num_dofs_total = 6.
// =========================================================================

namespace {

// Equation bouchon : K_local = Identite(n,n), F_local = Ones(n). Le
// scatter-add resultant est donc calculable a la main : chaque noeud
// partage par plusieurs elements accumule autant de contributions que
// d'elements auxquels il appartient, sur la diagonale de K et dans F ; les
// termes hors-diagonale de K_local etant nuls, K_global reste diagonale.
class ConstantIdentityEquation : public IElementEquation {
public:
    // --- TAGS POUR TRACY OBLIGATOIRES POUR QUE L'ASSEMBLEUR COMPILE ---
    static constexpr const char* ModuleName       = "TestDummy";
    static constexpr uint32_t    Color            = PROFILE_COLOR_SOLVE;
    static constexpr const char* ZoneAssemblySlow = "TestDummy::Assembly_Triplets";
    static constexpr const char* ZoneAssemblyFast = "TestDummy::Assembly_CSR";
    static constexpr const char* ZoneSolveCompute = "TestDummy::Solve_Compute";
    static constexpr const char* ZoneSolveApply   = "TestDummy::Solve_Apply";
    static constexpr const char* ZonePostProcess  = "TestDummy::PostProcessing";

    void compute_element_matrices(const ElementContext& ctx, const PhysicsState&,
                                  Eigen::Ref<Eigen::MatrixXd> K_local, Eigen::Ref<Eigen::VectorXd> F_local,
                                  std::vector<int>& global_dofs) const override
    {
        int n = ctx.elem.num_nodes;
        global_dofs.clear();
        for (int i = 0; i < n; ++i) {
            global_dofs.push_back(ctx.mesh.get_node_index(ctx.elem_idx, i));
            K_local(i, i) += 1.0;
            F_local(i) += 1.0;
        }
    }
};

Mesh make_two_element_mesh() {
    std::vector<Node> nodes = {
        {0.0, 0.0, 0}, {1.0, 0.0, 0}, {1.0, 1.0, 0}, {0.0, 1.0, 0},
        {2.0, 0.0, 0}, {2.0, 1.0, 0}
    };
    Element e0; e0.num_nodes = 4; e0.offset = 0;
    Element e1; e1.num_nodes = 4; e1.offset = 4;
    std::vector<Element> elems = {e0, e1};
    std::vector<int> connectivity = {0, 1, 2, 3,  1, 4, 5, 2};
    
    Mesh m(2.0, 1.0, 2, 1, ElementType::QUAD4, nodes, elems, connectivity);
    m.compute_coloring(); 
    return m;
}

} // namespace

TEST_CASE("SystemAssembler : scatter-add correct sur noeuds partages (premier assemblage, chemin triplets)") {
    Mesh mesh = make_two_element_mesh();
    ConstantIdentityEquation eq;
    PhysicsState state;
    std::vector<FlattenedBC> no_bcs;
    std::vector<long> csr_mapping; // Cache CSR requis par la nouvelle signature
    Eigen::SparseMatrix<double> K_global;
    Eigen::VectorXd F_global;

    // Suppression de l'argument "test"
    SystemAssembler::assemble(mesh.get_num_elements(), mesh.get_num_nodes(), /*is_first_assembly=*/true,
                            K_global, F_global, csr_mapping, eq, no_bcs, state, mesh);

    Eigen::MatrixXd K_dense = Eigen::MatrixXd(K_global);

    // Noeuds 1 et 2 (partages entre element0 et element1) doivent accumuler
    // 2 contributions ; les autres (0,3,4,5) une seule.
    Eigen::VectorXd expected_diag(6);
    expected_diag << 1.0, 2.0, 2.0, 1.0, 1.0, 1.0;

    for (int i = 0; i < 6; ++i) {
        CHECK_NEAR(K_dense(i, i), expected_diag(i), 1e-12);
    }

    // K_local etant l'identite (aucun terme hors-diagonale), K_global doit
    // rester purement diagonale.
    CHECK_TRUE((K_dense - K_dense.diagonal().asDiagonal().toDenseMatrix()).cwiseAbs().maxCoeff() < 1e-12);

    // F_global suit exactement la meme logique de comptage que la diagonale.
    for (int i = 0; i < 6; ++i) {
        CHECK_NEAR(F_global(i), expected_diag(i), 1e-12);
    }
}

TEST_CASE("SystemAssembler : le chemin rapide (CSR bypass) donne le meme resultat que le premier assemblage") {
    Mesh mesh = make_two_element_mesh();
    ConstantIdentityEquation eq;
    PhysicsState state;
    std::vector<FlattenedBC> no_bcs;
    std::vector<long> csr_mapping;
    Eigen::SparseMatrix<double> K_global;
    Eigen::VectorXd F_global;

    SystemAssembler::assemble(mesh.get_num_elements(), mesh.get_num_nodes(), true,
                            K_global, F_global, csr_mapping, eq, no_bcs, state, mesh);

    Eigen::MatrixXd K_first = Eigen::MatrixXd(K_global);
    Eigen::VectorXd F_first = F_global;

    SystemAssembler::assemble(mesh.get_num_elements(), mesh.get_num_nodes(), false,
                            K_global, F_global, csr_mapping, eq, no_bcs, state, mesh);

    Eigen::MatrixXd K_second = Eigen::MatrixXd(K_global);

    CHECK_TRUE((K_first - K_second).cwiseAbs().maxCoeff() < 1e-12);
    CHECK_TRUE((F_first - F_global).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("SystemAssembler : une CL Dirichlet impose value via penalite (diagonale et second membre)") {
    Mesh mesh = make_two_element_mesh();
    ConstantIdentityEquation eq;
    PhysicsState state;
    std::vector<long> csr_mapping;

    // CL Dirichlet sur le noeud 0 : valeur imposee 3.0.
    std::vector<FlattenedBC> bcs = { {0, 3.0, true} };
    Eigen::SparseMatrix<double> K_global;
    Eigen::VectorXd F_global;

    SystemAssembler::assemble(mesh.get_num_elements(), mesh.get_num_nodes(), true,
                            K_global, F_global, csr_mapping, eq, bcs, state, mesh);

    Eigen::MatrixXd K_dense = Eigen::MatrixXd(K_global);

    // max_diag avant CL = 2.0 (noeuds partages) -> penalty = 2.0 * 1e5.
    double expected_penalty = 2.0 * 1e5;

    // Diagonale du noeud 0 = contribution physique (1.0) + penalite.
    CHECK_NEAR(K_dense(0, 0), 1.0 + expected_penalty, 1e-6);

    // F(0) = contribution physique (1.0) + penalty * value.
    CHECK_NEAR(F_global(0), 1.0 + expected_penalty * 3.0, 1e-3);

    // Les autres degres de liberte ne doivent pas etre affectes par cette CL.
    CHECK_NEAR(K_dense(3, 3), 1.0, 1e-12);
    CHECK_NEAR(F_global(3), 1.0, 1e-12);
}