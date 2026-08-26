#include "Tests/include/TestFramework.h"
#include "Mesh/include/MeshGenerators.h"
#include <fstream>
#include <cstdio>
#include <stdexcept>

// =========================================================================
// Tests GmshReader : parsing de fichiers .msh (format Gmsh ASCII v2). C'est
// le chemin utilise en production quand mesh.calcul_mesh=false (maillage
// externe charge depuis input/rectangle.msh), donc particulierement
// sensible - contrairement a StructuredMeshGenerator, il n'etait couvert
// par aucun test avant celui-ci.
//
// Maillage de test (rectangle 2x1, 6 noeuds, IDs Gmsh volontairement NON
// contigus et non ordonnes pour verifier le mapping gmsh_id -> index local) :
//
//   30--20--60
//   |E0 |E1 |    physical tags bord : bottom=1, right=2, top=3, left=4
//   10--50--40    physical tag domaine (triangles) = 5
//
// Triangulation : E0=(10,50,30)+(50,20,30), E1=(50,40,20)+(40,60,20)
// simplifie ici a 2 triangles seulement pour rester lisible (cf fichier).
// =========================================================================

namespace {

// Ecrit un .msh minimal sur disque et retourne son chemin. Reutilise le
// pattern deja en place dans test_datafile.cpp (fichier temporaire,
// supprime juste apres construction de l'objet a tester).
std::string write_msh_basic(const std::string& path) {
    std::ofstream out(path);
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    out << "$PhysicalNames\n5\n";
    out << "1 1 \"bottom\"\n1 2 \"right\"\n1 3 \"top\"\n1 4 \"left\"\n2 5 \"domain\"\n";
    out << "$EndPhysicalNames\n";
    // 6 noeuds, IDs Gmsh non contigus et non ordonnes : 10,50,30,20,40,60.
    out << "$Nodes\n6\n";
    out << "10 0 0 0\n";
    out << "50 1 0 0\n";
    out << "30 0 1 0\n";
    out << "20 1 1 0\n";
    out << "40 2 0 0\n";
    out << "60 2 1 0\n";
    out << "$EndNodes\n";
    // Elements : segments de bord (type 1, 2 noeuds) + triangles (type 2, 3 noeuds).
    // Format : id elm_type n_tags [tags...] node_ids...
    // Coin (noeud 30) partage par 2 segments de bord differents (left=4 puis top=3) :
    // sert a documenter/verifier le comportement "dernier tag ecrit gagne".
    out << "$Elements\n6\n";
    out << "1 1 2 4 4 10 30\n";   // segment gauche (physical=4), noeuds 10,30
    out << "2 1 2 3 3 30 20\n";   // segment haut   (physical=3), noeuds 30,20 (30 re-tag ici)
    out << "3 1 2 1 1 10 50\n";   // segment bas    (physical=1), noeuds 10,50
    out << "4 2 2 5 5 10 50 30\n"; // triangle 1 (physical=5)
    out << "5 2 2 5 5 50 20 30\n"; // triangle 2 (physical=5)
    out << "6 2 2 5 5 50 40 20\n"; // triangle 3 (physical=5)
    out << "$EndElements\n";
    out.close();
    return path;
}

} // namespace

TEST_CASE("GmshReader : parse correctement les coordonnees des noeuds (IDs non contigus)") {
    std::string path = write_msh_basic("test_gmsh_basic.msh");
    Mesh mesh = MeshGenerators::load_from_gmsh(path, 2.0, 1.0);
    std::remove(path.c_str());

    CHECK_TRUE(mesh.get_num_nodes() == 6);
    // L'ordre des noeuds locaux suit l'ordre d'apparition dans $Nodes, pas
    // les IDs Gmsh eux-memes : noeud local 0 = Gmsh id 10 = (0,0), etc.
    const auto& nodes = mesh.get_nodes();
    CHECK_NEAR(nodes[0].x, 0.0, 1e-12); CHECK_NEAR(nodes[0].y, 0.0, 1e-12);
    CHECK_NEAR(nodes[1].x, 1.0, 1e-12); CHECK_NEAR(nodes[1].y, 0.0, 1e-12);
    CHECK_NEAR(nodes[4].x, 2.0, 1e-12); CHECK_NEAR(nodes[4].y, 0.0, 1e-12);
    CHECK_NEAR(nodes[5].x, 2.0, 1e-12); CHECK_NEAR(nodes[5].y, 1.0, 1e-12);
}

TEST_CASE("GmshReader : triangles construits avec le bon mapping gmsh_id -> index local") {
    std::string path = write_msh_basic("test_gmsh_triangles.msh");
    Mesh mesh = MeshGenerators::load_from_gmsh(path, 2.0, 1.0);
    std::remove(path.c_str());

    // 3 triangles dans le fichier -> 3 elements (les 3 segments de bord ne
    // doivent PAS etre comptes comme des elements du maillage).
    CHECK_TRUE(mesh.get_num_elements() == 3);
    const auto& elems = mesh.get_elements();
    for (const auto& e : elems) {
        CHECK_TRUE(e.get_num_nodes() == 3);
        CHECK_TRUE(e.ref_tag == 5); // physical tag "domain"
    }
    // Premier triangle du fichier : noeuds Gmsh (10,50,30) -> index locaux (0,1,2).
    CHECK_TRUE(mesh.get_node_index(0, 0) == 0);
    CHECK_TRUE(mesh.get_node_index(0, 1) == 1);
    CHECK_TRUE(mesh.get_node_index(0, 2) == 2);
}

TEST_CASE("GmshReader : les noeuds de bord recoivent le physical tag Gmsh du segment") {
    std::string path = write_msh_basic("test_gmsh_bc_tags.msh");
    Mesh mesh = MeshGenerators::load_from_gmsh(path, 2.0, 1.0);
    std::remove(path.c_str());

    const auto& nodes = mesh.get_nodes();
    // Noeud local 0 (Gmsh id 10) : uniquement sur le segment bas (physical=1)
    // et le segment gauche (physical=4). Le segment bas est ecrit APRES le
    // segment gauche dans le fichier -> dernier tag ecrit gagne (pas
    // d'accumulation / de priorite documentee dans GmshReader.cpp, le code
    // ecrase simplement ref_tag a chaque passage). On documente ici le
    // comportement REEL du code, pour detecter toute regression silencieuse
    // si cet ordre venait a changer.
    CHECK_TRUE(nodes[0].ref_tag == 1);
    // Noeud local 4 (Gmsh id 40) : n'apparait dans AUCUN segment de bord du
    // fichier de test -> doit garder sa valeur par defaut (0).
    CHECK_TRUE(nodes[4].ref_tag == 0);
}

TEST_CASE("GmshReader : les dimensions Lx/Ly passees au constructeur sont bien reportees sur le Mesh") {
    std::string path = write_msh_basic("test_gmsh_dims.msh");
    Mesh mesh = MeshGenerators::load_from_gmsh(path, 7.5, 3.25);
    std::remove(path.c_str());

    CHECK_NEAR(mesh.get_Lx(), 7.5, 1e-12);
    CHECK_NEAR(mesh.get_Ly(), 3.25, 1e-12);
}

TEST_CASE("GmshReader : fichier introuvable leve une exception explicite") {
    bool threw = false;
    try {
        Mesh mesh = MeshGenerators::load_from_gmsh("chemin_qui_nexiste_vraiment_pas.msh", 1.0, 1.0);
        (void)mesh;
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK_TRUE(threw);
}

TEST_CASE("GmshReader : type d'element Gmsh inconnu est ignore silencieusement") {
    std::string path = "test_gmsh_bad_type.msh";
    std::ofstream out(path);
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    out << "$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 0\n$EndNodes\n";
    
    // Type 99 = type inconnu. Le reader doit l'ignorer et ne pas lever d'exception.
    // (L'ancien test utilisait le type 3, qui est desormais supporte comme QUAD4).
    out << "$Elements\n1\n1 99 2 5 5 1 2 3 4\n$EndElements\n";
    out.close();
    
    bool parsed_without_crashing = true;
    int num_elements = -1;
    
    try {
        Mesh mesh = MeshGenerators::load_from_gmsh(path, 1.0, 1.0);
        num_elements = mesh.get_num_elements();
    } catch (const std::exception& e) {
        parsed_without_crashing = false;
    }
    
    std::remove(path.c_str());
    
    // On verifie que la fonction est allee au bout sans lever d'exception
    CHECK_TRUE(parsed_without_crashing);
    // On verifie que l'element inconnu a bien ete ignore (le maillage est donc vide)
    CHECK_TRUE(num_elements == 0);
}