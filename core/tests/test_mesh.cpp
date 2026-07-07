#include "Core/Mesh.h"
#include "Core/Quadrature.h"
#include "Core/ShapeFunctions.h"
#include "Utils/Logger.h"
#include "IO/Datafile.h"

int main() {
    // // 1. Créer un maillage 10x10 (domaine [0,1]x[0,1])
    // Mesh mesh(Datafile("config.txt"));
    // Logger::info("Nombre de nœuds : " + std::to_string(mesh.get_num_nodes()));
    // Logger::info("Nombre de quadrangles : " + std::to_string(mesh.get_num_elements()));
    // Logger::info("Nombre de triangles : " + std::to_string(mesh.get_num_elements()));

    // // 2. Vérifier les tags de bord
    // for (int edge = 1; edge <= 4; ++edge) {
    //     auto nodes = mesh.get_boundary_nodes(edge);
    //     std::string name = (edge==1)?"bas":(edge==2)?"droite":(edge==3)?"haut":"gauche";
    //     Logger::info("Bord " + name + " : " + std::to_string(nodes.size()) + " nœuds");
    // }

    // // 3. Vérifier l'aire de chaque élément (via intégration numérique)
    // double total_area = 0.0;
    // auto gauss_points = Quadrature::get_gauss_2x2();
    // for (int e = 0; e < mesh.get_num_elements(); ++e) {
    //     auto coords = mesh.get_element_coords(e);
    //     double area = 0.0;
    //     for (const auto& gp : gauss_points) {
    //         auto dN_xi_eta = ShapeFunctions::get_shape_function_gradients(gp.xi, gp.eta);
    //         auto [dN_xy, detJ] = ShapeFunctions::compute_physical_derivatives(coords, dN_xi_eta);
    //         area += gp.weight * detJ;
    //     }
    //     total_area += area;
    //     double expected_area = mesh.get_dx() * mesh.get_dy();
    //     if (std::abs(area - expected_area) > 1e-10) {
    //         Logger::warning("Élément " + std::to_string(e) + " : aire=" + std::to_string(area) + " (attendu " + std::to_string(expected_area) + ")");
    //     }
    // }
    // Logger::info("Aire totale = " + std::to_string(total_area) + " (attendu 1.0)");

    // // 4. Vérifier la somme des fonctions de forme
    // double xi = 0.3, eta = -0.4;
    // auto N = ShapeFunctions::get_shape_functions(xi, eta);
    // double sum_N = 0;
    // for (double v : N) sum_N += v;
    // if (std::abs(sum_N - 1.0) > 1e-12)
    //     Logger::error("Somme des fonctions de forme = " + std::to_string(sum_N) + " (≠1)");
    // else
    //     Logger::success("Somme des fonctions de forme ok.");

    // Logger::success("Tests de base réussis !");
    // return 0;
}