// ==============================================================================
// Maillage Gmsh pour PhaseFieldFerro - domaine rectangulaire 200 x 200
// Triangulation non structuree (Delaunay) - elimine le biais directionnel
// d'un maillage structure coupe en diagonale (fixe ou alternee).
// ==============================================================================

Lx = 200.0;
Ly = 200.0;
h  = 1.0;      // taille de maille cible, coherente avec le papier (h ~ 1)

// --- Points des 4 coins ---
Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {Lx,  0.0, 0.0, h};
Point(3) = {Lx,  Ly,  0.0, h};
Point(4) = {0.0, Ly,  0.0, h};

// --- Aretes (bottom, right, top, left) ---
Line(1) = {1, 2};   // bottom : y = 0
Line(2) = {2, 3};   // right  : x = Lx
Line(3) = {3, 4};   // top    : y = Ly
Line(4) = {4, 1};   // left   : x = 0

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Groupes physiques : IDs coherents avec les tags utilises dans
//     MeshGenerator::generate_structured_mesh (1=left, 2=right, 3=bottom, 4=top)
Physical Line("bottom") = {1};
Physical Line("right")  = {2};
Physical Line("top")    = {3};
Physical Line("left")   = {4};
Physical Surface("domain") = {1};

// --- Options de maillage ---
Mesh.Algorithm = 5;         // Delaunay (2D)
Mesh.ElementOrder = 1;      // triangles lineaires (P1), coherent avec ton solveur actuel
Mesh.RecombineAll = 0;      // NE PAS recombiner en quads : on veut des triangles Delaunay purs
Mesh.CharacteristicLengthMin = h * 0.8;
Mesh.CharacteristicLengthMax = h * 1.2;

// Pour generer le maillage en ligne de commande :
//   gmsh -2 rectangle.geo -format msh2 -o rectangle.msh
// (on force le format MSH v2.2 ASCII, plus simple a parser)
