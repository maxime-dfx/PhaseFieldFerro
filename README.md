# PhaseFieldFerro Solver

## Description
PhaseFieldFerro est un solveur multi-physique haute performance dédié à la simulation mésoscopique des matériaux ferroélectriques. 
Ce projet implémente un couplage fort entre :
- **Thermodynamique :** Évolution de la polarisation via l'équation de Ginzburg-Landau (schéma IMEX).
- **Électrostatique :** Résolution de l'équation de Poisson pour le champ électrique.
- **Mécanique :** Élasticité linéaire et électrostriction via une méthode d'éléments finis (FEM) sur grille structurée.

Le code est optimisé avec `Eigen` pour l'algèbre linéaire creuse et permet l'exportation multi-physique au format VTK pour une visualisation avancée sous ParaView.

## Fonctionnalités
- Solveur multi-physique couplé.
- Gestion flexible des conditions aux limites (Dirichlet homogène/non-homogène).
- Exportation unifiée des résultats (Polarisation, Champ Électrique, Déplacement, Contraintes).
- Architecture modulaire et performante (LDLT factorisation, stockage creux).

## Prérequis
- Compilateur C++ (supportant C++17).
- [CMake](https://cmake.org/) (version 3.10+).
- [Eigen3](https://eigen.tuxfamily.org/) (bibliothèque d'algèbre linéaire).

## Installation
```bash
# Cloner le dépôt
git clone [https://github.com/ton-pseudo/PhaseFieldFerro-Solver.git](https://github.com/ton-pseudo/PhaseFieldFerro-Solver.git)
cd PhaseFieldFerro-Solver

# Créer le répertoire de build
mkdir build && cd build

# Configurer et compiler
cmake ..
make -j4
