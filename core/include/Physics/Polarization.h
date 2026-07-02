#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Core/BoundaryManager.h"

class Math;
class Fracture;
class Mechanics;
class Electrostatics;

class Polarization {
    private:
        Eigen::VectorXd Px_current;
        Eigen::VectorXd Py_current;
        Eigen::VectorXd Px_prev_iter;
        Eigen::VectorXd Py_prev_iter;
        Eigen::VectorXd Px_n;
        Eigen::VectorXd Py_n;

        double Px_0;
        double Py_0;
        const Datafile& config;
        const Mesh& mesh;
        const BoundaryManager& bc_manager;

        // Solveur persistant : analyzePattern() (reordonnancement AMD, le plus
        // couteux) une seule fois pour toute la simulation -- le PATTERN de
        // sparsite est identique pour Px et Py (memes elements/connectivite),
        // meme si les CL Dirichlet different entre les deux composantes.
        // factorize() reste appele deux fois par pas (une par composante),
        // car les valeurs de la matrice different une fois les CL appliquees.
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
        bool pattern_analyzed_ = false;

        // Assemble la matrice de base A (commune a Px/Py, AVANT application
        // des CL Dirichlet) et les deux seconds membres b_x, b_y en une seule
        // passe sur les elements.
        void assemble_system(const Fracture& fracture, const Mechanics& mechanics,
                              const Electrostatics& electrostatics, const Math& math,
                              Eigen::SparseMatrix<double>& A_sparse,
                              Eigen::VectorXd& b_x, Eigen::VectorXd& b_y);

    public:
        // Initialization of fracture based on the configuration and mesh
        Polarization(const Datafile& config, const Mesh& mesh, const BoundaryManager& bc_manager);
        void apply_boundary_conditions(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const std::vector<NodeBC>& bcs);

        // Nouvelle methode : assemble la matrice de base une seule fois pour
        // Px et Py, partage un unique analyzePattern() (reordonnancement AMD),
        // et fait deux factorize()+solve() (un par composante, car les CL
        // Dirichlet Px/Py peuvent differer). A appeler UNE FOIS par pas/
        // iteration de Picard, a la place de deux appels separes
        // update_Px()+update_Py().
        void update_polarization(double time, const Fracture& fracture, const Mechanics& mechanics,
                                  const Electrostatics& electrostatics, const Math& math);

        // --- Conservees pour compatibilite ascendante ---
        // ATTENTION: si les deux sont appelees a la suite dans la boucle
        // appelante, le gain de factorisation partagee est perdu (2 resolutions
        // completes). Prefer update_polarization() directement.
        [[deprecated("Utiliser update_polarization() pour beneficier de la factorisation partagee Px/Py")]]
        void update_Px(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
            update_polarization(time, fracture, mechanics, electrostatics, math);
        }

        [[deprecated("Utiliser update_polarization() pour beneficier de la factorisation partagee Px/Py")]]
        void update_Py(double time, const Fracture& fracture, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math) {
            update_polarization(time, fracture, mechanics, electrostatics, math);
        }

        const Eigen::VectorXd& get_Px() const { return Px_current; }
        const Eigen::VectorXd& get_Py() const { return Py_current; }

        double calculate_error() {
            return (Px_current - Px_prev_iter).norm() + (Py_current - Py_prev_iter).norm();
        }

        void save_previous_iteration() {
            Px_prev_iter = Px_current;
            Py_prev_iter = Py_current;
        }

        void freeze_time_step() {
            Px_n = Px_current;
            Py_n = Py_current;
        }
        void set_state(const Eigen::VectorXd& Px, const Eigen::VectorXd& Py,
            const Eigen::VectorXd& Px_n_in, const Eigen::VectorXd& Py_n_in) {
            Px_current = Px;
            Py_current = Py;
            Px_n = Px_n_in;
            Py_n = Py_n_in;
            Px_prev_iter = Px;
            Py_prev_iter = Py;
        }

};