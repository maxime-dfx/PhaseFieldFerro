#pragma once
#include <vector>
#include <Eigen/Dense>
// Forward declarations pour éviter les inclusions circulaires
class Mesh;
class Math;
class Datafile;
class Electrostatics;
class Polarization;
class Mechanics;
class Fracture;

class Fields {
private:
    const Mesh& mesh;
    const Math& math;
    const Datafile& config;
    const Electrostatics& electrostatics;
    const Polarization& polarization;
    const Mechanics& mechanics;
    const Fracture& fracture;

public:
    Fields(const Mesh& m, const Math& ma, const Datafile& c, 
           const Electrostatics& e, const Polarization& p, 
           const Mechanics& me, const Fracture& f);

    // Ajout des types de retour manquants (const std::vector<double>&)
    const Eigen::VectorXd& get_Px() const;
    const Eigen::VectorXd& get_Py() const;
    const Eigen::VectorXd& get_Ex() const;
    const Eigen::VectorXd& get_Ey() const;
    const Eigen::VectorXd& get_v() const;
    const Eigen::VectorXd& get_ux() const;
    const Eigen::VectorXd& get_uy() const;
    const Eigen::VectorXd& get_strain() const;
    const Mesh& get_mesh() const { return mesh; }
    const Mechanics& get_mechanics() const { return mechanics; }
    const Electrostatics& get_electrostatics() const { return electrostatics; }
    const Math& get_math() const { return math; }

};