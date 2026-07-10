#ifndef PHYSICS_FRACTURE_H
#define PHYSICS_FRACTURE_H

#include <Eigen/Core>
#include "IO/Datafile.h"
#include "Core/Mesh.h"

class Polarization;
class Mechanics;
class Electrostatics;
class Math;

class Fracture {
private:
    const Datafile& config;
    const Mesh& mesh;

    Eigen::VectorXd v_current;
    Eigen::VectorXd v_prev_iter;
    Eigen::VectorXd v_n;
    Eigen::VectorXd v_backup;

public:
    Fracture(const Datafile& config, const Mesh& mesh);

    void update_v(double dt, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics, const Math& math);
    void enforce_physical_bounds(const Eigen::VectorXd& v_new);

    double calculate_error() const;
    void save_previous_iteration();
    void save_previous_state();
    void restore_previous_state();
    void update_history();

    const Eigen::VectorXd& get_v() const { return v_current; }
};

#endif // PHYSICS_FRACTURE_H