#include "Physics/Fracture.h"

Fracture::Fracture(const Datafile& config, const Mesh& mesh) : config(config), mesh(mesh) {
    size_t n_nodes = mesh.get_num_nodes();
    v_prev_iter.setZero(n_nodes);
    v_current.resize(n_nodes); 
    v_current.setZero();

    for (size_t i = 0; i < v_current.size(); ++i) {
        v_current[i] = 1.0;
    }    
}

void Fracture::update_v(double time, const Polarization& polarization, const Mechanics& mechanics, const Electrostatics& electrostatics) {
    for (size_t i = 0; i < v_current.size(); ++i) {
        v_current[i] = 1.0 + 0.1 * sin(time);
    }
}