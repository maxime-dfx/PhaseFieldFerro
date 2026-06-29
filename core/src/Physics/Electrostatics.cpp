#include "Physics/Electrostatics.h"
#include "Core/Mesh.h"
#include "IO/Datafile.h"
#include "Utils/Logger.h"
#include "Core/Quadrature.h"
#include <cmath>
#include <cstdlib>

Electrostatics::Electrostatics(const Datafile& config, const Mesh& mesh) 
    : config(config), mesh(mesh) 
{
    size_t n_nodes = mesh.get_num_nodes();
    Ex_prev_iter.setZero(n_nodes);
    Ey_prev_iter.setZero(n_nodes);
    Ex_current.resize(n_nodes); 
    Ex_current.setZero();
    Ey_current.resize(n_nodes);
    Ey_current.setZero();
    
    if (config.get_initial_polarization() == PolarizationInitializationType::UNIFORM) {
        Ex_0 = config.get_Ex_0();
        Ey_0 = config.get_Ey_0();
        for (size_t i = 0; i < n_nodes; ++i) {
            Ex_current[i] = Ex_0;
            Ey_current[i] = Ey_0;
        }
    } else if (config.get_initial_polarization() == PolarizationInitializationType::RANDOM) {
        for (size_t i = 0; i < n_nodes; ++i) {
            Ex_current[i] = static_cast<double>(rand()) / RAND_MAX; 
            Ey_current[i] = static_cast<double>(rand()) / RAND_MAX;
        }
    } else {
        Logger::error("Type d'initialisation de polarisation non reconnu.");
    }
}

void Electrostatics::update_Ex(double time, const Polarization& polarization, const Mechanics& mechanics, const Fracture& fracture) {
        for (size_t i = 0; i < Ex_current.size(); ++i) {
        Ex_current[i] = Ex_0 * std::sin(time);
    }
}

void Electrostatics::update_Ey(double time, const Polarization& polarization, const Mechanics& mechanics, const Fracture& fracture) {
        for (size_t i = 0; i < Ey_current.size(); ++i) {
        Ey_current[i] = Ey_0 * std::sin(time);
    }
}

// 2. Le nom doit correspondre EXACTEMENT à celui du .h
double Electrostatics::get_Ex_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    (void)elem;
    (void)gp;
    return 0.0; 
}

double Electrostatics::get_Ey_at_gp(const Element& elem, const GaussPoint2D& gp) const {
    (void)elem;
    (void)gp;
    return 0.0; 
}