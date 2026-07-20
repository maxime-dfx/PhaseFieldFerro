// include/Physics/IPhysicsModule.h
#pragma once

class IPhysicsModule {
public:
    virtual ~IPhysicsModule() = default;
    
    virtual void initialize() = 0;
    virtual int compute_step(double time, double dt) = 0;
    
    // Délégation de la gestion des états (remplace les méthodes globales de Simulation)
    virtual void save_previous_state() = 0;
    virtual void restore_previous_state() = 0;
    virtual void update_history() = 0;
};