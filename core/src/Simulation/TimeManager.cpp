#include "Simulation/TimeManager.h"
#include "Utils/Logger.h"
#include <stdexcept>

TimeManager::TimeManager(const Datafile& config) 
    : current_time(0.0), 
      current_dt(config.simulation.dt), 
      end_time(config.simulation.total_time), 
      current_step(0),
      dt_min(1e-8), 
      dt_max(config.simulation.dt * 5.0) 
{}

bool TimeManager::is_finished() const {
    return current_time >= end_time;
}

void TimeManager::advance_step() {
    current_time += current_dt;
    current_step++;
}

void TimeManager::adapt_dt_after_failure() {
    Logger::warning("Non-convergence a t=" + std::to_string(current_time + current_dt) + ". Reduction de dt...");
    current_dt *= 0.5;
    
    if (current_dt < dt_min) {
        throw std::runtime_error("Erreur fatale : dt est devenu trop petit (< dt_min) ! Rupture numerique.");
    }
}

void TimeManager::adapt_dt_after_success() {
    // Facultatif : on pourrait remonter le pas de temps s'il devient très stable, 
    // par exemple si le nombre d'itérations de la boucle de Picard est très bas.
    // current_dt = std::min(current_dt * 1.1, dt_max);
}