#include "../include/TimeManager.h"
#include "Utils/include/Logger.h"

TimeManager::TimeManager(const Datafile& config)
    : current_time(0.0),
      fixed_dt(config.simulation.dt),
      end_time(config.simulation.total_time),
      current_step(0)
{
}

bool TimeManager::is_finished() const {
    return current_time >= end_time;
}

void TimeManager::advance_step() {
    current_time += fixed_dt;
    current_step++;
}
