#pragma once
#include "IO/Datafile.h"

class TimeManager {
private:
    double current_time;
    double current_dt;
    double end_time;
    int current_step;
    
    double dt_min;
    double dt_max;

public:
    explicit TimeManager(const Datafile& config);

    bool is_finished() const;
    void advance_step();
    void adapt_dt_after_failure();
    void adapt_dt_after_success();

    double get_time() const { return current_time; }
    double get_dt() const { return current_dt; }
    int get_step() const { return current_step; }
};