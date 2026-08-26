#pragma once
#include "IO/include/Datafile.h"

// =========================================================================
// TimeManager - pas de temps FIXE, conforme a l'Algorithm 1 du papier
// (Abdollahi & Arias, 2011) : n = 100 increments de charge, Delta t^n = 3e-2,
// sans rejet de pas ni adaptation. La robustesse numerique pres du saut
// instable de propagation est geree par la boucle de Picard interne
// (tolerance / max_iter), pas par une reduction de dt.
// =========================================================================
class TimeManager {
private:
    double current_time;
    const double fixed_dt;
    double end_time;
    int current_step;

public:
    explicit TimeManager(const Datafile& config);

    bool is_finished() const;
    void advance_step();

    double get_time() const { return current_time; }
    double get_dt() const { return fixed_dt; }
    int get_step() const { return current_step; }
};