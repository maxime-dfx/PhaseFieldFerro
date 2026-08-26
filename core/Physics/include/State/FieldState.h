#pragma once
#include <Eigen/Core>
struct FieldState {
    Eigen::VectorXd current, prev_iter, n, backup;
    void resize_and_reset(int size) { current.setZero(size); prev_iter.setZero(size); n.setZero(size); backup.setZero(size); }
    void set(const Eigen::VectorXd& val) { current = val; n = val; }
    void save_iteration() { prev_iter = current; }
    void save_state() { backup = current; }
    void restore_state() { current = backup; }
    void update_history() { n = current; }
};
