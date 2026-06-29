#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

class ProgressBar {
public:
    ProgressBar(int max_steps, const std::string& prefix = "[CALCUL]");

    void update(int current_step, double current_delta = 0.0);

    void finish();

private:
    int m_max_steps;
    std::string m_prefix;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start_time;
};