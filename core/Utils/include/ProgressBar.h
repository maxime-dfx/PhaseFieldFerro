#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

class ProgressBar {
public:
    ProgressBar(int max_steps, const std::string& prefix)
        : m_max_steps(max_steps), m_prefix(prefix) 
    {
        m_start_time = std::chrono::high_resolution_clock::now();
    }

    void update(int current_step, double current_delta) {
        const char spinner[4] = {'|', '/', '-', '\\'};

        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - m_start_time;
        
        int eta_seconds = 0;
        if (current_step > 0) {
            eta_seconds = static_cast<int>((elapsed.count() / current_step) * (m_max_steps - current_step));
        }

        int h = eta_seconds / 3600;
        int m = (eta_seconds % 3600) / 60;
        int s = eta_seconds % 60;
        
        int percentage = (m_max_steps > 0) ? (int)((current_step * 100.0) / m_max_steps) : 0;
        
        std::cout << "\r\033[K" << m_prefix << " " << spinner[current_step % 4] 
                << " Etape : " << current_step << " / " << m_max_steps 
                << " [" << percentage << "%] "
                << "- ETA : " << std::setfill('0') << std::setw(2) << h << ":"
                            << std::setfill('0') << std::setw(2) << m << ":"
                            << std::setfill('0') << std::setw(2) << s << " ";
        
        if (current_delta > 0.0) {
            std::cout << "- Delta : " << std::scientific << std::setprecision(4) << current_delta;
        }
        
        std::cout << std::flush;
    }

    void finish() {
        std::cout << "\n";
    }

private:
    int m_max_steps;
    std::string m_prefix;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start_time;
};