#pragma once
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

class Chrono {
private:
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> end_time;
    bool is_running;

public:
    Chrono() : is_running(false) {}

    void start() {
        start_time = std::chrono::steady_clock::now();
        is_running = true;
    }

    std::string get_datetime_string() {
        auto now = std::chrono::system_clock::now();
        
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
        
        return ss.str();
    }

    void stop() {
        end_time = std::chrono::steady_clock::now();
        is_running = false;
    }

    double elapsed_ms() const {
        auto current_end = is_running ? std::chrono::steady_clock::now() : end_time;
        std::chrono::duration<double, std::milli> elapsed = current_end - start_time;
        return elapsed.count();
    }
};