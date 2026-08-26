#include <filesystem>
#include <string>

#include "Utils/include/Logger.h"

namespace fs = std::filesystem;

class EnsureDir {
public:
    static void prepareOutputDirectory(const std::string& directoryPath) {
        try {
            if (fs::create_directories(directoryPath)) {
                Logger::info("Directory created: ", directoryPath);
            } else if (!fs::exists(directoryPath)) {
                Logger::error("Failed to create directory: ", directoryPath);
            }
        } catch (const fs::filesystem_error& e) {
            Logger::error("Error: ", e.what());
        }
    }
};