#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

class EnsureDir {
public:
    static void prepareOutputDirectory(const std::string& directoryPath) {
        try {
            if (fs::create_directories(directoryPath)) {
                std::cout << "Directory created: " << directoryPath << std::endl;
            } else if (!fs::exists(directoryPath)) {
                std::cerr << "Failed to create directory: " << directoryPath << std::endl;
            }
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
};