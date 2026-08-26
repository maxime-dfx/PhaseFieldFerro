#pragma once
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

// Identifiant unique de run, utilise pour construire un dossier de sortie
// dedie a chaque execution et eviter que deux runs partageant le meme
// config.toml (meme output_dir) n'ecrivent dans les memes fichiers.
//
// Priorite a SLURM_JOB_ID quand il est present : ca fait correspondre
// directement le dossier de sortie au numero de job visible dans les logs
// SLURM ("Job SLURM 610 sur w102..."). Fallback sur un horodatage sinon
// (execution locale hors SLURM).
namespace RunId {

inline std::string make() {
    if (const char* slurm_id = std::getenv("SLURM_JOB_ID")) {
        return std::string("job") + slurm_id;
    }

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}

}
