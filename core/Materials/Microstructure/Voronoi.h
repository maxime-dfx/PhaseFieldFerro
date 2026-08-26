#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "Materials/Microstructure/Crystal.h"
#include "Utils/include/Profiling.h"

// VoronoiTessellation
// -----------------------------------------------------------------------
// Genere une microstructure polycristalline synthetique par tessellation
// de Voronoi : num_grains germes tires aleatoirement sur le domaine
// rectangulaire [0,Lx]x[0,Ly], chacun avec une orientation cristalline
// aleatoire independante. Le critere d'appartenance d'un point (x,y) a
// un grain est la distance euclidienne au germe le plus proche (Voronoi
// standard).
//
// Reproductibilite : la seed est fournie explicitement par l'appelant
// (cf. config.crystal.seed) et pilote un std::mt19937 deterministe -
// deux runs avec la meme seed produisent exactement la meme
// microstructure.
class VoronoiTessellation {
public:
    // Genere num_grains germes + orientations sur [0,Lx]x[0,Ly].
    // num_grains doit etre >= 1. seed pilote le tirage aleatoire (germes
    // ET orientations) de maniere reproductible.
    inline VoronoiTessellation(double Lx, double Ly, int num_grains, std::uint32_t seed) {
        PROFILE_ZONE_NC("VoronoiTessellation::VoronoiTessellation", PROFILE_COLOR_SEQUENTIAL);
        if (num_grains < 1) {
            throw std::invalid_argument("VoronoiTessellation: num_grains doit etre >= 1");
        }

        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> dist_x(0.0, Lx);
        std::uniform_real_distribution<double> dist_y(0.0, Ly);
        std::uniform_real_distribution<double> dist_theta(0.0, kTwoPi);

        grains_.reserve(static_cast<std::size_t>(num_grains));
        for (int g = 0; g < num_grains; ++g) {
            Crystal grain;
            grain.id = g;
            grain.seed_x = dist_x(rng);
            grain.seed_y = dist_y(rng);
            grain.theta = dist_theta(rng);
            grains_.push_back(grain);
        }
    }

    const std::vector<Crystal>& grains() const { return grains_; }

    // Renvoie l'id (0-based, index dans grains()) du grain dont le germe
    // est le plus proche du point (x,y).
    inline int nearest_grain(double x, double y) const {
        int best_id = -1;
        double best_dist_sq = std::numeric_limits<double>::max();

        for (const Crystal& grain : grains_) {
            double dx = x - grain.seed_x;
            double dy = y - grain.seed_y;
            double dist_sq = dx * dx + dy * dy;
            if (dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best_id = grain.id;
            }
        }
        return best_id;
    }

private:
    std::vector<Crystal> grains_;
    static constexpr double kTwoPi = 6.28318530717958647692;
};
