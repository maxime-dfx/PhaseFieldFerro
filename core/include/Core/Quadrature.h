#pragma once
#include <vector>
#include <array>
#include <cmath>

struct GaussPoint2D {
    double xi, eta, weight;
};

class Quadrature {
public:
    static std::vector<GaussPoint2D> get_gauss_2x2() {
        static const std::vector<GaussPoint2D> points = {
            {-1.0/std::sqrt(3.0), -1.0/std::sqrt(3.0), 1.0},
            { 1.0/std::sqrt(3.0), -1.0/std::sqrt(3.0), 1.0},
            { 1.0/std::sqrt(3.0),  1.0/std::sqrt(3.0), 1.0},
            {-1.0/std::sqrt(3.0),  1.0/std::sqrt(3.0), 1.0}
        };
        return points;
    }
    std::vector<GaussPoint2D> get_quad_gauss_points() {
        double pt = 1.0 / std::sqrt(3.0);
        return {
            {-pt, -pt, 1.0}, {pt, -pt, 1.0}, {pt, pt, 1.0}, {-pt, pt, 1.0}
        };
    }
};