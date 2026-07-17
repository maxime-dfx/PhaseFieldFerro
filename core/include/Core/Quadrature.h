#pragma once
#include <vector>
#include <array>
#include <cmath>

struct GaussPoint2D {
    double xi, eta, weight;
};

class Quadrature {
public:
    static std::vector<GaussPoint2D> get_gauss_3_points_tri() {
        return {
            {1.0/6.0, 1.0/6.0, 1.0/6.0},
            {2.0/3.0, 1.0/6.0, 1.0/6.0},
            {1.0/6.0, 2.0/3.0, 1.0/6.0}
        };
    }

    static std::vector<GaussPoint2D> get_gauss_3x3() {
        double pt = std::sqrt(3.0 / 5.0);
        double w_center = 8.0 / 9.0;
        double w_edge = 5.0 / 9.0;
        return {
            {-pt, -pt, w_edge*w_edge}, {0.0, -pt, w_center*w_edge}, {pt, -pt, w_edge*w_edge},
            {-pt,  0.0, w_edge*w_center}, {0.0,  0.0, w_center*w_center}, {pt,  0.0, w_edge*w_center},
            {-pt,  pt, w_edge*w_edge}, {0.0,  pt, w_center*w_edge}, {pt,  pt, w_edge*w_edge}
        };
    }
};