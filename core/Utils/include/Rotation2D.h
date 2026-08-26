#pragma once

#include <array>
#include <cmath>
#include <utility>

#include <Eigen/Core>

struct Rotation2D final {
    double theta = 0.0;
    double cos_t = 1.0;
    double sin_t = 0.0;

    explicit Rotation2D(double angle = 0.0)
        : theta(angle), cos_t(std::cos(angle)), sin_t(std::sin(angle)) {}

    Eigen::Vector2d to_local(const Eigen::Vector2d& v_global) const {
        return {cos_t * v_global.x() + sin_t * v_global.y(),
                -sin_t * v_global.x() + cos_t * v_global.y()};
    }

    Eigen::Vector2d to_global(const Eigen::Vector2d& v_local) const {
        return {cos_t * v_local.x() - sin_t * v_local.y(),
                sin_t * v_local.x() + cos_t * v_local.y()};
    }

    Eigen::Matrix2d to_local(const Eigen::Matrix2d& T_global) const {
        const Eigen::Matrix2d R = rotation_matrix();
        return R.transpose() * T_global * R;
    }

    Eigen::Matrix2d to_global(const Eigen::Matrix2d& T_local) const {
        const Eigen::Matrix2d R = rotation_matrix();
        return R * T_local * R.transpose();
    }

    Eigen::Vector3d to_local(const Eigen::Vector3d& voigt_global) const {
        return voigt_to_global_or_local(voigt_global, true);
    }

    Eigen::Vector3d to_global(const Eigen::Vector3d& voigt_local) const {
        return voigt_to_global_or_local(voigt_local, false);
    }

    Eigen::Matrix3d to_local(const Eigen::Matrix3d& C_global) const {
        return rotate_stiffness(C_global, true);
    }

    Eigen::Matrix3d to_global(const Eigen::Matrix3d& C_local) const {
        return rotate_stiffness(C_local, false);
    }

private:
    Eigen::Matrix2d rotation_matrix() const {
        Eigen::Matrix2d R;
        R << cos_t, -sin_t,
             sin_t,  cos_t;
        return R;
    }

    Eigen::Vector3d voigt_to_global_or_local(const Eigen::Vector3d& voigt, bool use_local_frame) const {
        Eigen::Matrix2d T;
        T << voigt(0), voigt(2),
             voigt(2), voigt(1);
        const Eigen::Matrix2d rotated = use_local_frame ? this->to_local(T) : this->to_global(T);
        return {rotated(0, 0), rotated(1, 1), rotated(0, 1)};
    }

    static int voigt_index(int i, int j) {
        return (i == j) ? i : 2;
    }

    Eigen::Matrix3d rotate_stiffness(const Eigen::Matrix3d& C_in, bool use_local_frame) const {
        std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> C4{};
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    for (int l = 0; l < 2; ++l) {
                        C4[i][j][k][l] = C_in(voigt_index(i, j), voigt_index(k, l));
                    }
                }
            }
        }

        const Eigen::Matrix2d R = rotation_matrix();
        const Eigen::Matrix2d Rt = R.transpose();
        const Eigen::Matrix2d& A = use_local_frame ? Rt : R;

        std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> C4_rot{};
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    for (int l = 0; l < 2; ++l) {
                        double sum = 0.0;
                        for (int a = 0; a < 2; ++a) {
                            for (int b = 0; b < 2; ++b) {
                                for (int c = 0; c < 2; ++c) {
                                    for (int d = 0; d < 2; ++d) {
                                        sum += A(i, a) * A(j, b) * A(k, c) * A(l, d) * C4[a][b][c][d];
                                    }
                                }
                            }
                        }
                        C4_rot[i][j][k][l] = sum;
                    }
                }
            }
        }

        Eigen::Matrix3d C_out;
        static const std::array<std::pair<int, int>, 3> pairs = {{{0, 0}, {1, 1}, {0, 1}}};
        for (int I = 0; I < 3; ++I) {
            for (int J = 0; J < 3; ++J) {
                auto [i, j] = pairs[I];
                auto [k, l] = pairs[J];
                C_out(I, J) = C4_rot[i][j][k][l];
            }
        }
        return C_out;
    }
};