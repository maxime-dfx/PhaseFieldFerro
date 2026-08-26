#pragma once

#include <Eigen/Core>

namespace MaterialMath {

inline double symmetrized_shear(const Eigen::Matrix2d& tensor) {
    return 0.5 * (tensor(0, 1) + tensor(1, 0));
}

inline Eigen::Matrix2d make_symmetric_tensor(double xx, double yy, double xy) {
    Eigen::Matrix2d tensor = Eigen::Matrix2d::Zero();
    tensor(0, 0) = xx;
    tensor(1, 1) = yy;
    tensor(0, 1) = xy;
    tensor(1, 0) = xy;
    return tensor;
}

} // namespace MaterialMath
