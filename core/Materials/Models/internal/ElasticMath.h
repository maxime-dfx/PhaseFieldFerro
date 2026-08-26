#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "IO/include/ConfigTypes.h"
#include "Materials/Models/internal/Math/Math.h"

struct ElasticTerms {
    double W_mec_energy = 0.0;
    double W_elas_energy = 0.0;
    Eigen::Matrix3d elastic_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d sigma_0 = Eigen::Vector3d::Zero();
};

class ElasticMath {
public:
    static double compute_mec_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain, const MaterialConfig& mat) {
        const double p1 = P(0);
        const double p2 = P(1);
        const double p1_2 = p1 * p1;
        const double p2_2 = p2 * p2;
        return -0.5 * (mat.b1 * strain(0,0) * p1_2 +
                       mat.b2 * strain(1,1) * p1_2 +
                       mat.b1 * strain(1,1) * p2_2 +
                       mat.b2 * strain(0,0) * p2_2 +
                       4.0 * mat.b3 * strain(0,1) * p1 * p2);
    }

    static double compute_elas_energy(const Eigen::Matrix2d& strain, const MaterialConfig& mat) {
        return 0.5 * mat.c1 * (strain(0,0) * strain(0,0) + strain(1,1) * strain(1,1))
             + mat.c2 * strain(0,0) * strain(1,1)
             + 0.5 * mat.c3 * (strain(0,1) * strain(0,1) + strain(1,0) * strain(1,0));
    }

    static double compute_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain, const MaterialConfig& mat) {
        return compute_mec_energy(P, strain, mat) + compute_elas_energy(strain, mat);
    }

    static Eigen::Matrix3d compute_elastic_matrix(const MaterialConfig& mat) {
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        C(0,0) = mat.c1;
        C(1,1) = mat.c1;
        C(0,1) = mat.c2;
        C(1,0) = mat.c2;
        C(2,2) = mat.c3 / 2.0;
        return C;
    }

    static Eigen::Vector3d compute_sigma_0(const Eigen::Vector2d& P, const MaterialConfig& mat) {
        Eigen::Vector3d sigma_0 = Eigen::Vector3d::Zero();
        const double p1_2 = P(0) * P(0);
        const double p2_2 = P(1) * P(1);
        sigma_0(0) = -0.5 * (mat.b1 * p1_2 + mat.b2 * p2_2);
        sigma_0(1) = -0.5 * (mat.b1 * p2_2 + mat.b2 * p1_2);
        sigma_0(2) = -mat.b3 * P(0) * P(1);
        return sigma_0;
    }
};

class ElasticMaterial {
private:
    const MaterialConfig& material;

public:
    explicit ElasticMaterial(const MaterialConfig& config)
        : material(config) {}

    double mechanical_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const {
        return ElasticMath::compute_mec_energy(P, strain, material);
    }

    double elastic_energy(const Eigen::Matrix2d& strain) const {
        return ElasticMath::compute_elas_energy(strain, material);
    }

    double total_energy(const Eigen::Vector2d& P, const Eigen::Matrix2d& strain) const {
        return ElasticMath::compute_energy(P, strain, material);
    }

    Eigen::Matrix3d elastic_matrix() const {
        return ElasticMath::compute_elastic_matrix(material);
    }

    Eigen::Vector3d sigma_0(const Eigen::Vector2d& P) const {
        return ElasticMath::compute_sigma_0(P, material);
    }
};
