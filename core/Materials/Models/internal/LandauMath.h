#pragma once

#include <Eigen/Core>
#include "IO/include/ConfigTypes.h"

// LandauMath
// -----------------------------------------------------------------------
// Toolbox mathématique pure pour l'énergie de Landau-Devonshire et ses
// dérivées (gradient/hessienne par rapport à P). Anciennement dupliquée
// entre internal/LandauDevonshire.h (jamais utilisée) et une copie inline
// dans internal/FerroelectricCore.h (effectivement utilisée par
// Models/Ferroelectric.h) : les deux implémentations étaient
// mathématiquement identiques, on les fusionne ici en une seule source de
// vérité.
struct LandauTerms {
    double energy = 0.0;
    Eigen::Vector2d dchi_dP = Eigen::Vector2d::Zero();
    Eigen::Matrix2d d2chi_dP2 = Eigen::Matrix2d::Zero();
};

class LandauMath {
public:
    static LandauTerms compute(const Eigen::Vector2d& P, const MaterialConfig& mat) {
        LandauTerms terms;

        const double a1 = mat.alpha_1, a11 = mat.alpha_11, a12 = mat.alpha_12;
        const double a111 = mat.alpha_111, a112 = mat.alpha_112;
        const double a1111 = mat.alpha_1111, a1112 = mat.alpha_1112, a1122 = mat.alpha_1122;

        const double p1 = P(0), p2 = P(1);
        const double p1_2 = p1 * p1, p1_3 = p1_2 * p1, p1_4 = p1_2 * p1_2, p1_5 = p1_4 * p1, p1_6 = p1_4 * p1_2, p1_7 = p1_6 * p1, p1_8 = p1_4 * p1_4;
        const double p2_2 = p2 * p2, p2_3 = p2_2 * p2, p2_4 = p2_2 * p2_2, p2_5 = p2_4 * p2, p2_6 = p2_4 * p2_2, p2_7 = p2_6 * p2, p2_8 = p2_4 * p2_4;

        // 1. Énergie
        terms.energy = mat.alpha_1    * (p1_2 + p2_2) +
                       mat.alpha_11   * (p1_4 + p2_4) +
                       mat.alpha_12   * (p1_2 * p2_2) +
                       mat.alpha_111  * (p1_6 + p2_6) +
                       mat.alpha_112  * (p1_2 * p2_4 + p2_2 * p1_4) +
                       mat.alpha_1111 * (p1_8 + p2_8) +
                       mat.alpha_1112 * (p1_6 * p2_2 + p2_6 * p1_2) +
                       mat.alpha_1122 * (p1_4 * p2_4);

        // 2. Gradient dchi / dP
        terms.dchi_dP(0) = 2.0*a1*p1 + 4.0*a11*p1_3 + 2.0*a12*p1*p2_2 + 6.0*a111*p1_5 + 4.0*a112*p1_3*p2_2 + 2.0*a112*p1*p2_4 + 8.0*a1111*p1_7 + 6.0*a1112*p1_5*p2_2 + 2.0*a1112*p1*p2_6 + 4.0*a1122*p1_3*p2_4;
        terms.dchi_dP(1) = 2.0*a1*p2 + 4.0*a11*p2_3 + 2.0*a12*p2*p1_2 + 6.0*a111*p2_5 + 4.0*a112*p2_3*p1_2 + 2.0*a112*p2*p1_4 + 8.0*a1111*p2_7 + 6.0*a1112*p2_5*p1_2 + 2.0*a1112*p2*p1_6 + 4.0*a1122*p2_3*p1_4;

        // 3. Hessienne d2chi / dP2
        terms.d2chi_dP2(0,0) = 2.0*a1 + 12.0*a11*p1_2 + 2.0*a12*p2_2 + 30.0*a111*p1_4 + 12.0*a112*p1_2*p2_2 + 2.0*a112*p2_4 + 56.0*a1111*p1_6 + 30.0*a1112*p1_4*p2_2 + 2.0*a1112*p2_6 + 12.0*a1122*p1_2*p2_4;
        terms.d2chi_dP2(1,1) = 2.0*a1 + 12.0*a11*p2_2 + 2.0*a12*p1_2 + 30.0*a111*p2_4 + 12.0*a112*p2_2*p1_2 + 2.0*a112*p1_4 + 56.0*a1111*p2_6 + 30.0*a1112*p2_4*p1_2 + 2.0*a1112*p1_6 + 12.0*a1122*p2_2*p1_4;
        terms.d2chi_dP2(0,1) = 4.0*a12*p1*p2 + 8.0*a112*p1*p2*(p1_2 + p2_2) + 12.0*a1112*p1*p2*(p1_4 + p2_4) + 16.0*a1122*p1_3*p2_3;
        terms.d2chi_dP2(1,0) = terms.d2chi_dP2(0,1);

        return terms;
    }
};
