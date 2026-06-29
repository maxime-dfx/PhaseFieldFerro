#pragma once
#include "IO/Datafile.h"

class Material {
    private :

    public :
        Material(const Datafile& config);
        double a0;
        double b1, b2, b3;
        double c1, c2, c3;
        double alpha_1, alpha_11, alpha_111, alpha_1111,
                        alpha_12, alpha_112, alpha_1112, 
                                            alpha_1122;
        double eps0, t, P0, mu_p, xi, c0, eta_k;
};
