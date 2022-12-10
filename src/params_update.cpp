//
// Created by Никольский Владимир on 06.12.2022.
//
#include <vector>
#include <cmath>

#include "objects/Atom.h"

void update_params(double &a_0, double &E_coh) {
    float h = 1;
    int m = 0;
    int i = 1;
    double E_0 = 0;

    while (fabs(h) > 0.00001) {
        vector<Atom> v;
        a_0 += h;
        calculate_heliocentric_lattice(v, 3, "Ag", a_0);
        E_coh = energy(v, a_0, 3, D) / double(v.size());
        if (E_coh > E_0) {
            if (m == 1 || i == 1) {
                a_0 -= h;
                h /= 10;
                m = 0;
            } else {
                a_0 -= h;
                h = -h;
                m = 0;
            }
        } else {
            m = 1;
            E_0 = E_coh;
        }
        i++;
    }
}