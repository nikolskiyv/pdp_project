//
// Created by Никольский Владимир on 07.12.2022.
//
#include <vector>

#include "../objects/Atom.h"

using namespace std;

void calculate_heliocentric_lattice(vector<Atom> &v, int d, const string& name, double a0_param) {
    // Расчет гелиоцентрической решетки

    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            for (int k = 0; k < d; k++) {
                v.emplace_back(i * a0_param, j * a0_param, k * a0_param, name);
                v.emplace_back(i * a0_param, (0.5 + j) * a0_param, (0.5 + k) * a0_param, name);
                v.emplace_back((0.5 + i) * a0_param, (0.5 + j) * a0_param, k * a0_param, name);
                v.emplace_back((0.5 + i) * a0_param, j * a0_param, (0.5 + k) * a0_param, name);

            }
        }
    }
}