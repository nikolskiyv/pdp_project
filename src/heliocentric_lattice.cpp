//
// Created by Никольский Владимир on 07.12.2022.
//
#include <vector>

#include "Atom.h"

using namespace std;

#include <vector>
#include <string>
#include "Atom.h"

using namespace std;

void init_heliocentric_lattice(vector<Atom> &atoms, const string& name, double a0_param) {
    /*
    * Функция расставляет атомы в узлах кристаллической решетки
    */

    for (int i = 0; i < CONST_D; i++) {
        for (int j = 0; j < CONST_D; j++) {
            for (int k = 0; k < CONST_D; k++) {
                Atom atom1(i * a0_param, j * a0_param, k * a0_param, name);
                atoms.push_back(atom1);

                Atom atom2(i * a0_param, (0.5 + j) * a0_param, (0.5 + k) * a0_param, name);
                atoms.push_back(atom2);

                Atom atom3((0.5 + i) * a0_param, (0.5 + j) * a0_param, k * a0_param, name);
                atoms.push_back(atom3);

                Atom atom4((0.5 + i) * a0_param, j * a0_param, (0.5 + k) * a0_param, name);
                atoms.push_back(atom4);
            }
        }
    }
}