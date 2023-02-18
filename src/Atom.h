//
// Created by Никольский Владимир on 06.12.2022.
//

#include <string>
#include <utility>

#ifndef PDP_PROJECT_ATOM_H
#define PDP_PROJECT_ATOM_H

using namespace std;

class Atom {
public:
    double x, y, z;
    bool isDim = false;

    Atom(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}

    bool operator==(const Atom &right) const {
        return x == right.x &&
               y == right.y &&
               z == right.z;
    }

};


#endif //PDP_PROJECT_ATOM_H
