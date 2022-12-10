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

    string name;
    double x, y, z;

    Atom(double _x, double _y, double _z, string _name)
            : x(_x), y(_y), z(_z), name(std::move(_name)) {}

};


#endif //PDP_PROJECT_ATOM_H
