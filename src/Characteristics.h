//
// Created by Никольский Владимир on 09.02.2023.
//

#ifndef PDP_PROJECT_CHARACTERISTICS_H
#define PDP_PROJECT_CHARACTERISTICS_H

#include <iostream>

using namespace std;

struct Characteristics {
    double B, C11, C12, C44;

    Characteristics() {};

    Characteristics(double _B, double _C11, double _C12, double _C44): B(_B), C11(_C11), C12(_C12), C44(_C44) {};

    void print() const {
        cout << "B: " << B << ' ';
        cout << "C11: " << C11 << ' ';
        cout << "C12: " << C12 << ' ';
        cout << "C44: " << C44 << endl;
    }

    void update(double _B, double _C11, double _C12, double _C44) {
        B = _B;
        C11 = _C11;
        C12 = _C12;
        C44 = _C44;
    }
};


#endif //PDP_PROJECT_CHARACTERISTICS_H
