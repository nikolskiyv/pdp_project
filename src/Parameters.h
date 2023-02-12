//
// Created by Никольский Владимир on 09.02.2023.
//

#ifndef PDP_PROJECT_PARAMETERS_H
#define PDP_PROJECT_PARAMETERS_H

#include <iostream>

using namespace std;

struct Parameters {
    double A1, A0, ksi, p, q, r0, a0;

    Parameters(double _A1, double _A0, double _ksi, double _p, double _q, double _r0, double _a0):
        A1(_A1), A0(_A0), ksi(_ksi), p(_p), q(_q), r0(_r0), a0(_a0) {}

    void update(double p_[7]) {
        A1 = p_[0];
        A0 = p_[1];
        ksi = p_[2];
        p = p_[3];
        q = p_[4];
        r0 = p_[5];
        a0 = p_[6];
    }

    void print() {
        cout << "A1: " << A1 << ' ';
        cout << "A0: " << A0 << ' ';
        cout << "ksi: " << ksi << ' ';
        cout << "p: " << p << ' ';
        cout << "q: " << q << ' ';
        cout << "r0: " << r0 << ' ';
        cout << "a0: " << a0 << endl;
    }
};


#endif //PDP_PROJECT_PARAMETERS_H
