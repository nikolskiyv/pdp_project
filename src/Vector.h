//
// Created by Никольский Владимир on 08.02.2023.
//

#ifndef PDP_PROJECT_VECTOR_H
#define PDP_PROJECT_VECTOR_H

#include <cmath>

class Vector {
public:
    double vec[7];

    Vector() {
        for(int i=0; i<7; i++) {
            vec[i] = 0;
        }
    }

    Vector(double init) {
        for(int i=0; i<7; i++) {
            vec[i] = init;
        }
    }

    Vector strict() {
        for(int i=0; i<7; i++) {
            vec[i] = fmax(vec[i], 0);
        }
        return *this;
    }

    Vector operator/=(const double& v) {
        for(int i=0; i<7; i++) {
            vec[i] /= v;
        }
        return *this;
    }

    Vector operator*(const double& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] * v;
        }
        return result;
    }


    Vector& operator+=(const Vector& v) {
        for(int i=0; i<7; i++) {
            vec[i] += v.vec[i];
        }
        return *this;
    }


    Vector operator-(const Vector& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] - v.vec[i];
        }
        return result;
    }

    Vector operator+(const Vector& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] + v.vec[i];
        }
        return result;
    }
};


#endif //PDP_PROJECT_VECTOR_H
