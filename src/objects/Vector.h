//
// Created by Никольский Владимир on 06.12.2022.
//
#include <cmath>

#ifndef PDP_PROJECT_VECTOR_H
#define PDP_PROJECT_VECTOR_H


class Vector {
public:
    double vec[7]{};

    Vector() {
        for(double & i : vec) { i = 0; }
    }

    explicit Vector(double init) {
        for(double & i : vec) { i = init; }
    }

    Vector strict() {
        for(double & i : vec) { i = fmax(i, 0); }
        return *this;
    }

    Vector operator/=(const double& v) {
        for(double & i : vec) { i /= v; }
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
