//
// Created by Никольский Владимир on 13.01.2023.
//

#include <iostream>
#include <cmath>

#ifndef PDP_PROJECT_METALVALUES_H
#define PDP_PROJECT_METALVALUES_H

using namespace std;

struct MetalValues {
    MetalValues(double i_a = 0.0,
                   double i_Ec = 0.0,
                   double i_B = 0.0,
                   double i_C11 = 0.0,
                   double i_C12 = 0.0,
                   double i_C44 = 0.0) :
            a0(i_a), Ec(i_Ec), B(i_B),
            C11(i_C11), C12(i_C12), C44(i_C44)
    {
    }

    double error(const MetalValues i_tableSpecs) const
    {
        const auto sqr = [](const double i_value, const double i_tableValue)
        {
            return (i_value - i_tableValue) * (i_value - i_tableValue) / (i_tableValue * i_tableValue);
        };

        return sqrt(
                (sqr(a0, i_tableSpecs.a0)
                + sqr(Ec, i_tableSpecs.Ec)
                + sqr(B, i_tableSpecs.B)
                + sqr(C11, i_tableSpecs.C11)
                + sqr(C12, i_tableSpecs.C12)
                + sqr(C44, i_tableSpecs.C44)
                ) / 6
                );
    }

    void print() const
    {
        cout << "a = " << a0 << ", Ec = "<< Ec <<
                  ", B = " << B << ", C11 = " << C11 << ", C12 = " << C12 <<
                  ", C44 = " << C44 << endl;
    }

    double a0;
    double Ec;
    double B;
    double C11;
    double C12;
    double C44;

};


#endif //PDP_PROJECT_METALVALUES_H
