//
// Created by Никольский Владимир on 13.01.2023.
//

#include <iostream>

#ifndef PDP_PROJECT_PARAMETERS_H
#define PDP_PROJECT_PARAMETERS_H

using namespace std;

struct Parameters
{
    Parameters(double i_r0 = 0.0,
               double i_A0 = 0.0,
               double i_A1 = 0.0,
               double i_p = 0.0,
               double i_q = 0.0,
               double i_ksi = 0.0) :
            r0(i_r0), A0(i_A0), A1(i_A1),
            p(i_p), q(i_q), ksi(i_ksi)
    {
    }

    Parameters(const Parameters& i_o) :
            r0(i_o.r0), A0(i_o.A0), A1(i_o.A1),
            p(i_o.p), q(i_o.q), ksi(i_o.ksi)
    {
    }

    double operator[] (int i_num) const
    {
        switch(i_num)
        {
            case 0:
                return r0;
            case 1:
                return A0;
            case 2:
                return A1;
            case 3:
                return p;
            case 4:
                return q;
            case 5:
                return ksi;
            default:
                throw runtime_error("Incorrect index");
        }
    }

    bool operator== (const Parameters& i_o) const
    {
        return r0 == i_o[0] && A0 == i_o[1] && A1 == i_o[2] &&
               p == i_o[3] && q == i_o[4] && ksi == i_o[5];
    }

    // Use only for comparing steps and limits!!!!
    bool operator> (const Parameters& i_o) const
    {
        return r0 > i_o[0] || A0 > i_o[1] || A1 > i_o[2] ||
               p > i_o[3] || q > i_o[4] || ksi > i_o[5];
    }


    void setValue(double i_value, int i_id)
    {
        switch(i_id)
        {
            case 0:
                r0 = i_value;
                return;
            case 1:
                A0 = i_value;
                return;
            case 2:
                A1 = i_value;
                return;
            case 3:
                p = i_value;
                return;
            case 4:
                q = i_value;
                return;
            case 5:
                ksi = i_value;
                return;
            default:
                return;
        }
    }

    void print() const
    {
        cout << "r0 = " << r0 << ", A0 = "<< A0 <<
                  ", A1 = " << A1 << ", p = " << p << ", q = " << q <<
                  ", ksi = " << ksi << endl;
    }

    double r0;
    double A0;
    double A1;
    double p;
    double q;
    double ksi;
};


#endif //PDP_PROJECT_PARAMETERS_H
