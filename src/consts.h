//
// Created by Никольский Владимир on 14.01.2023.
//

#ifndef PDP_PROJECT_CONSTS_H
#define PDP_PROJECT_CONSTS_H

double const_p = 0.8018993929636421;
double alpha_1 = 0.001;
double alpha_2 = 0.000001;

double const alpha_coefficient = 1;  // Коэффициент отражения
double const beta_coefficient = 0.5;  // Коэффициент сжатия
double const gamma_coefficient = 2;  // Коэффициент растяжения
double const sigma_coefficient = 0.5;  // Коэффициент глобального сжатия

const int CONST_D = 3;
const double  MATR[9] = {1, 0, 0,
                         0, 1, 0,
                         0, 0, 1 };

double A1 = 0, A0 = 0.1028, ksi = 1.178, p = 10.928, q = 3.139, r0 = 2.889, a_0 = 4.085;

void update_parameters(double p_[7]) {
    A1 = p_[0];
    A0 = p_[1];
    ksi = p_[2] ;
    p = p_[3] ;
    q = p_[4] ;
    r0 = p_[5] ;
    a_0 = p_[6];
}

void printParams(double p_[7]) {
    cout << "A1: " << p_[0] << ' ';
    cout << "A0: " << p_[1] << ' ';
    cout << "ksi: " << p_[2] << ' ';
    cout << "p: " << p_[3] << ' ';
    cout << "q: " << p_[4] << ' ';
    cout << "r0: " << p_[5] << ' ';
    cout << "a0: " << p_[6] << endl;
}

#endif //PDP_PROJECT_CONSTS_H
