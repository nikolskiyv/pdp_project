//
// Created by Никольский Владимир on 14.01.2023.
//

#ifndef PDP_PROJECT_CONSTS_H
#define PDP_PROJECT_CONSTS_H

double alpha_1 = 0.001;
double alpha_2 = 0.000001;

double const alpha_coefficient = 1;  // Коэффициент отражения
double const beta_coefficient = 0.5;  // Коэффициент сжатия
double const gamma_coefficient = 2;  // Коэффициент растяжения
double const sigma_coefficient = 0.5;  // Коэффициент глобального сжатия

const int CONST_D = 3;  // В базисе 4 атоме, в большой решетке 27. Всего 108. Итерируемся по
const double  MATR[9] = {1, 0, 0,
                         0, 1, 0,
                         0, 0, 1 };

double A1 = 0, A0 = 0.1028, ksi = 1.178, p = 10.928, q = 3.139, r0 = 2.889, a_0 = 4.085;

void update_parameters(const double new_parameters[7]) {
    A1 = new_parameters[0];
    A0 = new_parameters[1];
    ksi = new_parameters[2];
    p = new_parameters[3];
    q = new_parameters[4];
    r0 = new_parameters[5];
    a_0 = new_parameters[6];
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
