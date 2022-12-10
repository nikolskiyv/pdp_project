#include <vector>
#include <cmath>
#include <iostream>

#include "objects/Atom.h"
#include "objects/Vector.h"

using namespace std;

double A1 = 0, A0 = 0.03726, ksi = 1.070, p = 16.999, q = 1.189, r0 = 2.492;
double  D[9] = { 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1 };



double energy(const std::vector<Atom> &heliocentric_lattice, double a0_param, int d, const double d_trans[9]) {

    double E = 0;

#pragma omp parallel for reduction(+: E)
    // heliocentric_lattice[i] и heliocentric_lattice[j] тут i-й и j-й атомы соответственно
    // Перебираем все пары, исключая одинаковые атомы в паре (i != j)

    for (int i = 0; i < heliocentric_lattice.size(); i++) {

        double Er = 0, Eb = 0;
        double cutoff = 1.7 * a0_param; // радиус, в который должны попасть атомы, энергию которых мы считаем

#pragma omp parallel for reduction(+: Er, Eb)
        for (int j = 0; j < heliocentric_lattice.size(); j++) {

            // считаем в трехмерном пространстве
            for (int dx = -1; dx < 2; dx++)
                for (int dy = -1; dy < 2; dy++)
                    for (int dz = -1; dz < 2; dz++)


                        if (i != j || dx != 0 || dy != 0 || dz != 0) {

                            // условие периодичности?
                            double tmp_x = heliocentric_lattice[j].x + a0_param * d * dx;
                            double tmp_y = heliocentric_lattice[j].y + a0_param * d * dy;
                            double tmp_z = heliocentric_lattice[j].z + a0_param * d * dz;

                            tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1]
                                    + tmp_z * d_trans[2];
                            tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4]
                                    + tmp_z * d_trans[5];
                            tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7]
                                    + tmp_z * d_trans[8];
                            tmp_x -= heliocentric_lattice[i].x * d_trans[0] + heliocentric_lattice[i].y * d_trans[1]
                                     + heliocentric_lattice[i].z * d_trans[2];
                            tmp_y -= heliocentric_lattice[i].x * d_trans[3] + heliocentric_lattice[i].y * d_trans[4]
                                     + heliocentric_lattice[i].z * d_trans[5];
                            tmp_z -= heliocentric_lattice[i].x * d_trans[6] + heliocentric_lattice[i].y * d_trans[7]
                                     + heliocentric_lattice[i].z * d_trans[8];

                            // считаем расстояние между i и j атомами
                            double distance = sqrt(tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z);
                            if (distance < cutoff) {
                                Er = Er + (A1 * (distance - r0) + A0) * exp(-p * (distance / r0 - 1));
                                Eb = Eb + ksi * ksi * exp(-2 * q * (distance / r0 - 1));
                            }
                        }
        }
        Eb = -sqrt(Eb);
        E += Er + Eb;

    }
    return E;

}

void params(double a0_param, double E_p, double &B, double &c11, double &c12, double &c44) {
    // Функция подбора модуля векторного растяжения(B) и констант упругости(c11, c12, c44)

    double v0 = a0_param * a0_param * a0_param / 4;
    double const_p = 0.8018993929636421;
    double alpha = 0.001;
    double alpha2 = 0.000001;

    double d_b_plus[9] = {1 + alpha, 0, 0,
                          0, 1 + alpha, 0,
                          0, 0, 1 + alpha };
    double d_b_minus[9] = {1 - alpha, 0, 0,
                           0, 1 - alpha, 0,
                           0, 0, 1 - alpha };

    double d_c11_plus[9] = {1 + alpha, 0, 0,
                            0, 1 + alpha, 0,
                            0, 0, 1};
    double d_c11_minus[9] = {1 - alpha, 0, 0,
                             0, 1 - alpha, 0,
                             0, 0, 1 };

    double d_c12_plus[9] = {1 + alpha, 0, 0,
                            0, 1 - alpha, 0,
                            0, 0, 1 };
    double d_c12_minus[9] = {1 - alpha, 0, 0,
                             0, 1 + alpha, 0,
                             0, 0, 1 };

    double d_c44_plus[9] = {1, alpha, 0,
                            alpha, 1, 0,
                            0, 0, 1 / (1 - alpha2)} ;
    double d_c44_minus[9] = {1, -alpha, 0,
                             -alpha, 1, 0,
                             0, 0, 1 / (1 - alpha2) };

    vector<Atom> Vect_p, Vect_B_plus, Vect_B_minus, Vect_c11_plus, Vect_c11_minus,
            Vect_c12_plus, Vect_c12_minus, Vect_c44_plus, Vect_c44_minus;


    calculate_heliocentric_lattice(Vect_p, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_B_plus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_B_minus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c11_plus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c11_minus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c12_plus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c12_minus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c44_plus, 3, "Ag", a0_param);
    calculate_heliocentric_lattice(Vect_c44_minus, 3, "Ag", a0_param);


    double E_B_plus = energy(Vect_B_plus, a0_param, 3, d_b_plus) / double(Vect_p.size());
    double E_B_minus = energy(Vect_B_minus, a0_param, 3, d_b_minus) / double(Vect_p.size());
    double E_c11_plus = energy(Vect_c11_plus, a0_param, 3, d_c11_plus) / double(Vect_p.size());
    double E_c11_minus = energy(Vect_c11_minus, a0_param, 3, d_c11_minus) / double(Vect_p.size());
    double E_c12_plus = energy(Vect_c12_plus, a0_param, 3, d_c12_plus) / double(Vect_p.size());
    double E_c12_minus = energy(Vect_c12_minus, a0_param, 3, d_c12_minus) / double(Vect_p.size());
    double E_c44_plus = energy(Vect_c44_plus, a0_param, 3, d_c44_plus) / double(Vect_p.size());
    double E_c44_minus = energy(Vect_c44_minus, a0_param, 3, d_c44_minus) / double(Vect_p.size());
    double d2_E_B = (E_B_plus - 2 * E_p + E_B_minus) / alpha2;
    double d2_E_c11 = (E_c11_plus - 2 * E_p + E_c11_minus) / alpha2;
    double d2_E_c12 = (E_c12_plus - 2 * E_p + E_c12_minus) / alpha2;
    double d2_E_c44 = (E_c44_plus - 2 * E_p + E_c44_minus) / alpha2;

    B = (d2_E_B * 2 * const_p) / (9.0 * v0);
    c11 = ((d2_E_c11 + d2_E_c12) * const_p) / (2.0 * v0);
    c12 = ((d2_E_c11 - d2_E_c12) * const_p) / (2.0 * v0);
    c44 = (d2_E_c44 * const_p) / (2.0 * v0);
}

