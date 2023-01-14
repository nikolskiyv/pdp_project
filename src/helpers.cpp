#include <vector>
#include <cmath>
#include <iostream>

#include "objects/Atom.h"
#include "consts.h"

#include "calculation/heliocentric_lattice.cpp"

using namespace std;

double calculate_energy(const vector<Atom> &heliocentric_lattice, double a0_param, const double d_trans[9]) {

    double E = 0;

#pragma omp parallel for reduction(+: E)
    /*
     * heliocentric_lattice[i] и heliocentric_lattice[j] тут i-й и j-й атомы соответственно
     * Перебираем все пары, исключая одинаковые атомы в паре (i != j)
     **/
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
                            // условие периодичности
                            double tmp_x = heliocentric_lattice[j].x + a0_param * CONST_D * dx;
                            double tmp_y = heliocentric_lattice[j].y + a0_param * CONST_D * dy;
                            double tmp_z = heliocentric_lattice[j].z + a0_param * CONST_D * dz;

                            tmp_x = tmp_x * d_trans[0] + tmp_y * d_trans[1] + tmp_z * d_trans[2];
                            tmp_y = tmp_x * d_trans[3] + tmp_y * d_trans[4] + tmp_z * d_trans[5];
                            tmp_z = tmp_x * d_trans[6] + tmp_y * d_trans[7] + tmp_z * d_trans[8];
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

void calculate_parameters(double a0_param, double E_coh, double &B, double &C11, double &C12, double &C44) {
    /*
     * Функция подбора модуля векторного растяжения(B) и констант упругости(C11, C12, C44)
     **/

    double v0 = a0_param * a0_param * a0_param / 4;
    double const_p = 0.8018993929636421;
    double alpha = 0.001;
    double alpha2 = 0.000001;

    double D_B_pos[9] = {1 + alpha, 0, 0,
                         0, 1 + alpha, 0,
                         0, 0, 1 + alpha };
    double D_B_neg[9] = {1 - alpha, 0, 0,
                         0, 1 - alpha, 0,
                         0, 0, 1 - alpha };

    double D_C11_pos[9] = {1 + alpha, 0, 0,
                           0, 1 + alpha, 0,
                           0, 0, 1};
    double D_C11_neg[9] = {1 - alpha, 0, 0,
                           0, 1 - alpha, 0,
                           0, 0, 1 };

    double D_C12_pos[9] = {1 + alpha, 0, 0,
                           0, 1 - alpha, 0,
                           0, 0, 1 };
    double D_C12_neg[9] = {1 - alpha, 0, 0,
                           0, 1 + alpha, 0,
                           0, 0, 1 };

    double D_C44_pos[9] = {1, alpha, 0,
                           alpha, 1, 0,
                           0, 0, 1 / (1 - alpha2)} ;
    double D_C44_neg[9] = {1, -alpha, 0,
                           -alpha, 1, 0,
                           0, 0, 1 / (1 - alpha2) };

    vector<Atom> atoms_p, atoms_B_pos, atoms_B_neg, atoms_C11_pos, atoms_C11_neg, atoms_c12_pos, atoms_C12_neg,
        atoms_C44_pos, atoms_C44_neg;


    calculate_heliocentric_lattice(atoms_p, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_B_pos, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_B_neg, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_C11_pos, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_C11_neg, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_c12_pos, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_C12_neg, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_C44_pos, "Ni", a0_param);
    calculate_heliocentric_lattice(atoms_C44_neg, "Ni", a0_param);


    double E_B_pos = calculate_energy(atoms_B_pos, a0_param, D_B_pos) / double(atoms_p.size());
    double E_B_neg = calculate_energy(atoms_B_neg, a0_param, D_B_neg) / double(atoms_p.size());
    double E_C11_pos = calculate_energy(atoms_C11_pos, a0_param, D_C11_pos) / double(atoms_p.size());
    double E_C11_neg = calculate_energy(atoms_C11_neg, a0_param, D_C11_neg) / double(atoms_p.size());
    double E_C12_pos = calculate_energy(atoms_c12_pos, a0_param, D_C12_pos) / double(atoms_p.size());
    double E_C12_neg = calculate_energy(atoms_C12_neg, a0_param, D_C12_neg) / double(atoms_p.size());
    double E_C44_pos = calculate_energy(atoms_C44_pos, a0_param, D_C44_pos) / double(atoms_p.size());
    double E_C44_neg = calculate_energy(atoms_C44_neg, a0_param, D_C44_neg) / double(atoms_p.size());

    double d2_E_B = (E_B_pos - 2 * E_coh + E_B_neg) / alpha2;
    double d2_E_C11 = (E_C11_pos - 2 * E_coh + E_C11_neg) / alpha2;
    double d2_E_C12 = (E_C12_pos - 2 * E_coh + E_C12_neg) / alpha2;
    double d2_E_C44 = (E_C44_pos - 2 * E_coh + E_C44_neg) / alpha2;

    B = (d2_E_B * 2 * const_p) / (9.0 * v0);
    C11 = ((d2_E_C11 + d2_E_C12) * const_p) / (2.0 * v0);
    C12 = ((d2_E_C11 - d2_E_C12) * const_p) / (2.0 * v0);
    C44 = (d2_E_C44 * const_p) / (2.0 * v0);
}

void print_params(double B, double C11, double C12, double C44) {
    cout << "B = " << B << ", C11 = " << C11 << ", C12 = " << C12 << ", C44 = " << C44 <<endl;
}
