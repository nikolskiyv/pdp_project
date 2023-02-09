#include <vector>
#include <cmath>
#include <iostream>

#include "Atom.h"
#include "consts.h"

#include "heliocentric_lattice.cpp"

using namespace std;

double calculate_energy(const vector<Atom> &heliocentric_lattice, double a0_param, const double trans_matr[9]) {
    /*
     *
     */

    double E = 0;

    // heliocentric_lattice[i] и heliocentric_lattice[j] тут i-й и j-й атомы соответственно
    // Перебираем все пары, исключая одинаковые атомы в паре (i != j)
// #pragma omp parallel for reduction(+: E)
    for (int i = 0; i < heliocentric_lattice.size(); i++) {
        double Er = 0, Eb = 0;
        double cutoff = 1.7 * a0_param;  // Радиус, в который должны попасть атомы, энергию которых мы считаем

// #pragma omp parallel for reduction(+: Er, Eb)
        for (int j = 0; j < heliocentric_lattice.size(); j++) {

            // считаем в трехмерном пространстве
            for (int dx = -1; dx < 2; dx++)
                for (int dy = -1; dy < 2; dy++)
                    for (int dz = -1; dz < 2; dz++)
                        if (i != j || dx != 0 || dy != 0 || dz != 0) {

                            // Условие периодичности
                            double tmp_x = heliocentric_lattice[j].x + a0_param * CONST_D * dx;
                            double tmp_y = heliocentric_lattice[j].y + a0_param * CONST_D * dy;
                            double tmp_z = heliocentric_lattice[j].z + a0_param * CONST_D * dz;

                            tmp_x = tmp_x * trans_matr[0] + tmp_y * trans_matr[1] + tmp_z * trans_matr[2];
                            tmp_y = tmp_x * trans_matr[3] + tmp_y * trans_matr[4] + tmp_z * trans_matr[5];
                            tmp_z = tmp_x * trans_matr[6] + tmp_y * trans_matr[7] + tmp_z * trans_matr[8];
                            tmp_x -= heliocentric_lattice[i].x * trans_matr[0] + heliocentric_lattice[i].y * trans_matr[1]
                                     + heliocentric_lattice[i].z * trans_matr[2];
                            tmp_y -= heliocentric_lattice[i].x * trans_matr[3] + heliocentric_lattice[i].y * trans_matr[4]
                                     + heliocentric_lattice[i].z * trans_matr[5];
                            tmp_z -= heliocentric_lattice[i].x * trans_matr[6] + heliocentric_lattice[i].y * trans_matr[7]
                                     + heliocentric_lattice[i].z * trans_matr[8];

                            // Считаем расстояние между i и j атомами
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

void calculate_characteristics(double a0_param, double E_coh, double &B, double &C11, double &C12, double &C44) {
    /*
     * Функция подбора модуля векторного растяжения B и констант упругости C11/C12/C44
     **/

    double v0 = a0_param * a0_param * a0_param / 4;

    double D_B_pos[9] = {1 + alpha_1, 0, 0,
                         0, 1 + alpha_1, 0,
                         0, 0, 1 + alpha_1 };

    double D_B_neg[9] = {1 - alpha_1, 0, 0,
                         0, 1 - alpha_1, 0,
                         0, 0, 1 - alpha_1 };

    double D_C11_pos[9] = {1 + alpha_1, 0, 0,
                           0, 1 + alpha_1, 0,
                           0, 0, 1};

    double D_C11_neg[9] = {1 - alpha_1, 0, 0,
                           0, 1 - alpha_1, 0,
                           0, 0, 1 };

    double D_C12_pos[9] = {1 + alpha_1, 0, 0,
                           0, 1 - alpha_1, 0,
                           0, 0, 1 };

    double D_C12_neg[9] = {1 - alpha_1, 0, 0,
                           0, 1 + alpha_1, 0,
                           0, 0, 1 };

    double D_C44_pos[9] = {1, alpha_1, 0,
                           alpha_1, 1, 0,
                           0, 0, 1 / (1 - alpha_2) };

    double D_C44_neg[9] = {1, -alpha_1, 0,
                           -alpha_1, 1, 0,
                           0, 0, 1 / (1 - alpha_2) };

    vector<Atom> atoms_p, atoms_B_pos, atoms_B_neg, atoms_C11_pos, atoms_C11_neg, atoms_c12_pos, atoms_C12_neg,
        atoms_C44_pos, atoms_C44_neg;


    init_heliocentric_lattice(atoms_p, "Ni", a0_param);
    init_heliocentric_lattice(atoms_B_pos, "Ni", a0_param);
    init_heliocentric_lattice(atoms_B_neg, "Ni", a0_param);
    init_heliocentric_lattice(atoms_C11_pos, "Ni", a0_param);
    init_heliocentric_lattice(atoms_C11_neg, "Ni", a0_param);
    init_heliocentric_lattice(atoms_c12_pos, "Ni", a0_param);
    init_heliocentric_lattice(atoms_C12_neg, "Ni", a0_param);
    init_heliocentric_lattice(atoms_C44_pos, "Ni", a0_param);
    init_heliocentric_lattice(atoms_C44_neg, "Ni", a0_param);


    double E_B_pos = calculate_energy(atoms_B_pos, a0_param, D_B_pos) / double(atoms_p.size());
    double E_B_neg = calculate_energy(atoms_B_neg, a0_param, D_B_neg) / double(atoms_p.size());
    double E_C11_pos = calculate_energy(atoms_C11_pos, a0_param, D_C11_pos) / double(atoms_p.size());
    double E_C11_neg = calculate_energy(atoms_C11_neg, a0_param, D_C11_neg) / double(atoms_p.size());
    double E_C12_pos = calculate_energy(atoms_c12_pos, a0_param, D_C12_pos) / double(atoms_p.size());
    double E_C12_neg = calculate_energy(atoms_C12_neg, a0_param, D_C12_neg) / double(atoms_p.size());
    double E_C44_pos = calculate_energy(atoms_C44_pos, a0_param, D_C44_pos) / double(atoms_p.size());
    double E_C44_neg = calculate_energy(atoms_C44_neg, a0_param, D_C44_neg) / double(atoms_p.size());

    double d2_E_B = (E_B_pos - 2 * E_coh + E_B_neg) / alpha_2;
    double d2_E_C11 = (E_C11_pos - 2 * E_coh + E_C11_neg) / alpha_2;
    double d2_E_C12 = (E_C12_pos - 2 * E_coh + E_C12_neg) / alpha_2;
    double d2_E_C44 = (E_C44_pos - 2 * E_coh + E_C44_neg) / alpha_2;

    B = (d2_E_B * 2 * const_p) / (9.0 * v0);
    C11 = ((d2_E_C11 + d2_E_C12) * const_p) / (2.0 * v0);
    C12 = ((d2_E_C11 - d2_E_C12) * const_p) / (2.0 * v0);
    C44 = (d2_E_C44 * const_p) / (2.0 * v0);
}

double error(Vector a, double E_coh,   double &B, double &C11, double &C12, double &C44) {
    vector<Atom> Vect;

    init_heliocentric_lattice(Vect, "Ag", a.vec[6]);
    update_parameters(a.vec);
    //printParams(a.vec);

    double E_coh_f, B_f, C11_f, C12_f, C44_f;

    E_coh_f = calculate_energy(Vect, a.vec[6], MATR) / double(Vect.size());

    calculate_characteristics(a.vec[6], E_coh_f, B_f, C11_f, C12_f, C44_f);

    double err =
            (E_coh - E_coh_f) * (E_coh - E_coh_f) +
            (B - B_f) * (B - B_f) +
            (C11 - C11_f) * (C11 - C11_f) +
            (C12 - C12_f) * (C12 - C12_f) +
            (C44 - C44_f) * (C44 - C44_f);

    return err;
}

multimap<double, Vector> nelder_mead_method(double E_coh, double &B, double &c11, double &c12, double &c44) {
    multimap<double, Vector> simplex;
    double f_h;  // Наибольшее значение целевой функции
    double f_g;  // Второе по величине значение целевой функции
    double f_l;  // Наименьшее значение целевой функции
    Vector x_h, x_l;

    for(int atom_idx=0; atom_idx<8; atom_idx++) {
        Vector tmp(1);
        if(atom_idx != 7) {
            for(double & i : tmp.vec) {
                i = rand() % 10000 / 1000.0;
            }
        }
        simplex.insert({error(tmp, E_coh, B, c11, c12, c44), tmp});
    }

    while(simplex.begin()->first > 10e-10) {
        f_h = simplex.rbegin()->first;
        f_g = (++simplex.rbegin())->first;
        f_l = simplex.begin()->first;

        x_h = simplex.rbegin()->second;
        x_l = simplex.begin()->second;

        // cout << "BEST ERROR: " << f_l << endl;
        // cout << "WORST ERROR: " << f_h << endl;
        // printParams(x_h.vec);

        Vector x_c;

        for(auto it = ++simplex.rbegin(); it != simplex.rend(); ++it) {
            x_c += it->second;
        }
        x_c /= 7;
        Vector x_r = (x_c * (1 + alpha_coefficient) - x_h * alpha_coefficient).strict();
        double f_r = error(x_r, E_coh, B, c11, c12, c44);

        if (f_r < f_l) {
            Vector x_e = (x_c * (1 - gamma_coefficient) + x_r * gamma_coefficient).strict();
            double f_e = error(x_e, E_coh, B, c11, c12, c44);

            simplex.erase(--simplex.end());
            if (f_e < f_r) {
                simplex.insert({f_e, x_e});
            } else {
                simplex.insert({f_r, x_r});
            }
        } else if (f_r < f_g) {
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});
        } else if (f_r < f_h) {
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});

            Vector x_s = (x_c * (1 - beta_coefficient) + x_h * beta_coefficient).strict();
            double f_s = error(x_s, E_coh, B, c11, c12, c44);
            if (f_s < f_h) {
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for (auto& x_i : tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * sigma_coefficient;
                    simplex.insert({error(new_vec, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        } else {
            Vector x_s = x_c * (1 - beta_coefficient) + x_h * beta_coefficient;
            double f_s = error(x_s, E_coh, B, c11, c12, c44);
            if (f_s < f_h) {
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for (auto& x_i : tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * sigma_coefficient;
                    simplex.insert({error(new_vec, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        }
    }

    return simplex;
}