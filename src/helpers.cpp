#include <vector>
#include <cmath>
#include <iostream>
#include <functional>

#include "Atom.h"
#include "consts.h"

#include "heliocentric_lattice.cpp"

using namespace std;

double calculate_energy(const vector<Atom>& heliocentric_lattice, double a0_param, const double trans_matr[9]) {
    double E = 0.0;
    double cutoff = 1.7 * a0_param;

    /*
     * Для распараллеливания этого кода используется директива #pragma omp parallel for перед внешним циклом for,
     * который перебирает все пары атомов.
     * Эта директива позволяет распараллелить цикл между несколькими потоками.
     * При выполнении этого цикла каждый поток получает набор итераций цикла, которые должен выполнить.
     * Потоки работают независимо друг от друга и параллельно вычисляют Er и Eb для каждой пары атомов.
     * Когда поток заканчивает свою работу, значение E для этого потока добавляется в общую сумму E,
     * используя #pragma omp reduction.
     * Это позволяет предотвратить гонку за ресурсами, когда несколько потоков пытаются изменять одну и ту же
     * переменную E одновременно.
     **/
#pragma omp parallel for reduction(+:E)
    for (int i = 0; i < heliocentric_lattice.size(); i++) {
        double Er = 0.0, Eb = 0.0;
        const auto& atom_i = heliocentric_lattice[i];

        for (int j = 0; j < heliocentric_lattice.size(); j++) {
            if (i == j) {
                continue;
            }
            const auto& atom_j = heliocentric_lattice[j];
            double dx = atom_j.x - atom_i.x;
            double dy = atom_j.y - atom_i.y;
            double dz = atom_j.z - atom_i.z;

            for (int k = -1; k <= 1; k++) {
                for (int l = -1; l <= 1; l++) {
                    for (int m = -1; m <= 1; m++) {
                        double x = dx + a0_param * CONST_D * k;
                        double y = dy + a0_param * CONST_D * l;
                        double z = dz + a0_param * CONST_D * m;
                        double tmp_x = x * trans_matr[0] + y * trans_matr[1] + z * trans_matr[2];
                        double tmp_y = x * trans_matr[3] + y * trans_matr[4] + z * trans_matr[5];
                        double tmp_z = x * trans_matr[6] + y * trans_matr[7] + z * trans_matr[8];
                        double distance = sqrt(tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z);
                        if (distance < cutoff) {
                            Er += (A1 * (distance - r0) + A0) * exp(-p * (distance / r0 - 1));
                            Eb += ksi * ksi * exp(-2 * q * (distance / r0 - 1));
                        }
                    }
                }
            }
        }

        E += Er - sqrt(Eb);
    }

    return E;
}


void calculate_characteristics(double a0_param, double E_coh, double &B, double &C11, double &C12, double &C44) {
    /*
     * Функция подбора модуля векторного растяжения B и констант упругости C11/C12/C44
     **/

    double v0 = a0_param * a0_param * a0_param / 4;  // равновесный объем

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

double error(Vector a, double E_coh_target, double &B_target, double &C11_target, double &C12_target, double &C44_target) {
    vector<Atom> atoms;

    init_heliocentric_lattice(atoms, "Ni", a.vec[6]);
    update_parameters(a.vec);

    double E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual;

    // Считаем энергию, которую имеем на данный момент
    E_coh_actual = calculate_energy(atoms, a.vec[6], MATR) / double(atoms.size());

    // Находим характеристики, которые имеем на данный момент
    calculate_characteristics(a.vec[6], E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual);

    // Считаем ошибку как сумму квадратов отклонений
    double error = (E_coh_target - E_coh_actual) * (E_coh_target - E_coh_actual) +
                   (B_target - B_actual) * (B_target - B_actual) +
                   (C11_target - C11_actual) * (C11_target - C11_actual) +
                   (C12_target - C12_actual) * (C12_target - C12_actual) +
                   (C44_target - C44_actual) * (C44_target - C44_actual);

    return error;
}


multimap<double, Vector> nelder_mead_method(double E_coh, double &B, double &c11, double &c12, double &c44) {
    multimap<double, Vector> simplex;
    double f_h, f_g, f_l;
    Vector x_h, x_l;

    // Инициализация метода
    for (int atom_idx = 0; atom_idx < 8; ++atom_idx) {
        Vector tmp(1);
        if (atom_idx != 7) {
            for (auto &i: tmp.vec) {
                i = rand() % 10000 / 1000.0;
            }
        }
        simplex.emplace(error(tmp, E_coh, B, c11, c12, c44), tmp);
    }

    while (simplex.begin()->first > 10e-10) {  // Проверка сходимости. Смотрим на падение дисперсии ниже некоторого порога
        f_h = simplex.rbegin()->first;  // Наибольшее значение целевой функции
        f_g = (++simplex.rbegin())->first;  // Второе по величине значение целевой функции
        f_l = simplex.begin()->first;  // Наименьшее значение целевой функции
        x_h = simplex.rbegin()->second;
        x_l = simplex.begin()->second;

        Vector x_c;

        // Находим центр тяжести
        for (auto it = ++simplex.rbegin(); it != simplex.rend(); ++it) {
            x_c += it->second;
        }
        x_c /= 7;

        // Отразим точку x_h относительно точки x_c с коэффициентом alpha
        Vector x_r = (x_c * (1 + alpha_coefficient) - x_h * alpha_coefficient).strict();
        double f_r = error(x_r, E_coh, B, c11, c12, c44);  // Значение ЦФ в точке x_r

        if (f_r < f_l) {
            // Направление выбрано удачно, пробуем увеличить шаг
            Vector x_e = (x_c * (1 - gamma_coefficient) + x_r * gamma_coefficient).strict();
            double f_e = error(x_e, E_coh, B, c11, c12, c44);

            simplex.erase(--simplex.end());
            if (f_e < f_r) {
                // Из точек x_r и x_e выбираем наилучшую и заменяем ею x_h
                simplex.insert({f_e, x_e});
            } else {
                simplex.insert({f_r, x_r});
            }
        } else if (f_r < f_g) {
            // Новая точка улучшает ответ, заменим ею x_h
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});
        } else if (f_r < f_h) {
            // Новая точка улучшает ответ, но слабо, заменим ею x_h и проведем операцию сжатия
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});

            Vector x_s = (x_c * (1 - beta_coefficient) + x_h * beta_coefficient).strict();
            double f_s = error(x_s, E_coh, B, c11, c12, c44);
            if (f_s < f_h) {
                // Заменяем вершину x_h точкой x_s
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                // Сжимаем весь симплекс к точке с наименьшим значением x_i
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for (auto &x_i: tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * sigma_coefficient;
                    simplex.insert({error(new_vec, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        } else {
            // Новая точка не улучшает ответ, проведем операцию сжатия
            Vector x_s = x_c * (1 - beta_coefficient) + x_h * beta_coefficient;
            double f_s = error(x_s, E_coh, B, c11, c12, c44);
            if (f_s < f_h) {
                // Заменяем вершину x_h точкой x_s
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                // Сжимаем весь симплекс к точке с наименьшим значением x_i
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for (auto &x_i: tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * sigma_coefficient;
                    simplex.insert({error(new_vec, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        }
    }

    return simplex;
}