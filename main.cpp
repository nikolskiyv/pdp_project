#include <iostream>
#include <vector>
#include <map>

#include "src/Atom.h"
#include "src/Vector.h"
#include "src/helpers.cpp"

using namespace std;

int main() {
    clock_t tStart = clock();

    cout << endl << endl;
    cout << "#################################################################################################" << endl;
    cout << "------------------------------------[ЧАСТЬ 1. ПРЯМАЯ ЗАДАЧА.]------------------------------------" << endl;
    cout << "#################################################################################################" << endl;
    // Постановка: пусть даны параметры потенциала, необходимо вычислить характеристики исследуемого материала.

    cout << "Параметры потенциала: A1=" << A1 << ", A0=" << A0 << ", ksi=" << ksi << ", p=" << p << ", q=" << q
        << ", r0=" << r0 << ", a0=" << a_0 << endl << endl;

    double E_coh;  // Когезионная энергия системы

    // Гелиоцентрическая решетка
    vector<Atom> heliocentric_lattice;

    // Задаем решетку из 3 атомов в кристаллической решетке (кубик)
    init_heliocentric_lattice(heliocentric_lattice, a_0);

    // Находим когезионную энергию системы
    E_coh = calculate_full_energy(heliocentric_lattice, a_0, MATR) / double(heliocentric_lattice.size());
/*
 * Отрисовка графика
 *
    for (int a0 = 30; a0 < 100; a0++)
    {
        double a0_iter = a0 / 10.0;

        heliocentric_lattice.clear();
        init_heliocentric_lattice(heliocentric_lattice, a0_iter);
        E_coh = calculate_sol_energy(heliocentric_lattice, a0_iter, MATR);
        cout << a0_iter << " " << E_coh << endl;
    }

    cout << "ECOH" << E_coh << endl;
*/

    // Поиск характеристик
    double B, c11, c12, c44;
    calculate_characteristics(a_0, E_coh, B, c11, c12, c44);

    cout << "----------------------------------------------[Ag]-----------------------------------------------" << endl;
    cout << "Ecoh = " << E_coh << endl;
    cout << "B = " << B << ", C11 = " << c11 << ", C12 = " << c12 << ", C44 = " << c44 << endl;
    cout << "-------------------------------------------------------------------------------------------------" << endl;

    cout << endl << endl;
    cout << "#################################################################################################" << endl;
    cout << "-----------------------------------[ЧАСТЬ 2. ОБРАТНАЯ ЗАДАЧА.]-----------------------------------" << endl;
    cout << "#################################################################################################" << endl;
    // Постановка: пусть даны характеристики материала, необходимо найти набор параметров, который позволит их получить.

    cout << endl;

    cout << "----------------------------------------------[Ag-Ag]-----------------------------------------------" << endl;
    multimap<double, Vector> simplexBB;

    simplexBB = nelder_mead_method(E_coh, B, c11, c12, c44, '1');

    printParams(simplexBB.begin()->second.vec);
    update_parameters(simplexBB.begin()->second.vec);

    double E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual;

    vector<Atom> Vect1;
    double a0_f;
    a0_f = simplexBB.begin()->second.vec[6];

    init_heliocentric_lattice(Vect1, a0_f);

    E_coh_actual = calculate_full_energy(Vect1, a0_f, MATR) / double(Vect1.size());
    calculate_characteristics(
            a0_f, E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual);

    cout << "E_coh_actual = " << E_coh_actual << endl;
    cout << "B_actual = " << B_actual << ", C11_actual = " << C11_actual << ", C12_actual = " << C12_actual
        << ", C44_actual = " << C44_actual << endl;
    cout << "-------------------------------------------------------------------------------------------------" << endl;

    cout << endl;

    cout << "----------------------------------------------[Ag-Ni]--------------------------------------------" << endl;

    multimap<double, Vector> simplexAB;

    B = B_true;
    c11 = C11_true;
    c12 = C12_true;
    c44 = C44_true;

    simplexAB = nelder_mead_method(E_sol_true, B, c11, c12, c44, '2');
    cout << endl;

    printParams(simplexAB.begin()->second.vec);
    update_parameters(simplexAB.begin()->second.vec);

    vector<Atom> Vect1_AB;
    double a0_f_AB;
    a0_f_AB = simplexAB.begin()->second.vec[6];

    init_heliocentric_lattice(Vect1_AB, a0_f_AB);

    double E_sol_actual = calculate_sol_energy(Vect1_AB, a0_f_AB, MATR);

    cout << "E_sol_actual = " << E_sol_actual << endl;

    cout << "-------------------------------------------------------------------------------------------------" << endl;

    cout << endl;

    printf("Execution time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

