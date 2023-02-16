#include <iostream>
#include <vector>
#include <map>

#include "src/Atom.h"
#include "src/Vector.h"
#include "src/helpers.cpp"

using namespace std;

int main() {
    cout << "#################################################################################################" << endl;
    cout << "ЧАСТЬ 1. ПРЯМАЯ ЗАДАЧА." << endl;
    cout << "#################################################################################################" << endl;
    // Постановка: пусть даны параметры потенциала, необходимо вычислить характеристики исследуемого материала.

    cout << "Параметры потенциала: A1=" << A1 << ", A0=" << A0 << ", ksi=" << ksi << ", p=" << p << ", q=" << q
        << ", r0=" << r0 << ", a0=" << a_0 << endl << endl;

    double E_coh;  // Когезионная энергия системы

    // Гелиоцентрическая решетка
    vector<Atom> heliocentric_lattice;

    // Задаем решетку из 3 атомов в кристаллической решетке (кубик)
    init_heliocentric_lattice(heliocentric_lattice, "Ni", a_0);

    // Находим когезионную энергию системы
    E_coh = calculate_energy(heliocentric_lattice, a_0, MATR) / double(heliocentric_lattice.size());

    // Поиск характеристик
    double B, c11, c12, c44;
    calculate_characteristics(a_0, E_coh, B, c11, c12, c44);

    cout << "Ecoh = " << E_coh << endl;
    cout << "B = " << B << ", C11 = " << c11 << ", C12 = " << c12 << ", C44 = " << c44 << endl;

    cout << "#################################################################################################" << endl;
    cout << "ЧАСТЬ 2. ОБРАТНАЯ ЗАДАЧА." << endl;
    cout << "#################################################################################################" << endl;
    // Постановка: пусть даны характеристики материала, необходимо найти набор параметров, который позволит их получить.

    multimap<double, Vector> simplex;

    simplex = nelder_mead_method(E_coh, B, c11, c12, c44);

    printParams(simplex.begin()->second.vec);
    update_parameters(simplex.begin()->second.vec);

    double E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual;

    vector<Atom> Vect1;
    double a0_f;
    a0_f = simplex.begin()->second.vec[6];

    init_heliocentric_lattice(Vect1, "Ni", a0_f);

    E_coh_actual = calculate_energy(Vect1, a0_f, MATR) / double(Vect1.size());
    calculate_characteristics(
            a0_f, E_coh_actual, B_actual, C11_actual, C12_actual, C44_actual);

    cout << "E_coh_actual = " << E_coh_actual << endl;
    cout << "B_actual = " << B_actual << ", C11_actual = " << C11_actual << ", C12_actual = " << C12_actual
        << ", C44_actual = " << C44_actual << endl;

    return 0;
}

