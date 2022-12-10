#include <iostream>
#include <vector>
#include <fstream>

#include "src/objects/Vector.h"
#include "src/objects/Atom.h"
#include "src/helpers.cpp"
#include "src/params_update.cpp"
#include "src/calculation/heliocentric_lattice.cpp"

using namespace std;


int main() {
    double E_p;
    double E_coh;

    double a_0 = 4.085;  // Параметр решетки

    // Подбираем параметры a_0 и E_p
    update_params(a_0, E_p);

    vector<Atom> heliocentric_lattice;  // Гелиоцентрическая решетка

    // Задаем решетку из 3 атомов в кристаллической решетке (кубик)
    calculate_heliocentric_lattice(heliocentric_lattice, 3, "Ag", a_0);

    // Находим когезионную энергию системы
    E_coh = energy(heliocentric_lattice, a_0, 3, D) / double(heliocentric_lattice.size());

    ofstream e_coh_out("E_coh_result.txt");
    e_coh_out << a_0 << " " << E_coh << endl;
    e_coh_out.close();

    cout << "a0 = " << a_0 << endl << "E_coh = " << E_coh << endl;

    // Поиск параметров
    double B, c11, c12, c44;
    params(a_0, E_coh, B, c11, c12, c44);
    cout << "B = " << B << endl << "c11 = " << c11 << endl
              << "c12 = " << c12 << endl << "c44 = " << c44 << endl;

    return 0;

}

