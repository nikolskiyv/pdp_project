#include <iostream>
#include <vector>
#include "src/objects/Atom.h"
#include "src/helpers.cpp"

using namespace std;

int main() {
    double E_coh;

    // Гелиоцентрическая решетка
    vector<Atom> heliocentric_lattice;

    // Задаем решетку из 3 атомов в кристаллической решетке (кубик)
    calculate_heliocentric_lattice(heliocentric_lattice, "Ni", a_0);

    // Находим когезионную энергию системы
    E_coh = calculate_energy(heliocentric_lattice, a_0, MATR) / double(heliocentric_lattice.size());

    cout << "#################################################################################################" << endl;
    cout << "| Ag | ";
    cout << "Ecoh = " << E_coh << ", a0 = " << a_0 << ", ";

    // Поиск параметров
    double B, c11, c12, c44;
    calculate_parameters(a_0, E_coh, B, c11, c12, c44);

    print_params(B, c11, c12, c44);
    cout << "#################################################################################################" << endl;

    return 0;

}

