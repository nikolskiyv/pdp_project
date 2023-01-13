#include <iostream>
#include <vector>

#include "src/objects/Vector.h"
#include "src/objects/Atom.h"
#include "src/objects/MetalValues.h"

#include "src/helpers.cpp"
#include "src/params_update.cpp"

using namespace std;

double a_0 = 4.085;
const MetalValues metal_values(4.085, -2.96, 1.08, 1.32, 0.97, 0.51);

int main() {
    // ЧАСТЬ 1. НАХОЖДЕНИЕ ЭНЕРГИИ

    double E_p;
    double E_coh;

    // Подбираем параметры a_0 и E_p
    update_params(a_0, E_p);

    // Гелиоцентрическая решетка
    vector<Atom> heliocentric_lattice;

    // Задаем решетку из 3 атомов в кристаллической решетке (кубик)
    calculate_heliocentric_lattice(heliocentric_lattice, 3, "Ag", a_0);

    // Находим когезионную энергию системы
    E_coh = energy(heliocentric_lattice, a_0, 3, D) / double(heliocentric_lattice.size());

    cout << "#################################################################################################" << endl;
    cout << "| Ni | ";
    cout << "-Ec = " << E_coh << ", a0 = " << a_0 << ", ";

    // Поиск параметров
    double B, c11, c12, c44;
    params(a_0, E_coh, B, c11, c12, c44);

    print_params(B, c11, c12, c44);
    cout << "#################################################################################################" << endl;

    // ЧАСТЬ 2. ПОДГОНКА ПАРАМЕТРОВ

    // Find A-A Parameters
    const Parameters xAA = random_parameters();
    const Parameters AAParams = HookeJeevesMethod(xAA, metal_values);
    const MetalValues AASpecs = computeCharasteristics(AAParams);

    std::cout << std::endl << "A-A parameters are found. Let's find A-B parameters." << std::endl << std::endl;

    // Find A-B Parameters
    const Parameters xAB = random_parameters(AAParams);
    const Parameters ABParams = HookeJeevesMethod(xAB, metal_values, AAParams);
    const MetalValues ABSpecs = params(AAParams, ABParams);

    std::cout << std::endl << "A-B parameters are found. Let's find B-B parameters." << std::endl << std::endl;

    // Find A-B Parameters
    const Parameters xBB = random_parameters(ABParams);
    const Parameters BBParams = HookeJeevesMethod(xAB, metal_values, AAParams, ABParams);
    const MetalValues BBSpecs = computeCharasteristics(AAParams, ABParams, BBParams);

    std::cout << std::endl;

    std::cout << "A-A parameters: " << std::endl;
    AAParams.print();
    AASpecs.print();
    std::cout << "Start Error = " << computeCharasteristics(xAA).error(metal_values) <<
              ", Final Error = " << AASpecs.error(metal_values) << std::endl;
    std::cout << std::endl;

    std::cout << "A-B parameters: " << std::endl;
    ABParams.print();
    ABSpecs.print();
    std::cout << "Start Error = " << computeCharasteristics(AAParams, xAB).error(metal_values) <<
              ", Final Error = " << ABSpecs.error(metal_values) << std::endl;
    std::cout << std::endl;

    std::cout << "B-B parameters: " << std::endl;
    BBParams.print();
    BBSpecs.print();
    std::cout << "Start Error = " << computeCharasteristics(AAParams, ABParams, xBB).error(metal_values) <<
              ", Final Error = " << BBSpecs.error(metal_values) << std::endl;
    std::cout << std::endl;

    return 0;

}

