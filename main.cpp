#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>


#include "src/objects/Atom.h"
#include "src/helpers.cpp"

using namespace std;

class Vector {
public:
    double vec[7];

    Vector() {
        for(int i=0; i<7; i++) {
            vec[i] = 0;
        }
    }

    Vector(double init) {
        for(int i=0; i<7; i++) {
            vec[i] = init;
        }
    }

    Vector strict() {
        for(int i=0; i<7; i++) {
            vec[i] = fmax(vec[i], 0);
        }
        return *this;
    }

    Vector operator/=(const double& v) {
        for(int i=0; i<7; i++) {
            vec[i] /= v;
        }
        return *this;
    }

    Vector operator*(const double& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] * v;
        }
        return result;
    }


    Vector& operator+=(const Vector& v) {
        for(int i=0; i<7; i++) {
            vec[i] += v.vec[i];
        }
        return *this;
    }


    Vector operator-(const Vector& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] - v.vec[i];
        }
        return result;
    }

    Vector operator+(const Vector& v) const {
        Vector result;
        for(int i=0; i<7; i++) {
            result.vec[i] = vec[i] + v.vec[i];
        }
        return result;
    }
};

void passParams(double p_[7]) {
    A1 = p_[0];
    A0 = p_[1];
    ksi = p_[2] ;
    p = p_[3] ;
    q = p_[4] ;
    r0 = p_[5] ;
    a_0 = p_[6];
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

double error(Vector a,double a0, double E_coh,   double &B, double &C11, double &C12, double &C44) {
    vector<Atom> Vect;

    calculate_heliocentric_lattice(Vect, "Ag", a.vec[6]);
    passParams(a.vec);
    //printParams(a.vec);



    double E_coh_f, B_f, C11_f, C12_f, C44_f;

    E_coh_f = calculate_energy(Vect, a.vec[6], MATR) / double(Vect.size());

    calculate_parameters(a.vec[6],E_coh_f, B_f, C11_f, C12_f, C44_f);

    double err =
            (E_coh - E_coh_f) * (E_coh - E_coh_f) +
            (B - B_f) * (B - B_f) +
            (C11 - C11_f) * (C11 - C11_f) +
            (C12 - C12_f) * (C12 - C12_f) +
            (C44 - C44_f) * (C44 - C44_f);

    return err;
}

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

    multimap<double, Vector> simplex;

    for(int i=0; i<8; i++) {
        Vector tmp(1);
        if(i != 7) {
            for(int i=0; i < 7; i++) {
                tmp.vec[i] = rand() % 10000 / 1000.0;
            }
        }
        simplex.insert({error(tmp, a_0, E_coh, B, c11, c12, c44), tmp});
    }

    double const reflection_coeff = 1;
    double const compression_coeff = .5;
    double const strectching_coeff = 2;
    double const global_strectching_coeff = .5;

    while(simplex.begin()->first > 10e-10) {
        Vector x_h = simplex.rbegin()->second;
        double f_h = simplex.rbegin()->first;
        Vector x_g = (++simplex.rbegin())->second;
        double f_g = (++simplex.rbegin())->first;
        Vector x_l = simplex.begin()->second;
        double f_l = simplex.begin()->first;

        cout << "BEST ERROR: " << f_l << endl;
        //cout << "WORST ERROR: " << f_h << endl;
        //printParams(x_h.vec);


        Vector x_c;

        for(auto it = ++simplex.rbegin(); it != simplex.rend(); ++it) {
            x_c += it->second;
        }
        x_c /= 7;
        Vector x_r = (x_c*(1+reflection_coeff) - x_h*reflection_coeff).strict();
        double f_r = error(x_r, a_0, E_coh,  B, c11, c12, c44);

        if(f_r < f_l) {
            Vector x_e = (x_c * (1 - strectching_coeff) + x_r * strectching_coeff).strict();
            double f_e = error(x_e, a_0, E_coh, B, c11, c12, c44);

            simplex.erase(--simplex.end());
            if(f_e < f_r) {
                simplex.insert({f_e, x_e});
            } else {
                simplex.insert({f_r, x_r});
            }
        } else if(f_r < f_g) {
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});
        } else if(f_r < f_h) {
            simplex.erase(--simplex.end());
            simplex.insert({f_r, x_r});

            Vector x_s = (x_c * (1 - compression_coeff) + x_h * compression_coeff).strict();
            double f_s = error(x_s, a_0,E_coh, B, c11, c12, c44);
            if(f_s < f_h) {
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for(auto& x_i : tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * global_strectching_coeff;
                    simplex.insert({error(new_vec, a_0, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        } else {

            Vector x_s = x_c * (1 - compression_coeff) + x_h * compression_coeff;
            double f_s = error(x_s, a_0, E_coh, B, c11, c12, c44);
            if(f_s < f_h) {
                simplex.erase(--simplex.end());
                simplex.insert({f_s, x_s});
            } else {
                multimap<double, Vector> tmp(simplex);
                simplex.clear();

                for(auto& x_i : tmp) {
                    Vector new_vec = x_l + (x_i.second - x_l) * global_strectching_coeff;
                    simplex.insert({error(new_vec, a_0, E_coh, B, c11, c12, c44), new_vec});
                }
            }
        }
    }

    printParams(simplex.begin()->second.vec);
    passParams(simplex.begin()->second.vec);

    double E_coh_f, B_f, C11_f, C12_f, C44_f;

    vector<Atom> Vect1;
    calculate_heliocentric_lattice(Vect1, "Ag", simplex.begin()->second.vec[6]);
    E_coh_f = calculate_energy(Vect1, simplex.begin()->second.vec[6], MATR) / double(Vect1.size());

    calculate_parameters(simplex.begin()->second.vec[6], E_coh_f, B_f, C11_f, C12_f, C44_f);


    cout << "E_coh_f: " << E_coh_f << endl;
    cout << "B_f: " << B_f << endl;
    cout << "C11_f: " << C11_f << endl;
    cout << "C12_f: " << C12_f << endl;
    cout << "C44_f: " << C44_f << endl;

    return 0;

}

