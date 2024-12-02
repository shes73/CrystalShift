#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"

Element periodic_table[118] = {
    {"H", 1}, {"He", 2},
    {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
    {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26},
    {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
    {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
    {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54},
    {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62},
    {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78},
    {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86},
    {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94},
    {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102},
    {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
    {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}
};

// // Function to get an element by index
// Element get_element(int index) {
//     if (index < 0 || index >= 118) {
//         Element invalid = {"", -1};
//         return invalid;
//     }
//     return periodic_table[index];
// }

// // Function to get atomic number by symbol
// int get_atomic_number(const char* symbol) {
//     for (int i = 0; i < 118; i++) {
//         if (strcmp(periodic_table[i].symbol, symbol) == 0) {
//             return periodic_table[i].atomic_number;
//         }
//     }
//     return -1; // Symbol not found
// }

double normalize_coordinate(double coord) {
    if (coord < 0.0) {
        coord += 50.0; 
        return fmod(coord, 1.0);
    } else return fmod(coord, 1.0);
}

// to check duplicates need to create normalized structure
Structure create_normalized_structure(const Structure* original) {
    Structure normalized;
    normalized.atoms = malloc(sizeof(Atom) * original->atom_count);
    normalized.atom_count = original->atom_count;

    for (int i = 0; i < original->atom_count; i++) {
        strcpy(normalized.atoms[i].element, original->atoms[i].element);
        normalized.atoms[i].x = normalize_coordinate(original->atoms[i].x);
        normalized.atoms[i].y = normalize_coordinate(original->atoms[i].y);
        normalized.atoms[i].z = normalize_coordinate(original->atoms[i].z);
    }

    return normalized;
}

int check_for_duplicates(const Structure* normalized) {
    for (int i = 0; i < normalized->atom_count; i++) {
        for (int j = i + 1; j < normalized->atom_count; j++) {
            if (strcmp(normalized->atoms[i].element, normalized->atoms[j].element) == 0 &&
                fabs(normalize_coordinate(normalized->atoms[i].x) - normalize_coordinate(normalized->atoms[j].x)) < 1e-6 &&
                fabs(normalize_coordinate(normalized->atoms[i].y) - normalize_coordinate(normalized->atoms[j].y)) < 1e-6 &&
                fabs(normalize_coordinate(normalized->atoms[i].z) - normalize_coordinate(normalized->atoms[j].z)) < 1e-6) {
                return 1; // there are duplicates
            }
        }
    }
    return 0; // there are no duplicates
}

void remove_duplicate_atoms(Structure* original, const Structure* normalized) {
    int new_atom_count = 0;
    Atom* unique_atoms = malloc(sizeof(Atom) * original->atom_count);

    for (int i = 0; i < normalized->atom_count; i++) {
        int is_duplicate = 0;
        for (int j = 0; j < new_atom_count; j++) {
            if (strcmp(normalized->atoms[i].element, unique_atoms[j].element) == 0 &&
                fabs(normalize_coordinate(normalized->atoms[i].x) - normalize_coordinate(unique_atoms[j].x)) < 1e-6 &&
                fabs(normalize_coordinate(normalized->atoms[i].y) - normalize_coordinate(unique_atoms[j].y)) < 1e-6 &&
                fabs(normalize_coordinate(normalized->atoms[i].z) - normalize_coordinate(unique_atoms[j].z)) < 1e-6) {
                is_duplicate = 1;
                break;
            }
        }
        if (!is_duplicate) {
            unique_atoms[new_atom_count++] = original->atoms[i];
        }
    }

    free(original->atoms);
    original->atoms = unique_atoms;
    original->atom_count = new_atom_count;
}
