#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "basis_changer.h"

#define M_PI 3.14159265358979323846

// BASIS CHANGER FUNCTIONS: B = C^-1 * A * C, where A - old basis, B - new basis, C - transformation matrix

void get_the_transformation_matrix(double transformation_matrix[3][3]) {
    printf("Enter the vectors of the basis transformation matrix.\n");
    for (int i = 0; i < 3; i++) {
        printf("%c vector:\n", 'a' + i);
        for (int j = 0; j < 3; j++) {
            if (scanf("%lf", &transformation_matrix[i][j]) != 1) {
                printf("Error! Please, enter float numbers.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
}

void swap_rows(double matrix[3][3], double inverse[3][3], int row1, int row2) {
    for (int j = 0; j < 3; j++) {
        double temp = matrix[row1][j];
        matrix[row1][j] = matrix[row2][j];
        matrix[row2][j] = temp;
        
        temp = inverse[row1][j];
        inverse[row1][j] = inverse[row2][j];
        inverse[row2][j] = temp;
    }
}

void inverse_transformation_matrix(double transformation_matrix[3][3], double inversed_transformation_matrix[3][3]) {
    double temporary_matrix[3][3];
    double ratio;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inversed_transformation_matrix[i][j] = (i == j) ? 1.0 : 0.0;
            temporary_matrix[i][j] = transformation_matrix[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        if (temporary_matrix[i][i] == 0.0) {
            int found = 0;
            for (int j = i + 1; j < 3; j++) {
                if (temporary_matrix[j][i] != 0.0) {
                    swap_rows(temporary_matrix, inversed_transformation_matrix, i, j);
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("Matrix is singular and cannot be inverted.\n");
                exit(EXIT_FAILURE);
            }
        }

        for (int j = 0; j < 3; j++) {
            if (i != j) {
                ratio = temporary_matrix[j][i] / temporary_matrix[i][i];
                for (int k = 0; k < 3; k++) {
                    temporary_matrix[j][k] -= ratio * temporary_matrix[i][k];
                    inversed_transformation_matrix[j][k] -= ratio * inversed_transformation_matrix[i][k];
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        ratio = temporary_matrix[i][i];
        for (int j = 0; j < 3; j++) {
            temporary_matrix[i][j] /= ratio;
            inversed_transformation_matrix[i][j] /= ratio;
        }
    }
}

void calculate_new_basis_matrix(double transformation_matrix[3][3], double new_basis_matrix[3][3], double new_cell_lengths[3], double new_cell_angles[3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            new_basis_matrix[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                new_basis_matrix[i][j] += transformation_matrix[i][k] * structure.lattice[k][j];
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        new_cell_lengths[i] = sqrt(new_basis_matrix[i][0] * new_basis_matrix[i][0] + new_basis_matrix[i][1] * new_basis_matrix[i][1]
                                  + new_basis_matrix[i][2] * new_basis_matrix[i][2]);
    }

    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        new_cell_angles[i] = (acos((new_basis_matrix[j][0] * new_basis_matrix[k][0] +
                                    new_basis_matrix[j][1] * new_basis_matrix[k][1] +
                                    new_basis_matrix[j][2] * new_basis_matrix[k][2]) /
                                   (new_cell_lengths[j] * new_cell_lengths[k]))) * 180.0 / M_PI;
    }

    for (int i = 0; i < 3; i++) {
        structure.lengths[i] = new_cell_lengths[i];
        structure.angles[i] = new_cell_angles[i];
        for (int j = 0; j < 3; j++) {
            structure.lattice[i][j] = new_basis_matrix[i][j];
        }
    }
}

void calculate_new_coords(double inversed_transformation_matrix[3][3]) {
    for (int i = 0; i < structure.atom_count; i++) {
        double new_x = inversed_transformation_matrix[0][0] * structure.atoms[i].x + inversed_transformation_matrix[1][0] * structure.atoms[i].y + inversed_transformation_matrix[2][0] * structure.atoms[i].z;
        double new_y = inversed_transformation_matrix[0][1] * structure.atoms[i].x + inversed_transformation_matrix[1][1] * structure.atoms[i].y + inversed_transformation_matrix[2][1] * structure.atoms[i].z;
        double new_z = inversed_transformation_matrix[0][2] * structure.atoms[i].x + inversed_transformation_matrix[1][2] * structure.atoms[i].y + inversed_transformation_matrix[2][2] * structure.atoms[i].z;

        structure.atoms[i].x = new_x;
        structure.atoms[i].y = new_y;
        structure.atoms[i].z = new_z;
    }
}

void create_supercell_if_needed(double transformation_matrix[3][3]) {
    int supercell_factor[3] = {1, 1, 1};

    if (transformation_matrix[0][0] != 0 && transformation_matrix[0][1] == 0 && transformation_matrix[0][2] == 0) {
        supercell_factor[0] = (int)transformation_matrix[0][0];
    }
    if (transformation_matrix[1][1] != 0 && transformation_matrix[1][0] == 0 && transformation_matrix[1][2] == 0) {
        supercell_factor[1] = (int)transformation_matrix[1][1];
    }
    if (transformation_matrix[2][2] != 0 && transformation_matrix[2][0] == 0 && transformation_matrix[2][1] == 0) {
        supercell_factor[2] = (int)transformation_matrix[2][2];
    }

    int coef_x = supercell_factor[0];
    int coef_y = supercell_factor[1];
    int coef_z = supercell_factor[2];

    if (supercell_factor[0] > 1 || supercell_factor[1] > 1 || supercell_factor[2] > 1) {
        char response[3];
        printf("Supercell detected. Do you want to add new atom coordinates in the supercell? (yes/no): ");
        scanf("%3s", response);
        if (strcmp(response, "yes") == 0) {
            printf("Supercell factors: [%d, %d, %d]\n", coef_x, coef_y, coef_z);
            printf("Original atom count: %d\n", structure.atom_count);
            int new_atom_count = structure.atom_count * coef_x * coef_y * coef_z;

            Atom *old_atoms = malloc(sizeof(Atom) * structure.atom_count);
            memcpy(old_atoms, structure.atoms, sizeof(Atom) * structure.atom_count);

            Atom *new_atoms = malloc(sizeof(Atom) * new_atom_count);
            int index = 0;

            for (int dx = 0; dx < coef_x; dx++) {
                for (int dy = 0; dy < coef_y; dy++) {
                    for (int dz = 0; dz < coef_z; dz++) {
                        for (int i = 0; i < structure.atom_count; i++) {
                            new_atoms[index].x = old_atoms[i].x + dx;
                            new_atoms[index].y = old_atoms[i].y + dy;
                            new_atoms[index].z = old_atoms[i].z + dz;
                            strcpy(new_atoms[index].element, old_atoms[i].element);
                            index++;
                        }
                    }
                }
            }

            free(structure.atoms);

            // update initial structure with old and new atoms
            structure.atoms = new_atoms;
            structure.atom_count = new_atom_count;
            free(old_atoms);

            printf("Updated atom count: %d\n", structure.atom_count);
        }
    }
}
