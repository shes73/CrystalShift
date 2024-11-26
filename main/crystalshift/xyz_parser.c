#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "xyz_parser.h"

#define M_PI 3.14159265358979323846

Structure structure;

void inverse_matrix(double transformation_matrix[3][3], double inversed_transformation_matrix[3][3]) {
    double temporary_matrix[3][3];
    double ratio;

    // initialization of all matrices before Gauss-Jordan inversion
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inversed_transformation_matrix[i][j] = (i == j) ? 1.0 : 0.0;
            temporary_matrix[i][j] = transformation_matrix[i][j];
        }
    }

    // reducing a matrix to upper triangular form
    for (int i = 0; i < 3; i++) {
        if (temporary_matrix[i][i] == 0.0) {
            printf("Matrix is singular and cannot be inverted.\n");
            exit(EXIT_FAILURE);
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

    // reducing a matrix to E
    for (int i = 0; i < 3; i++) {
        ratio = temporary_matrix[i][i];
        for (int j = 0; j < 3; j++) {
            temporary_matrix[i][j] /= ratio;
            inversed_transformation_matrix[i][j] /= ratio;
        }
    }
}

int parse_xyz(const char *filename, int xyz_format) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }

    char line[256];
    int number_of_atoms_xyz;

    if (fgets(line, sizeof(line), file) == NULL || sscanf(line, "%d", &number_of_atoms_xyz) != 1) {
        fprintf(stderr, "Error parsing number of atoms!\n");
        fclose(file);
        return 1;
    }

    structure.atoms = malloc(sizeof(Atom) * structure.atom_capacity);
    if (structure.atoms == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        return 1;
    }

    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error reading file!\n");
        free(structure.atoms);
        fclose(file);
        return 1;
    }

    if (xyz_format == 2 && sscanf(line, "Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\"",
                   &structure.lattice[0][0], &structure.lattice[0][1], &structure.lattice[0][2],
                   &structure.lattice[1][0], &structure.lattice[1][1], &structure.lattice[1][2],
                   &structure.lattice[2][0], &structure.lattice[2][1], &structure.lattice[2][2]) != 9) {
        fprintf(stderr, "Error parsing lattice parameters!\n");
        free(structure.atoms);
        fclose(file);
        return 1;
    }

    double inverse_lattice[3][3];
    inverse_matrix(structure.lattice, inverse_lattice);

    int current_atom_index = 0;

    for (int i = 0; i < number_of_atoms_xyz; i++) {
        if (current_atom_index >= structure.atom_capacity) {
            structure.atom_capacity *= 2;
            Atom *tmp = realloc(structure.atoms, sizeof(Atom) * structure.atom_capacity);
            if (tmp == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                free(structure.atoms);
                fclose(file);
                return 1;
            }
            structure.atoms = tmp;
        }

        if (fgets(line, sizeof(line), file) == NULL) {
            fprintf(stderr, "Error reading atom data!\n");
            free(structure.atoms);
            fclose(file);
            return 1;
        }

        if (sscanf(line, "%s %lf %lf %lf", structure.atoms[current_atom_index].element,
                   &structure.atoms[current_atom_index].x, &structure.atoms[current_atom_index].y,
                   &structure.atoms[current_atom_index].z) != 4) {
            fprintf(stderr, "Error parsing atom data!\n");
            free(structure.atoms);
            fclose(file);
            return 1;
        }

        double x = structure.atoms[current_atom_index].x;
        double y = structure.atoms[current_atom_index].y;
        double z = structure.atoms[current_atom_index].z;

        double fx = x * inverse_lattice[0][0] + y * inverse_lattice[1][0] + z * inverse_lattice[2][0];
        double fy = x * inverse_lattice[0][1] + y * inverse_lattice[1][1] + z * inverse_lattice[2][1];
        double fz = x * inverse_lattice[0][2] + y * inverse_lattice[1][2] + z * inverse_lattice[2][2];

        structure.atoms[current_atom_index].x = fx;
        structure.atoms[current_atom_index].y = fy;
        structure.atoms[current_atom_index].z = fz;

        current_atom_index++;
    }

    structure.atom_count = current_atom_index;
    fclose(file);
    return 0;
}
