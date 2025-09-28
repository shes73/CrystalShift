#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "poscar_parser.h"

#define M_PI 3.14159265358979323846

void calculate_lengths_and_angles(Structure *structure) {
    double (*lattice)[3] = structure->lattice;

    // calculate lattices vectors lenths
    structure->lengths[0] = sqrt(lattice[0][0] * lattice[0][0] + lattice[0][1] * lattice[0][1] + lattice[0][2] * lattice[0][2]);
    structure->lengths[1] = sqrt(lattice[1][0] * lattice[1][0] + lattice[1][1] * lattice[1][1] + lattice[1][2] * lattice[1][2]);
    structure->lengths[2] = sqrt(lattice[2][0] * lattice[2][0] + lattice[2][1] * lattice[2][1] + lattice[2][2] * lattice[2][2]);

    double a = structure->lengths[0];
    double b = structure->lengths[1];
    double c = structure->lengths[2];

    // lattice angles
    structure->angles[0] = acos((lattice[1][0] * lattice[2][0] + lattice[1][1] * lattice[2][1] + lattice[1][2] * lattice[2][2]) / (b * c)) * 180.0 / M_PI;
    structure->angles[1] = acos((lattice[0][0] * lattice[2][0] + lattice[0][1] * lattice[2][1] + lattice[0][2] * lattice[2][2]) / (a * c)) * 180.0 / M_PI;
    structure->angles[2] = acos((lattice[0][0] * lattice[1][0] + lattice[0][1] * lattice[1][1] + lattice[0][2] * lattice[1][2]) / (a * b)) * 180.0 / M_PI;
}

int parse_poscar(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }

    char line[256];
    structure.atoms = malloc(sizeof(Atom) * structure.atom_capacity);
    if (structure.atoms == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        return 1;
    }

    // Skip the first line (comment)
    fgets(line, sizeof(line), file);

    // Read the scaling factor
    double scaling_factor;
    fgets(line, sizeof(line), file);
    sscanf(line, "%lf", &scaling_factor);

    // Read the lattice vectors
    for (int i = 0; i < 3; i++) {
        fgets(line, sizeof(line), file);
        sscanf(line, "%lf %lf %lf", &structure.lattice[i][0], &structure.lattice[i][1], &structure.lattice[i][2]);
        structure.lattice[i][0] *= scaling_factor;
        structure.lattice[i][1] *= scaling_factor;
        structure.lattice[i][2] *= scaling_factor;
    }

    // Read the element symbols
    fgets(line, sizeof(line), file);
    char elements[256][3];
    int unique_count = 0;
    char *token = strtok(line, " \n");
    while (token != NULL) {
        strncpy(elements[unique_count], token, sizeof(elements[unique_count]) - 1);
        elements[unique_count][sizeof(elements[unique_count]) - 1] = '\0';
        unique_count++;
        token = strtok(NULL, " \n");
    }

    // Read the element counts
    fgets(line, sizeof(line), file);
    int counts[256];
    token = strtok(line, " \n");
    int total_atoms = 0;
    for (int i = 0; i < unique_count; i++) {
        counts[i] = atoi(token);
        total_atoms += counts[i];
        token = strtok(NULL, " \n");
    }

    // Skip the "Direct" or "Cartesian" line
    fgets(line, sizeof(line), file);

    // Read the atom coordinates
    int current_atom_index = 0;
    for (int i = 0; i < unique_count; i++) {
        for (int j = 0; j < counts[i]; j++) {
            if (current_atom_index >= structure.atom_capacity) {
                structure.atom_capacity *= 2;
                structure.atoms = realloc(structure.atoms, sizeof(Atom) * structure.atom_capacity);
                if (structure.atoms == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    fclose(file);
                    return 1;
                }
            }

            fgets(line, sizeof(line), file);
            sscanf(line, "%lf %lf %lf", &structure.atoms[current_atom_index].x, &structure.atoms[current_atom_index].y, &structure.atoms[current_atom_index].z);
            strcpy(structure.atoms[current_atom_index].element, elements[i]);
            current_atom_index++;
        }
    }

    structure.atom_count = current_atom_index;
    calculate_lengths_and_angles(&structure);
    fclose(file);

    return 0;

}
