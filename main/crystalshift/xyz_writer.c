#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "xyz_writer.h"

#define M_PI 3.14159265358979323846

int write_xyz(const char *filename, const Structure *structure, int xyz_format) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file for writing!\n");
        return 1;
    }

    fprintf(file, "%d\n", structure->atom_count);

    if (xyz_format == 1) {
        fprintf(file, "Generated by CrystalShift\n"); // general case
    } else if (xyz_format == 2) {                     // extended case: note that there're also may be charges, velocities, etc. CHECK IT MANUALLY!!!
        fprintf(file, "Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\"\n",
                structure->lattice[0][0], structure->lattice[0][1], structure->lattice[0][2],
                structure->lattice[1][0], structure->lattice[1][1], structure->lattice[1][2],
                structure->lattice[2][0], structure->lattice[2][1], structure->lattice[2][2]);
    }

    Structure normalized_structure_xyz = create_normalized_structure(structure);

    if (check_for_duplicates(&normalized_structure_xyz)) {
        printf("Warning: duplicate coordinates found! Please, check it manually.\n");
    }

    free(normalized_structure_xyz.atoms);

    for (int i = 0; i < structure->atom_count; i++) {
        double x = structure->atoms[i].x * structure->lattice[0][0] +
                   structure->atoms[i].y * structure->lattice[1][0] +
                   structure->atoms[i].z * structure->lattice[2][0];
        double y = structure->atoms[i].x * structure->lattice[0][1] +
                   structure->atoms[i].y * structure->lattice[1][1] +
                   structure->atoms[i].z * structure->lattice[2][1];
        double z = structure->atoms[i].x * structure->lattice[0][2] +
                   structure->atoms[i].y * structure->lattice[1][2] +
                   structure->atoms[i].z * structure->lattice[2][2];

        fprintf(file, "%s %lf %lf %lf\n",
                structure->atoms[i].element,
                x,
                y,
                z);
    }

    fclose(file);
    return 0;
}
