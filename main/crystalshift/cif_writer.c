#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "cif_writer.h"

#define M_PI 3.14159265358979323846

void remove_duplicate_atoms() {
    int new_atom_count = 0;
    for (int i = 0; i < structure.atom_count; i++) {
        int is_duplicate = 0;
        for (int j = 0; j < new_atom_count; j++) {
            if (strcmp(structure.atoms[i].element, structure.atoms[j].element) == 0 &&
                structure.atoms[i].x == structure.atoms[j].x &&
                structure.atoms[i].y == structure.atoms[j].y &&
                structure.atoms[i].z == structure.atoms[j].z) {
                is_duplicate = 1;
                break;
            }
        }
        if (!is_duplicate) {
            structure.atoms[new_atom_count] = structure.atoms[i];
            new_atom_count++;
        }
    }
    structure.atom_count = new_atom_count;
    structure.atoms = realloc(structure.atoms, sizeof(Atom) * new_atom_count);
}

void write_cif(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file for writing!\n");
        return;
    }

    fprintf(file, "%s_generated_by_CrystalShift\n", filename);
    fprintf(file, "_symmetry_cell_setting triclinic\n");
    fprintf(file, "_symmetry_space_group_name_H-M 'P 1'\n");

    fprintf(file, "loop_\n");
    fprintf(file, "_symmetry_equiv_pos_site_id\n");
    fprintf(file, "_symmetry_equiv_pos_as_xyz\n");
    fprintf(file, "1 x,y,z\n");
    fprintf(file, "_cell_length_a %lf\n", structure.lengths[0]);
    fprintf(file, "_cell_length_b %lf\n", structure.lengths[1]);
    fprintf(file, "_cell_length_c %lf\n", structure.lengths[2]);
    fprintf(file, "_cell_angle_alpha %lf\n", structure.angles[0]);
    fprintf(file, "_cell_angle_beta %lf\n", structure.angles[1]);
    fprintf(file, "_cell_angle_gamma %lf\n", structure.angles[2]);

    double lattice_volume = structure.lattice[0][0] * (structure.lattice[1][1] * structure.lattice[2][2] - structure.lattice[1][2] * structure.lattice[2][1]) -
                            structure.lattice[0][1] * (structure.lattice[1][0] * structure.lattice[2][2] - structure.lattice[1][2] * structure.lattice[2][0]) +
                            structure.lattice[0][2] * (structure.lattice[1][0] * structure.lattice[2][1] - structure.lattice[1][1] * structure.lattice[2][0]);

    fprintf(file, "_cell_volume %lf\n", lattice_volume);

    fprintf(file, "loop_\n");
    fprintf(file, "_atom_site_label\n");
    fprintf(file, "_atom_site_type_symbol\n");
    fprintf(file, "_atom_site_fract_x\n");
    fprintf(file, "_atom_site_fract_y\n");
    fprintf(file, "_atom_site_fract_z\n");

    remove_duplicate_atoms();

    for (int i = 0; i < structure.atom_count; i++) {
        fprintf(file, "%s%d %s %lf %lf %lf\n", structure.atoms[i].element, i+1, structure.atoms[i].element, structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z);
    }

    fprintf(file, "\n#END\n");

    fclose(file);
}