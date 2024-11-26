#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "cif_parser.h"

#define M_PI 3.14159265358979323846

char *strdup(const char *s) {
    size_t len = strlen(s) + 1;
    char *d = malloc(len);
    if (d == NULL) return NULL;
    return (char *)memcpy(d, s, len);
}

void calculate_lattice_matrix() {
    double a = structure.lengths[0];
    double b = structure.lengths[1];
    double c = structure.lengths[2];
    double alpha = structure.angles[0] * M_PI / 180.0;
    double beta = structure.angles[1] * M_PI / 180.0;
    double gamma = structure.angles[2] * M_PI / 180.0;

    structure.lattice[0][0] = a;
    structure.lattice[0][1] = 0.0;
    structure.lattice[0][2] = 0.0;

    structure.lattice[1][0] = b * cos(gamma);
    structure.lattice[1][1] = b * sin(gamma);
    structure.lattice[1][2] = 0.0;

    structure.lattice[2][0] = c * cos(beta);
    structure.lattice[2][1] = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
    structure.lattice[2][2] = c * sqrt(1.0 - cos(beta) * cos(beta) - pow((cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma), 2));
}

void remove_errors(char *line) {
    char *src = line, *dst = line;
    while (*src) {
        if (*src == '(') { // removing experimental errors from cif to avoid errors with parsing coords
            while (*src && *src != ')') src++;
            if (*src) src++;
        } else {
            *dst++ = *src++;
        }
    }
    *dst = '\0';
}

int parse_cif(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }

    char line[256];
    int atom_index = 0;
    structure.atoms = malloc(sizeof(Atom) * structure.atom_capacity);
    if (structure.atoms == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        return 1;
    }

    int lines_with_coords[3];
    int line_with_atom_type;
    int parsing_atoms = 0;
    int *line_number = NULL; // pointer to an array
    int line_count = 0;
    int capacity = 0;
    int current_line = 0;

    while (fgets(line, sizeof(line), file)) {
        current_line++;

        if (sscanf(line, "_cell_length_a %lf", &structure.lengths[0]) == 1) continue;
        if (sscanf(line, "_cell_length_b %lf", &structure.lengths[1]) == 1) continue;
        if (sscanf(line, "_cell_length_c %lf", &structure.lengths[2]) == 1) continue;
        if (sscanf(line, "_cell_angle_alpha %lf", &structure.angles[0]) == 1) continue;
        if (sscanf(line, "_cell_angle_beta %lf", &structure.angles[1]) == 1) continue;
        if (sscanf(line, "_cell_angle_gamma %lf", &structure.angles[2]) == 1) continue;

        if (strstr(line, "_atom_site") != NULL) {
            if (line_count >= capacity) {
                capacity = (capacity == 0) ? 1 : capacity * 2;
                line_number = (int *)realloc(line_number, capacity * sizeof(int));
                if (line_number == NULL) {
                    perror("Error allocating memory!");
                    fclose(file);
                    free(structure.atoms);
                    return 1;
                }
            }
            line_number[line_count++] = current_line; // array with numbers of lines = amount of elements in atomic coords line, after (last element + 1) we may parse coords
        }

        if (strstr(line, "_atom_site_type_symbol") != NULL) {
            line_with_atom_type = current_line;
        }

        if (strstr(line, "_atom_site_fract_x") != NULL) {
            lines_with_coords[0] = current_line;
        }

        if (strstr(line, "_atom_site_fract_y") != NULL) {
            lines_with_coords[1] = current_line;
        }

        if (strstr(line, "_atom_site_fract_z") != NULL) {
            lines_with_coords[2] = current_line;
        }
    }

    current_line = 0;

    rewind(file);

    while (fgets(line, sizeof(line), file)) {
        current_line++;
        if (line_count > 0 && current_line == line_number[line_count - 1]) {
            parsing_atoms = 1;
            break;
        }
    }

    while (fgets(line, sizeof(line), file)) {
        if (parsing_atoms) {
            if (line[0] == '\n' || line[0] == '#') {
                break;
            }

            remove_errors(line);

            char element[3];
            if (sscanf(line, "%2s", element) != 1 || strlen(element) == 0) {
                continue;
            }

            if (atom_index >= structure.atom_capacity) {
                structure.atom_capacity *= 2;
                structure.atoms = realloc(structure.atoms, sizeof(Atom) * structure.atom_capacity);
                if (structure.atoms == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    fclose(file);
                    free(line_number);
                    return 1;
                }
            }

            char *inputCopy = strdup(line);
            if (inputCopy == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                free(line_number);
                free(structure.atoms);
                return 1;
            }

            const char *delimiter = " ";
            char *token;
            int count = 0;

            token = strtok(inputCopy, delimiter);
            while (token != NULL) {
                count++;
                if (count == (line_with_atom_type - line_number[0] + 1)) {
                    strncpy(structure.atoms[atom_index].element, token, sizeof(structure.atoms[atom_index].element) - 1);
                    structure.atoms[atom_index].element[sizeof(structure.atoms[atom_index].element) - 1] = '\0';
                }

                if (count == (lines_with_coords[0] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].x);
                }

                if (count == (lines_with_coords[1] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].y);
                }

                if (count == (lines_with_coords[2] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].z);
                }

                token = strtok(NULL, delimiter);
            }

            free(inputCopy);
            atom_index++;
        }
    }

    structure.atom_count = atom_index;
    fclose(file);
    free(line_number);

    calculate_lattice_matrix();

    return 0;
}
