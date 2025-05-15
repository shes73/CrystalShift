#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "cif_parser.h"

#define M_PI 3.14159265358979323846
#define MAX_OPERATIONS 36

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

void parse_symmetry_operations_from_file(FILE *file) {
    char line[256];
    int in_loop = 0;
    int symmetry_column_index = -1;
    int column_count = 0;
    int data_started = 0;
    int count = 0;
    char *operations[192]; // Убедитесь, что этот массив достаточно велик для всех операций

    rewind(file);

    const char *symmetry_keys[] = {
        "_symmetry_equiv_pos_as_xyz",
        "_space_group_symop_operation_xyz"
    };

    while (fgets(line, sizeof(line), file)) {
        // Start of loop
        if (strstr(line, "loop_")) {
            in_loop = 1;
            symmetry_column_index = -1;
            column_count = 0;
            data_started = 0;
            continue;
        }

        // Inside loop: header lines
        if (in_loop && line[0] == '_') {
            for (int i = 0; i < sizeof(symmetry_keys) / sizeof(symmetry_keys[0]); i++) {
                if (strstr(line, symmetry_keys[i])) {
                    symmetry_column_index = column_count;
                }
            }
            column_count++;
            continue;
        }

        // Start reading data lines
        if (in_loop && line[0] != '_' && symmetry_column_index != -1) {
            data_started = 1;
        }

        if (data_started) {
            if (line[0] == '\n' || line[0] == '\r' || line[0] == '#') break;

            char *start = line;
            while (*start == ' ' || *start == '\t') start++;

            // SHELX-style: quoted single string
            if (*start == '\'' || *start == '"') {
                start++;
                char *end = strchr(start, '\'');
                if (!end) end = strchr(start, '"');
                if (!end) end = strchr(start, '\n');
                if (end) *end = '\0';
                if (count < 192) {
                    operations[count++] = strdup(start);
                }
            } else {
                // Mercury-style: tokenized columns
                int col = 0;
                char *token = strtok(start, " \t\n\r");
                while (token) {
                    if (col == symmetry_column_index) {
                        if (count < 192) {
                            operations[count++] = strdup(token);
                        }
                        break;
                    }
                    token = strtok(NULL, " \t\n\r");
                    col++;
                }
            }
        }
    }

    if (count > 0) {
        structure.symmetry_operations = malloc(sizeof(char *) * count);
        for (int i = 0; i < count; i++) {
            structure.symmetry_operations[i] = operations[i];
        }
        structure.symmetry_operations_count = count;
    } else {
        structure.symmetry_operations = NULL;
        structure.symmetry_operations_count = 0;
    }
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
        free(structure.symmetry_operations);
        return 1;
    }

    int lines_with_coords[3];
    int line_with_atom_type;
    int parsing_atoms = 0;
    int *line_number = NULL; // pointer to an array
    int line_count = 0;
    int capacity = 0;
    int current_line = 0;
    int *matching_lines_atom_site = NULL; // pointer to the array with number of lines starting from "_atom_site"
    int count = 0; // counter for the number of elements in the array matching_lines

    char *last_atom_site_line = NULL; // Variable to store the last line starting with "_atom_site_"
    int last_atom_site_line_number = 0; // Variable to store the line number of the last "_atom_site_" line
    int first_atom_site_line_number = 0; // Variable to store the line number of the first "_atom_site_" line
    int in_atom_site_section = 0; // Flag to indicate if we are in the section with "_atom_site_" lines

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
                    free(structure.symmetry_operations);
                    return 1;
                }
            }
            line_number[line_count++] = current_line;
        }

        if (strncmp(line, "_atom_site_", 11) == 0) {
            in_atom_site_section = 1; // We are in the section with "_atom_site_" lines
            count++;
            matching_lines_atom_site = realloc(matching_lines_atom_site, count * sizeof(int));
            if (matching_lines_atom_site == NULL) {
                perror("Ошибка при выделении памяти");
                fclose(file);
                free(structure.atoms);
                free(structure.symmetry_operations);
                return EXIT_FAILURE;
            }

            matching_lines_atom_site[count - 1] = current_line; // Сохраняем номер текущей строки

            // Сохраняем последнюю строку, начинающуюся с "_atom_site_"
            if (last_atom_site_line != NULL) {
                free(last_atom_site_line);
            }
            last_atom_site_line = strdup(line);
            last_atom_site_line_number = current_line; // Сохраняем номер строки

            // Сохраняем номер первой строки, начинающейся с "_atom_site_"
            if (first_atom_site_line_number == 0) {
                first_atom_site_line_number = current_line;
            }
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

        // Check for empty line after "_atom_site_" lines
        if (in_atom_site_section && line[0] == '\n') {
            break; // Stop searching after the first empty line following "_atom_site_" lines
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
                Atom *temp_atoms = realloc(structure.atoms, sizeof(Atom) * structure.atom_capacity);
                if (temp_atoms == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    fclose(file);
                    free(line_number);
                    free(structure.symmetry_operations);
                    return 1;
                }
                structure.atoms = temp_atoms;
            }

            char *inputCopy = strdup(line);
            if (inputCopy == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                free(line_number);
                free(structure.atoms);
                free(structure.symmetry_operations);
                return 1;
            }

            const char *delimiter = " ";
            char *token;
            int count = 0;

            token = strtok(inputCopy, delimiter);
            while (token != NULL) {
                count++;
                if (count == (line_with_atom_type - first_atom_site_line_number + 1)) {
                    strncpy(structure.atoms[atom_index].element, token, sizeof(structure.atoms[atom_index].element) - 1);
                    structure.atoms[atom_index].element[sizeof(structure.atoms[atom_index].element) - 1] = '\0';
                }

                if (count == (lines_with_coords[0] - first_atom_site_line_number + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].x);
                }

                if (count == (lines_with_coords[1] - first_atom_site_line_number + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].y);
                }

                if (count == (lines_with_coords[2] - first_atom_site_line_number + 1)) {
                    sscanf(token, "%lf", &structure.atoms[atom_index].z);
                }

                token = strtok(NULL, delimiter);
            }

            free(inputCopy);
            atom_index++;
        }
    }

    structure.atom_count = atom_index;

    calculate_lattice_matrix();
    parse_symmetry_operations_from_file(file);

    fclose(file);
    free(line_number);
    free(matching_lines_atom_site);

    return 0;
}
