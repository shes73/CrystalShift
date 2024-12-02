#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef struct {
    char element[3];
    double x, y, z;
} Atom;

typedef struct {
    double lengths[3];
    double angles[3];
    double lattice[3][3];
    Atom *atoms;
    int atom_count;
    int atom_capacity;
} Structure;

Structure structure = {.atom_capacity = 100}; // Начальная емкость 100 атомов

typedef struct {
    char symbol[3];
    int atomic_number;
} Element;

Element periodic_table[] = {
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

int get_atomic_number(const char *symbol) {
    for (unsigned int i = 0; i < sizeof(periodic_table) / sizeof(periodic_table[0]); i++) {
        if (strncmp(symbol, periodic_table[i].symbol, 2) == 0) {
            return periodic_table[i].atomic_number;
        }
    }
    return -1;
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

    for (int i = 0; i < structure.atom_count; i++) {
        fprintf(file, "%s%d %s %lf %lf %lf\n", structure.atoms[i].element, i+1, structure.atoms[i].element, structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z);
    }

    fprintf(file, "\n#END\n");

    fclose(file);
}

int main() {
    char filename[256];
    char xyz_filename[256];
    char new_cif_filename[256];
    int xyz_format;

    printf("Enter the path to the CIF file: ");
    scanf("%s", filename);

    if (parse_cif(filename) != 0) {
        return 1;
    }

    printf("Choose the format for writing xyz file:\n");
    printf("1. Standard xyz format;\n");
    printf("2. Extended xyz format (with lattice parameters in the comment line);\n");
    printf("Enter your choice: ");
    scanf("%d", &xyz_format);
    getchar();

    if (xyz_format <= 0 || xyz_format > 2) {
        printf("Error! There's no such option!\n");
        return 1;
    }

    printf("Enter the path to the xyz file: ");
    scanf("%s", xyz_filename);
    write_xyz(xyz_filename, &structure, xyz_format);

    free(structure.atoms);
    structure.atoms = NULL;
    structure.atom_count = 0;
    structure.atom_capacity = 100;

    if (parse_xyz(xyz_filename, xyz_format) != 0) {
        return 1;
    }

    printf("Enter the path to the new CIF file: ");
    scanf("%s", new_cif_filename);
    write_cif(new_cif_filename);

    free(structure.atoms);

    return 0;
}
