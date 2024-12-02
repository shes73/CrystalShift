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

Structure structure1 = {.atom_capacity = 100};
Structure structure2 = {.atom_capacity = 100};

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
    for (int i = 0; i < sizeof(periodic_table) / sizeof(periodic_table[0]); i++) {
        if (strncmp(symbol, periodic_table[i].symbol, 2) == 0) {
            return periodic_table[i].atomic_number;
        }
    }
    return -1;
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

void calculate_lattice_matrix(Structure *structure) {
    double a = structure->lengths[0];
    double b = structure->lengths[1];
    double c = structure->lengths[2];
    double alpha = structure->angles[0] * M_PI / 180.0;
    double beta = structure->angles[1] * M_PI / 180.0;
    double gamma = structure->angles[2] * M_PI / 180.0;

    structure->lattice[0][0] = a;
    structure->lattice[0][1] = 0.0;
    structure->lattice[0][2] = 0.0;

    structure->lattice[1][0] = b * cos(gamma);
    structure->lattice[1][1] = b * sin(gamma);
    structure->lattice[1][2] = 0.0;

    structure->lattice[2][0] = c * cos(beta);
    structure->lattice[2][1] = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
    structure->lattice[2][2] = c * sqrt(1.0 - cos(beta) * cos(beta) - pow((cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma), 2));
}

int parse_cif(const char *filename, Structure *structure) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }

    char line[256];
    int atom_index = 0;
    structure->atoms = malloc(sizeof(Atom) * structure->atom_capacity);
    if (structure->atoms == NULL) {
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

        if (sscanf(line, "_cell_length_a %lf", &structure->lengths[0]) == 1) continue;
        if (sscanf(line, "_cell_length_b %lf", &structure->lengths[1]) == 1) continue;
        if (sscanf(line, "_cell_length_c %lf", &structure->lengths[2]) == 1) continue;
        if (sscanf(line, "_cell_angle_alpha %lf", &structure->angles[0]) == 1) continue;
        if (sscanf(line, "_cell_angle_beta %lf", &structure->angles[1]) == 1) continue;
        if (sscanf(line, "_cell_angle_gamma %lf", &structure->angles[2]) == 1) continue;

        if (strstr(line, "_atom_site") != NULL) {
            if (line_count >= capacity) {
                capacity = (capacity == 0) ? 1 : capacity * 2;
                line_number = (int *)realloc(line_number, capacity * sizeof(int));
                if (line_number == NULL) {
                    perror("Error allocating memory!");
                    fclose(file);
                    free(structure->atoms);
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

            if (atom_index >= structure->atom_capacity) {
                structure->atom_capacity *= 2;
                structure->atoms = realloc(structure->atoms, sizeof(Atom) * structure->atom_capacity);
                if (structure->atoms == NULL) {
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
                free(structure->atoms);
                return 1;
            }

            const char *delimiter = " ";
            char *token;
            int count = 0;

            token = strtok(inputCopy, delimiter);
            while (token != NULL) {
                count++;
                if (count == (line_with_atom_type - line_number[0] + 1)) {
                    strncpy(structure->atoms[atom_index].element, token, sizeof(structure->atoms[atom_index].element) - 1);
                    structure->atoms[atom_index].element[sizeof(structure->atoms[atom_index].element) - 1] = '\0';
                }

                if (count == (lines_with_coords[0] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure->atoms[atom_index].x);
                }

                if (count == (lines_with_coords[1] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure->atoms[atom_index].y);
                }

                if (count == (lines_with_coords[2] - line_number[0] + 1)) {
                    sscanf(token, "%lf", &structure->atoms[atom_index].z);
                }

                token = strtok(NULL, delimiter);
            }

            free(inputCopy);
            atom_index++;
        }
    }

    structure->atom_count = atom_index;
    fclose(file);
    free(line_number);

    calculate_lattice_matrix(structure);

    return 0;
}

int compare_atoms_by_coordinates(const void *a, const void *b) {
    Atom *atom1 = (Atom *)a;
    Atom *atom2 = (Atom *)b;

    if (atom1->x != atom2->x) {
        return (atom1->x > atom2->x) - (atom1->x < atom2->x);
    }
    if (atom1->y != atom2->y) {
        return (atom1->y > atom2->y) - (atom1->y < atom2->y);
    }
    return (atom1->z > atom2->z) - (atom1->z < atom2->z);
}

void sort_atoms_by_coordinates(Structure *structure) {
    qsort(structure->atoms, structure->atom_count, sizeof(Atom), compare_atoms_by_coordinates);
}


void calculate_statistics(Structure *structure1, Structure *structure2) {
    if (structure1->atom_count != structure2->atom_count) {
        fprintf(stderr, "Error: The number of atoms in the two structures is not the same.\n");
        return;
    }

    double sum_sq_x = 0.0, sum_sq_y = 0.0, sum_sq_z = 0.0;

    for (int i = 0; i < structure1->atom_count; i++) {
        double diff_x = structure1->atoms[i].x - structure2->atoms[i].x;
        double diff_y = structure1->atoms[i].y - structure2->atoms[i].y;
        double diff_z = structure1->atoms[i].z - structure2->atoms[i].z;

        sum_sq_x += diff_x * diff_x;
        sum_sq_y += diff_y * diff_y;
        sum_sq_z += diff_z * diff_z;
    }

    double rmsd_x = sqrt(sum_sq_x / structure1->atom_count);
    double rmsd_y = sqrt(sum_sq_y / structure1->atom_count);
    double rmsd_z = sqrt(sum_sq_z / structure1->atom_count);

    printf("RMSD (x y z): %lf %lf %lf\n", rmsd_x, rmsd_y, rmsd_z);
}

int main() {
    char filename1[256];
    char filename2[256];

    printf("Enter the path to the first CIF file: \n");
    scanf("%s", filename1);

    if (parse_cif(filename1, &structure1) != 0) {
        return 1;
    }

    sort_atoms_by_coordinates(&structure1);

    printf("Enter the path to the second CIF file: \n");
    scanf("%s", filename2);

    if (parse_cif(filename2, &structure2) != 0) {
        return 1;
    }

    sort_atoms_by_coordinates(&structure2);

    calculate_statistics(&structure1, &structure2);

    free(structure1.atoms);
    free(structure2.atoms);

    return 0;
}
