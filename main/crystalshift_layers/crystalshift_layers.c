#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define M_PI 3.14159265358979323846
#define MIN_PTS 2 // min number of atoms in a molecule
#define DIM 3 // dimensionality of the space

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
} Structure;

typedef struct {
    char symbol[3];
    int atomic_number;
} Element;

Structure structure;

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

// Typical bond lengths (slightly enlarged for common atomic combinations)
const double bond_lengths[118][118] = {
    [0][0] = 0.85, // H-H bond length
    [0][5] = 1.35, // H-C bond length
    [0][6] = 1.50, // H-N bond length
    [0][7] = 1.11, // H-O bond length
    [0][8] = 1.10, // H-F bond length
    [0][16] = 1.30, // H-Cl bond length
    [0][34] = 1.50, // H-Br bond length
    [0][52] = 1.70, // H-I bond length
    [5][0] = 1.35, // C-H bond length
    [5][5] = 1.54, // C-C bond length
    [5][6] = 2.10, // C-N bond length
    [5][7] = 1.50, // C-O bond length
    [5][8] = 1.34, // C-F bond length
    [5][16] = 1.76, // C-Cl bond length
    [5][34] = 1.93, // C-Br bond length
    [5][15] = 2.55, // C-S bond length
    [5][33] = 2.71, // C-Se bond length
    [6][0] = 1.50, // N-H bond length
    [6][6] = 1.75, // N-N bond length
    [6][5] = 2.10, // N-C bond length
    [6][7] = 1.90, // N-O bond length
    [7][0] = 1.11, // O-H bond length
    [7][5] = 1.50, // O-C bond length
    [7][6] = 1.90, // O-N bond length
    [7][7] = 1.55, // O-O bond length
    [8][0] = 1.10, // F-H bond length
    [8][5] = 1.34, // F-C bond length
    [16][0] = 1.30, // Cl-H bond length
    [16][5] = 1.76, // Cl-C bond length
    [34][0] = 1.50, // Br-H bond length
    [34][5] = 1.93, // Br-C bond length
    [52][0] = 1.70, // I-H bond length
};

int get_atomic_number(const char *symbol) {
    for (unsigned int i = 0; i < sizeof(periodic_table) / sizeof(periodic_table[0]); i++) {
        if (strcmp(symbol, periodic_table[i].symbol) == 0) {
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

void remove_errors(char *line) {
    char *src = line, *dst = line;
    while (*src) {
        if (*src == '(') {
            while (*src && *src != ')') src++;
            if (*src) src++;
        } else {
            *dst++ = *src++;
        }
    }
    *dst = '\0';
}

int parse_cif(const char *filename, int pizdec) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        return pizdec;
    }

    char line[256];
    int atom_index = 0;
    structure.atoms = NULL;

    int parsing_atoms = 0;

    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "_cell_length_a %lf", &structure.lengths[0]) == 1) continue;
        if (sscanf(line, "_cell_length_b %lf", &structure.lengths[1]) == 1) continue;
        if (sscanf(line, "_cell_length_c %lf", &structure.lengths[2]) == 1) continue;
        if (sscanf(line, "_cell_angle_alpha %lf", &structure.angles[0]) == 1) continue;
        if (sscanf(line, "_cell_angle_beta %lf", &structure.angles[1]) == 1) continue;
        if (sscanf(line, "_cell_angle_gamma %lf", &structure.angles[2]) == 1) continue;

        if (strstr(line, "_atom_site_fract_z") != NULL) {
            parsing_atoms = 1;
            continue;
        }

        if (parsing_atoms) {
            if (line[0] == '\n' || line[0] == '#') {
                break;
            }

            remove_errors(line);

            char element[3];
            if (sscanf(line, "%2s", element) != 1 || strlen(element) == 0) {
                continue;
            }

            structure.atoms = realloc(structure.atoms, sizeof(Atom) * (atom_index + 1));
            if (structure.atoms == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                return pizdec;
            }

            sscanf(line, "%*s %s %lf %lf %lf",
                   structure.atoms[atom_index].element,
                   &structure.atoms[atom_index].x,
                   &structure.atoms[atom_index].y,
                   &structure.atoms[atom_index].z);
            atom_index++;
        }
    }

    structure.atom_count = atom_index;
    fclose(file);

    calculate_lattice_matrix();

    return !pizdec;
}

void convert_fractional_to_direct() {
    for (int i = 0; i < structure.atom_count; i++) {
        double x = structure.atoms[i].x;
        double y = structure.atoms[i].y;
        double z = structure.atoms[i].z;

        structure.atoms[i].x = x * structure.lattice[0][0] + y * structure.lattice[1][0] + z * structure.lattice[2][0];
        structure.atoms[i].y = x * structure.lattice[0][1] + y * structure.lattice[1][1] + z * structure.lattice[2][1];
        structure.atoms[i].z = x * structure.lattice[0][2] + y * structure.lattice[1][2] + z * structure.lattice[2][2];
    }
}

void calculate_covariance_matrix(double **data, int rows, int cols, double **cov) {
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < cols; j++) {
            cov[i][j] = 0.0;
            for (int k = 0; k < rows; k++) {
                cov[i][j] += data[k][i] * data[k][j];
            }
            cov[i][j] /= rows;
        }
    }
}

void transpose_matrix(double **matrix, int rows, int cols, double **transposed) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transposed[j][i] = matrix[i][j];
        }
    }
}

void power_iteration(double **matrix, int size, double *eigenvalue, double *eigenvector, int max_iterations, double tolerance) {
    double *b_k = (double *)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        b_k[i] = 1.0;
    }

    double *b_k1 = (double *)malloc(size * sizeof(double));
    double norm;

    for (int it = 0; it < max_iterations; it++) {
        for (int i = 0; i < size; i++) {
            b_k1[i] = 0;
            for (int j = 0; j < size; j++) {
                b_k1[i] += matrix[i][j] * b_k[j];
            }
        }

        norm = 0.0;
        for (int i = 0; i < size; i++) {
            norm += b_k1[i] * b_k1[i];
        }
        norm = sqrt(norm);

        for (int i = 0; i < size; i++) {
            b_k[i] = b_k1[i] / norm;
        }

        *eigenvalue = 0;
        for (int i = 0; i < size; i++) {
            double temp = 0;
            for (int j = 0; j < size; j++) {
                temp += matrix[i][j] * b_k[j];
            }
            *eigenvalue += b_k[i] * temp;
        }

        if (norm < tolerance) {
            break;
        }
    }

    memcpy(eigenvector, b_k, size * sizeof(double));
    free(b_k);
    free(b_k1);
}

typedef struct KDNode {
    Atom atom;
    struct KDNode* left;
    struct KDNode* right;
} KDNode;

double distance(Atom a1, Atom a2) {
    return sqrt(pow(a1.x - a2.x, 2) + pow(a1.y - a2.y, 2) + pow(a1.z - a2.z, 2));
}

double get_bond_length(Atom a1, Atom a2) {
    int atomic_number1 = get_atomic_number(a1.element);
    int atomic_number2 = get_atomic_number(a2.element);
    if (atomic_number1 == -1 || atomic_number2 == -1) {
        return 3.0; // default bond length if elements are not found
    }
    return bond_lengths[atomic_number1][atomic_number2];
}

KDNode* createNode(Atom atom) {
    KDNode* node = (KDNode*)malloc(sizeof(KDNode));
    if (node == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    node->atom = atom;
    node->left = node->right = NULL;
    return node;
}

KDNode* insert(KDNode* root, Atom atom, int depth) {
    if (root == NULL) return createNode(atom);

    int cd = depth % DIM;
    if ((cd == 0 && atom.x < root->atom.x) ||
        (cd == 1 && atom.y < root->atom.y) ||
        (cd == 2 && atom.z < root->atom.z)) {
        root->left = insert(root->left, atom, depth + 1);
    } else {
        root->right = insert(root->right, atom, depth + 1);
    }

    return root;
}

void rangeSearch(KDNode* root, Atom target, double eps, int depth, Atom** neighbors, int* neighbor_count) {
    if (root == NULL) return;

    if (distance(root->atom, target) < eps) {
        (*neighbors)[(*neighbor_count)++] = root->atom;
    }

    int cd = depth % DIM;
    if ((cd == 0 && target.x - eps < root->atom.x) ||
        (cd == 1 && target.y - eps < root->atom.y) ||
        (cd == 2 && target.z - eps < root->atom.z)) {
        rangeSearch(root->left, target, eps, depth + 1, neighbors, neighbor_count);
    }
    if ((cd == 0 && target.x + eps >= root->atom.x) ||
        (cd == 1 && target.y + eps >= root->atom.y) ||
        (cd == 2 && target.z + eps >= root->atom.z)) {
        rangeSearch(root->right, target, eps, depth + 1, neighbors, neighbor_count);
    }
}

void expandCluster(Atom* atoms, int n, int* labels, int cluster_id, int index, KDNode* root) {
    labels[index] = cluster_id;

    Atom* neighbors = (Atom*)malloc(n * sizeof(Atom));
    if (neighbors == NULL) {
        fprintf(stderr, "Memory allocation failed for neighbors\n");
        return;
    }

    int neighbor_count = 0;
    double eps = get_bond_length(atoms[index], atoms[index]);
    rangeSearch(root, atoms[index], eps, 0, &neighbors, &neighbor_count);

    for (int i = 0; i < neighbor_count; i++) {
        int neighbor_index = -1;
        for (int j = 0; j < n; j++) {
            if (strcmp(atoms[j].element, neighbors[i].element) == 0 &&
                atoms[j].x == neighbors[i].x &&
                atoms[j].y == neighbors[i].y &&
                atoms[j].z == neighbors[i].z) {
                neighbor_index = j;
                break;
            }
        }
        if (neighbor_index == -1) continue;

        if (labels[neighbor_index] == -1) {
            labels[neighbor_index] = cluster_id;
        }
        if (labels[neighbor_index] == 0) {
            labels[neighbor_index] = cluster_id;
            Atom* new_neighbors = (Atom*)malloc(n * sizeof(Atom));
            if (new_neighbors == NULL) {
                fprintf(stderr, "Memory allocation failed for new_neighbors\n");
                free(neighbors);
                return;
            }
            int new_neighbor_count = 0;
            eps = get_bond_length(atoms[neighbor_index], atoms[neighbor_index]);
            rangeSearch(root, atoms[neighbor_index], eps, 0, &new_neighbors, &new_neighbor_count);
            if (new_neighbor_count >= MIN_PTS) {
                expandCluster(atoms, n, labels, cluster_id, neighbor_index, root);
            }
            free(new_neighbors);
        }
    }
    free(neighbors);
}

void dbscan(Atom* atoms, int n, int* labels) {
    KDNode* root = NULL;
    for (int i = 0; i < n; i++) {
        root = insert(root, atoms[i], 0);
    }

    int cluster_id = 0;
    for (int i = 0; i < n; i++) {
        if (labels[i] != 0) continue;

        Atom* neighbors = (Atom*)malloc(n * sizeof(Atom));
        if (neighbors == NULL) {
            fprintf(stderr, "Memory allocation failed for neighbors\n");
            continue;
        }

        int neighbor_count = 0;
        double eps = get_bond_length(atoms[i], atoms[i]);
        rangeSearch(root, atoms[i], eps, 0, &neighbors, &neighbor_count);

        if (neighbor_count < MIN_PTS) {
            labels[i] = -1; // mark as noise -- probably the disorder
            free(neighbors);
            continue;
        }

        cluster_id++;
        expandCluster(atoms, n, labels, cluster_id, i, root);
        free(neighbors);
    }
}

void create_supercell(Structure* structure, int nx, int ny, int nz) {
    int new_atom_count = structure->atom_count * nx * ny * nz;
    Atom* new_atoms = (Atom*)malloc(new_atom_count * sizeof(Atom));
    if (new_atoms == NULL) {
        fprintf(stderr, "Memory allocation failed for new_atoms\n");
        return;
    }
    int new_atom_index = 0;

    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
                for (int i = 0; i < structure->atom_count; i++) {
                    new_atoms[new_atom_index].x = structure->atoms[i].x + ix * structure->lattice[0][0] + iy * structure->lattice[1][0] + iz * structure->lattice[2][0];
                    new_atoms[new_atom_index].y = structure->atoms[i].y + ix * structure->lattice[0][1] + iy * structure->lattice[1][1] + iz * structure->lattice[2][1];
                    new_atoms[new_atom_index].z = structure->atoms[i].z + ix * structure->lattice[0][2] + iy * structure->lattice[1][2] + iz * structure->lattice[2][2];
                    strcpy(new_atoms[new_atom_index].element, structure->atoms[i].element);
                    new_atom_index++;
                }
            }
        }
    }

    free(structure->atoms);
    structure->atoms = new_atoms;
    structure->atom_count = new_atom_count;
}

void calculate_molecule_boundaries(Atom* atoms, int n, int* labels, int cluster_id, double* min_coords, double* max_coords) {
    min_coords[0] = min_coords[1] = min_coords[2] = INFINITY;
    max_coords[0] = max_coords[1] = max_coords[2] = -INFINITY;

    for (int i = 0; i < n; i++) {
        if (labels[i] == cluster_id) {
            if (atoms[i].x < min_coords[0]) min_coords[0] = atoms[i].x;
            if (atoms[i].y < min_coords[1]) min_coords[1] = atoms[i].y;
            if (atoms[i].z < min_coords[2]) min_coords[2] = atoms[i].z;
            if (atoms[i].x > max_coords[0]) max_coords[0] = atoms[i].x;
            if (atoms[i].y > max_coords[1]) max_coords[1] = atoms[i].y;
            if (atoms[i].z > max_coords[2]) max_coords[2] = atoms[i].z;
        }
    }
}

void calculate_center_of_mass(Atom* atoms, int n, int* labels, int cluster_id, double* com) {
    com[0] = com[1] = com[2] = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (labels[i] == cluster_id) {
            com[0] += atoms[i].x;
            com[1] += atoms[i].y;
            com[2] += atoms[i].z;
            count++;
        }
    }

    com[0] /= count;
    com[1] /= count;
    com[2] /= count;
}

double* projections; // Global variable to hold projections

int compare_projections(const void* a, const void* b) {
    int cluster_a = *(int*)a;
    int cluster_b = *(int*)b;
    double projection_a = projections[cluster_a - 1];
    double projection_b = projections[cluster_b - 1];
    if (projection_a < projection_b) {
        return -1;
    } else if (projection_a > projection_b) {
        return 1;
    } else {
        return 0;
    }
}

void calculate_reciprocal_lattice(double lattice[3][3], double reciprocal_lattice[3][3]) {
    double volume = lattice[0][0] * (lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1]) -
                   lattice[0][1] * (lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0]) +
                   lattice[0][2] * (lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0]);

    reciprocal_lattice[0][0] = (lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1]) / volume;
    reciprocal_lattice[0][1] = (lattice[0][2] * lattice[2][1] - lattice[0][1] * lattice[2][2]) / volume;
    reciprocal_lattice[0][2] = (lattice[0][1] * lattice[1][2] - lattice[0][2] * lattice[1][1]) / volume;

    reciprocal_lattice[1][0] = (lattice[1][2] * lattice[2][0] - lattice[1][0] * lattice[2][2]) / volume;
    reciprocal_lattice[1][1] = (lattice[0][0] * lattice[2][2] - lattice[0][2] * lattice[2][0]) / volume;
    reciprocal_lattice[1][2] = (lattice[0][2] * lattice[1][0] - lattice[0][0] * lattice[1][2]) / volume;

    reciprocal_lattice[2][0] = (lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0]) / volume;
    reciprocal_lattice[2][1] = (lattice[0][1] * lattice[2][0] - lattice[0][0] * lattice[2][1]) / volume;
    reciprocal_lattice[2][2] = (lattice[0][0] * lattice[1][1] - lattice[0][1] * lattice[1][0]) / volume;
}

void find_closest_planes(double direction[3], double reciprocal_lattice[3][3], int miller_indices[3]) {
    double min_distance = INFINITY;
    int closest_indices[3] = {0, 0, 0};

    for (int h = -1; h <= 1; h++) {
        for (int k = -1; k <= 1; k++) {
            for (int l = -1; l <= 1; l++) {
                if (h == 0 && k == 0 && l == 0) continue;

                double plane_normal[3] = {h * reciprocal_lattice[0][0] + k * reciprocal_lattice[1][0] + l * reciprocal_lattice[2][0],
                                        h * reciprocal_lattice[0][1] + k * reciprocal_lattice[1][1] + l * reciprocal_lattice[2][1],
                                        h * reciprocal_lattice[0][2] + k * reciprocal_lattice[1][2] + l * reciprocal_lattice[2][2]};

                double distance = fabs(plane_normal[0] * direction[0] + plane_normal[1] * direction[1] + plane_normal[2] * direction[2]);

                if (distance < min_distance) {
                    min_distance = distance;
                    closest_indices[0] = h;
                    closest_indices[1] = k;
                    closest_indices[2] = l;
                }
            }
        }
    }

    miller_indices[0] = closest_indices[0];
    miller_indices[1] = closest_indices[1];
    miller_indices[2] = closest_indices[2];
}

// Function to calculate distance between two atoms
double atom_distance(Atom a1, Atom a2) {
    return sqrt(pow(a1.x - a2.x, 2) + pow(a1.y - a2.y, 2) + pow(a1.z - a2.z, 2));
}

void find_nearest_molecules(Atom* atoms, int n, int* labels, int num_clusters, double** distances) {
    for (int i = 0; i < num_clusters; i++) {
        for (int j = i + 1; j < num_clusters; j++) {
            double min_distance = INFINITY;
            for (int k = 0; k < n; k++) {
                if (labels[k] == i) {
                    for (int l = 0; l < n; l++) {
                        if (labels[l] == j) {
                            double dist = atom_distance(atoms[k], atoms[l]);
                            if (dist < min_distance) {
                                min_distance = dist;
                            }
                        }
                    }
                }
            }
            distances[i][j] = min_distance;
            distances[j][i] = min_distance;
        }
    }
}

void analyze_molecular_layers(Atom* atoms, int n, int* labels, int num_clusters) {
    if (num_clusters <= 0) {
        fprintf(stderr, "Invalid number of clusters: %d\n", num_clusters);
        return;
    }

    double* coms = (double*)malloc(num_clusters * 3 * sizeof(double));
    if (coms == NULL) {
        fprintf(stderr, "Memory allocation failed for coms\n");
        return;
    }

    for (int i = 0; i < num_clusters; i++) {
        calculate_center_of_mass(atoms, n, labels, i + 1, &coms[i * 3]);
    }

    double** data = (double**)malloc(3 * sizeof(double*));
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed for data\n");
        free(coms);
        return;
    }
    for (int i = 0; i < 3; i++) {
        data[i] = (double*)malloc(num_clusters * sizeof(double));
        if (data[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for data[%d]\n", i);
            free(coms);
            for (int j = 0; j < i; j++) {
                free(data[j]);
            }
            free(data);
            return;
        }
    }

    for (int i = 0; i < num_clusters; i++) {
        data[0][i] = coms[i * 3];
        data[1][i] = coms[i * 3 + 1];
        data[2][i] = coms[i * 3 + 2];
    }

    double** data_t = (double**)malloc(num_clusters * sizeof(double*));
    if (data_t == NULL) {
        fprintf(stderr, "Memory allocation failed for data_t\n");
        free(coms);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        return;
    }
    for (int i = 0; i < num_clusters; i++) {
        data_t[i] = (double*)malloc(3 * sizeof(double));
        if (data_t[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for data_t[%d]\n", i);
            free(coms);
            for (int j = 0; j < 3; j++) {
                free(data[j]);
            }
            free(data);
            for (int j = 0; j < i; j++) {
                free(data_t[j]);
            }
            free(data_t);
            return;
        }
    }

    transpose_matrix(data, 3, num_clusters, data_t);

    double** cov = (double**)malloc(3 * sizeof(double*));
    if (cov == NULL) {
        fprintf(stderr, "Memory allocation failed for cov\n");
        free(coms);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        return;
    }
    for (int i = 0; i < 3; i++) {
        cov[i] = (double*)malloc(3 * sizeof(double));
        if (cov[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for cov[%d]\n", i);
            free(coms);
            for (int j = 0; j < 3; j++) {
                free(data[j]);
            }
            free(data);
            for (int j = 0; j < num_clusters; j++) {
                free(data_t[j]);
            }
            free(data_t);
            for (int j = 0; j < i; j++) {
                free(cov[j]);
            }
            free(cov);
            return;
        }
    }

    calculate_covariance_matrix(data_t, num_clusters, 3, cov);

    double eigenvalues[3];
    double **eigenvectors = (double**)malloc(3 * sizeof(double*));
    if (eigenvectors == NULL) {
        fprintf(stderr, "Memory allocation failed for eigenvectors\n");
        free(coms);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        for (int i = 0; i < 3; i++) {
            free(cov[i]);
        }
        free(cov);
        return;
    }
    for (int i = 0; i < 3; i++) {
        eigenvectors[i] = (double*)malloc(3 * sizeof(double));
        if (eigenvectors[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for eigenvectors[%d]\n", i);
            free(coms);
            for (int j = 0; j < 3; j++) {
                free(data[j]);
            }
            free(data);
            for (int j = 0; j < num_clusters; j++) {
                free(data_t[j]);
            }
            free(data_t);
            for (int j = 0; j < 3; j++) {
                free(cov[j]);
            }
            free(cov);
            for (int j = 0; j < i; j++) {
                free(eigenvectors[j]);
            }
            free(eigenvectors);
            return;
        }
    }

    power_iteration(cov, 3, eigenvalues, eigenvectors[0], 1000, 1e-6);
    power_iteration(cov, 3, eigenvalues + 1, eigenvectors[1], 1000, 1e-6);
    power_iteration(cov, 3, eigenvalues + 2, eigenvectors[2], 1000, 1e-6);

    // Calculate the reciprocal lattice vectors
    double reciprocal_lattice[3][3];
    calculate_reciprocal_lattice(structure.lattice, reciprocal_lattice);

    // Find the closest planes to the direction vectors
    int miller_indices[3][3];
    find_closest_planes(eigenvectors[0], reciprocal_lattice, miller_indices[0]);
    find_closest_planes(eigenvectors[1], reciprocal_lattice, miller_indices[1]);
    find_closest_planes(eigenvectors[2], reciprocal_lattice, miller_indices[2]);

    for (int i = 0; i < 3; i++) {
        fprintf(stdout, "Miller indices of the layer plane %d: (%d, %d, %d)\n", i + 1, miller_indices[i][0], miller_indices[i][1], miller_indices[i][2]);
    }

    // Project the centers of mass onto the direction vectors
    projections = (double*)malloc(num_clusters * sizeof(double));
    if (projections == NULL) {
        fprintf(stderr, "Memory allocation failed for projections\n");
        free(coms);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        for (int i = 0; i < 3; i++) {
            free(cov[i]);
        }
        free(cov);
        for (int i = 0; i < 3; i++) {
            free(eigenvectors[i]);
        }
        free(eigenvectors);
        return;
    }
    for (int i = 0; i < num_clusters; i++) {
        projections[i] = coms[i * 3] * eigenvectors[0][0] + coms[i * 3 + 1] * eigenvectors[0][1] + coms[i * 3 + 2] * eigenvectors[0][2];
    }

    // Sort the clusters based on their projections
    int* sorted_clusters = (int*)malloc(num_clusters * sizeof(int));
    if (sorted_clusters == NULL) {
        fprintf(stderr, "Memory allocation failed for sorted_clusters\n");
        free(coms);
        free(projections);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        for (int i = 0; i < 3; i++) {
            free(cov[i]);
        }
        free(cov);
        for (int i = 0; i < 3; i++) {
            free(eigenvectors[i]);
        }
        free(eigenvectors);
        return;
    }
    for (int i = 0; i < num_clusters; i++) {
        sorted_clusters[i] = i + 1;
    }
    qsort(sorted_clusters, num_clusters, sizeof(int), compare_projections);

    // Identify the molecular layers
    int* layer_indices = (int*)malloc(num_clusters * sizeof(int));
    if (layer_indices == NULL) {
        fprintf(stderr, "Memory allocation failed for layer_indices\n");
        free(coms);
        free(sorted_clusters);
        free(projections);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        for (int i = 0; i < 3; i++) {
            free(cov[i]);
        }
        free(cov);
        for (int i = 0; i < 3; i++) {
            free(eigenvectors[i]);
        }
        free(eigenvectors);
        return;
    }
    int num_layers = 0;
    double layer_threshold = 3.0;
    layer_indices[0] = 0;
    for (int i = 1; i < num_clusters; i++) {
        double dz = fabs(projections[i] - projections[i - 1]);
        if (dz > layer_threshold) {
            num_layers++;
        }
        layer_indices[i] = num_layers;
    }
    num_layers++;

    // Print the molecular layers and the distance between neighboring layers
    fprintf(stdout, "Molecular layers:\n");
    for (int i = 0; i < num_layers; i++) {
        fprintf(stdout, "Layer %d: ", i + 1);
        for (int j = 0; j < num_clusters; j++) {
            if (layer_indices[j] == i) {
                fprintf(stdout, "%d ", sorted_clusters[j]);
            }
        }
        fprintf(stdout, "\n");
    }

    // Check if the structure is layered
    bool is_layered = true;
    for (int i = 1; i < num_layers; i++) {
        int j = 0;
        while (j < num_clusters && layer_indices[j] != i) {
            j++;
        }
        if (j >= num_clusters) {
            fprintf(stderr, "Error: j index out of bounds\n");
            free(coms);
            free(sorted_clusters);
            free(projections);
            free(layer_indices);
            for (int k = 0; k < 3; k++) {
                free(data[k]);
            }
            free(data);
            for (int k = 0; k < num_clusters; k++) {
                free(data_t[k]);
            }
            free(data_t);
            for (int k = 0; k < 3; k++) {
                free(cov[k]);
            }
            free(cov);
            for (int k = 0; k < 3; k++) {
                free(eigenvectors[k]);
            }
            free(eigenvectors);
            return;
        }
        int k = j + 1;
        while (k < num_clusters) {
            k++;
        }
        if (k > num_clusters) {
            fprintf(stderr, "Error: k index out of bounds\n");
            fprintf(stderr, "k index = %d\n", k);
            fprintf(stderr, "layer_indices[%d] = %d\n", k, i - 1);
            free(coms);
            free(sorted_clusters);
            free(projections);
            free(layer_indices);
            for (int l = 0; l < 3; l++) {
                free(data[l]);
            }
            free(data);
            for (int l = 0; l < num_clusters; l++) {
                free(data_t[l]);
            }
            free(data_t);
            for (int l = 0; l < 3; l++) {
                free(cov[l]);
            }
            free(cov);
            for (int l = 0; l < 3; l++) {
                free(eigenvectors[l]);
            }
            free(eigenvectors);
            return;
        }
        double dz = fabs(projections[j] - projections[k]);
        if (dz < layer_threshold) {
            is_layered = false;
            break;
        }
    }

    if (is_layered) {
        fprintf(stdout, "The structure is layered.\n");
    } else {
        fprintf(stdout, "The structure is not layered.\n");
    }

    // Find nearest molecules and calculate distances
    double** distances = (double**)malloc(num_clusters * sizeof(double*));
    if (distances == NULL) {
        fprintf(stderr, "Memory allocation failed for distances\n");
        free(coms);
        free(sorted_clusters);
        free(projections);
        free(layer_indices);
        for (int i = 0; i < 3; i++) {
            free(data[i]);
        }
        free(data);
        for (int i = 0; i < num_clusters; i++) {
            free(data_t[i]);
        }
        free(data_t);
        for (int i = 0; i < 3; i++) {
            free(cov[i]);
        }
        free(cov);
        for (int i = 0; i < 3; i++) {
            free(eigenvectors[i]);
        }
        free(eigenvectors);
        return;
    }
    for (int i = 0; i < num_clusters; i++) {
        distances[i] = (double*)malloc(num_clusters * sizeof(double));
        if (distances[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for distances[%d]\n", i);
            free(coms);
            free(sorted_clusters);
            free(projections);
            free(layer_indices);
            for (int j = 0; j < 3; j++) {
                free(data[j]);
            }
            free(data);
            for (int j = 0; j < num_clusters; j++) {
                free(data_t[j]);
            }
            free(data_t);
            for (int j = 0; j < 3; j++) {
                free(cov[j]);
            }
            free(cov);
            for (int j = 0; j < 3; j++) {
                free(eigenvectors[j]);
            }
            free(eigenvectors);
            for (int j = 0; j < i; j++) {
                free(distances[j]);
            }
            free(distances);
            return;
        }
    }

    find_nearest_molecules(atoms, n, labels, num_clusters, distances);

    // Print the distances between neighboring layers
    fprintf(stdout, "Distances between neighboring layers:\n");
    for (int i = 0; i < num_layers - 1; i++) {
        double min_distance = INFINITY;
        for (int j = 0; j < num_clusters; j++) {
            if (layer_indices[j] == i) {
                for (int k = 0; k < num_clusters; k++) {
                    if (layer_indices[k] == i + 1) {
                        if (distances[j][k] < min_distance) {
                            min_distance = distances[j][k];
                        }
                    }
                }
            }
        }
        fprintf(stdout, "Distance between layer %d and layer %d: %f\n", i + 1, i + 2, min_distance);
    }

    // Check if the structure is layered based on the distances
    double factor = 1.15;
    bool is_layered_by_distance = true;
    for (int i = 0; i < num_layers - 1; i++) {
        double min_distance = INFINITY;
        for (int j = 0; j < num_clusters; j++) {
            if (layer_indices[j] == i) {
                for (int k = 0; k < num_clusters; k++) {
                    if (layer_indices[k] == i + 1) {
                        if (distances[j][k] < min_distance) {
                            min_distance = distances[j][k];
                        }
                    }
                }
            }
        }
        if (min_distance > factor * min_distance) {
            is_layered_by_distance = false;
            break;
        }
    }

    if (is_layered_by_distance) {
        fprintf(stdout, "The structure is layered based on distances.\n");
    } else {
        fprintf(stdout, "The structure is not layered based on distances.\n");
    }

    free(coms);
    free(sorted_clusters);
    free(projections);
    free(layer_indices);
    for (int i = 0; i < 3; i++) {
        free(data[i]);
    }
    free(data);
    for (int i = 0; i < num_clusters; i++) {
        free(data_t[i]);
    }
    free(data_t);
    for (int i = 0; i < 3; i++) {
        free(cov[i]);
    }
    free(cov);
    for (int i = 0; i < 3; i++) {
        free(eigenvectors[i]);
    }
    free(eigenvectors);
    for (int i = 0; i < num_clusters; i++) {
        free(distances[i]);
    }
    free(distances);
}

int main() {
    char filename[256];
    int pizdec = 1;

    printf("Enter the path to the CIF file: ");
    scanf("%s", filename);

    FILE *output_file = fopen("output.txt", "w");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening output file!\n");
        return 1;
    }

    // Redirect stdout to the output file
    if (freopen("output.txt", "w", stdout) == NULL) {
        fprintf(stderr, "Error redirecting stdout to file!\n");
        fclose(output_file);
        return 1;
    }

    if (parse_cif(filename, pizdec) != 0) {
        return !pizdec;
    }

    convert_fractional_to_direct(); // for clear bond lengths checking

    create_supercell(&structure, 3, 3, 3);

    int* labels = (int*)malloc(structure.atom_count * sizeof(int));
    if (labels == NULL) {
        fprintf(stderr, "Memory allocation failed for labels\n");
        return 1;
    }

    for (int i = 0; i < structure.atom_count; i++) {
        labels[i] = 0;
    }

    dbscan(structure.atoms, structure.atom_count, labels);

    for (int i = 0; i < structure.atom_count; i++) {
        fprintf(stdout, "Atom %d (%s): Cluster %d\n", i, structure.atoms[i].element, labels[i]);
    }

    int num_clusters = 0;
    for (int i = 0; i < structure.atom_count; i++) {
        if (labels[i] > num_clusters) {
            num_clusters = labels[i];
        }
    }

    fprintf(stdout, "Number of clusters: %d\n", num_clusters);

    if (num_clusters > structure.atom_count) {
        fprintf(stderr, "Invalid number of clusters: %d\n", num_clusters);
        free(labels);
        free(structure.atoms);
        fclose(stdout);
        return 1;
    }

    analyze_molecular_layers(structure.atoms, structure.atom_count, labels, num_clusters);

    free(labels);
    free(structure.atoms);
    fclose(stdout);

    return 0;
}
