#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define M_PI 3.14159265358979323846
#define MIN_PTS 2 // минимальное количество атомов в молекуле
#define DIM 3 // размерность пространства

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

typedef struct {
    double x, y, z;
} Vector3;

typedef struct {
    int h, k, l;
} MillerIndex;

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

// typical bond lengths (slightly increased for convinience)
const double bond_lengths[119][119] = {
    [1][1] = 0.85, // H-H bond length
    [1][6] = 1.35, // H-C bond length
    [1][7] = 1.50, // H-N bond length
    [1][8] = 1.11, // H-O bond length
    [1][9] = 1.10, // H-F bond length
    [1][17] = 1.30, // H-Cl bond length
    [1][35] = 1.50, // H-Br bond length
    [1][53] = 1.70, // H-I bond length
    [6][1] = 1.35, // C-H bond length
    [6][6] = 1.54, // C-C bond length
    [6][7] = 2.10, // C-N bond length
    [6][8] = 1.50, // C-O bond length
    [6][9] = 1.34, // C-F bond length
    [6][17] = 1.76, // C-Cl bond length
    [6][35] = 2.10, // C-Br bond length
    [6][16] = 2.55, // C-S bond length
    [6][34] = 2.71, // C-Se bond length
    [7][1] = 1.50, // N-H bond length
    [7][7] = 1.75, // N-N bond length
    [7][6] = 2.10, // N-C bond length
    [7][8] = 1.90, // N-O bond length
    [8][1] = 1.11, // O-H bond length
    [8][6] = 1.50, // O-C bond length
    [8][7] = 1.90, // O-N bond length
    [8][8] = 1.55, // O-O bond length
    [9][1] = 1.10, // F-H bond length
    [9][6] = 1.34, // F-C bond length
    [17][1] = 1.30, // Cl-H bond length
    [17][6] = 1.76, // Cl-C bond length
    [35][1] = 1.50, // Br-H bond length
    [35][6] = 2.10, // Br-C bond length
    [53][1] = 1.70, // I-H bond length
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

typedef struct KDNode {
    Atom atom;
    struct KDNode* left;
    struct KDNode* right;
} KDNode;

double distance(Atom a1, Atom a2) {
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double get_bond_length(Atom a1, Atom a2) {
    int atomic_number1 = get_atomic_number(a1.element);
    int atomic_number2 = get_atomic_number(a2.element);
    if (atomic_number1 == -1 || atomic_number2 == -1) {
        return 2.8; // default bond length if elements are not found
    }

    // check if the bond length is defined
    if (bond_lengths[atomic_number1][atomic_number2] > 0) {
        return bond_lengths[atomic_number1][atomic_number2];
    }
    if (bond_lengths[atomic_number2][atomic_number1] > 0) {
        return bond_lengths[atomic_number2][atomic_number1];
    }

    return 2.8;
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

void rangeSearch(KDNode* root, Atom target, double eps, int depth, Atom** neighbors, int* neighbor_count, double lattice[3][3]) {
    if (root == NULL) return;

    if (distance(root->atom, target) < eps) {
        (*neighbors)[(*neighbor_count)++] = root->atom;
    }

    int cd = depth % DIM;
    if ((cd == 0 && target.x - eps < root->atom.x) ||
        (cd == 1 && target.y - eps < root->atom.y) ||
        (cd == 2 && target.z - eps < root->atom.z)) {
        rangeSearch(root->left, target, eps, depth + 1, neighbors, neighbor_count, lattice);
    }
    if ((cd == 0 && target.x + eps >= root->atom.x) ||
        (cd == 1 && target.y + eps >= root->atom.y) ||
        (cd == 2 && target.z + eps >= root->atom.z)) {
        rangeSearch(root->right, target, eps, depth + 1, neighbors, neighbor_count, lattice);
    }
}

void expandCluster(Atom* atoms, int n, int* labels, int cluster_id, int index, KDNode* root, double lattice[3][3]) {
    labels[index] = cluster_id;

    Atom* neighbors = (Atom*)malloc(n * sizeof(Atom));
    if (neighbors == NULL) {
        fprintf(stderr, "Memory allocation failed for neighbors\n");
        return;
    }

    int neighbor_count = 0;
    rangeSearch(root, atoms[index], 3.5, 0, &neighbors, &neighbor_count, lattice);

    for (int i = 0; i < neighbor_count; i++) {
        int neighbor_index = -1;
        for (int j = 0; j < n; j++) {
            if (strcmp(atoms[j].element, neighbors[i].element) == 0 &&
                distance(atoms[j], neighbors[i]) < 0.01) { // Используем PBC-расстояние
                neighbor_index = j;
                break;
            }
        }
        if (neighbor_index == -1) continue;

        double eps = get_bond_length(atoms[index], atoms[neighbor_index]);

        if (distance(atoms[index], atoms[neighbor_index]) <= eps) {
            if (labels[neighbor_index] == -1) {
                labels[neighbor_index] = cluster_id;
            }
            if (labels[neighbor_index] == 0) {
                labels[neighbor_index] = cluster_id;
                expandCluster(atoms, n, labels, cluster_id, neighbor_index, root, lattice);
            }
        }
    }
    free(neighbors);
}

void dbscan(Atom* atoms, int n, int* labels, double lattice[3][3]) {
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
        rangeSearch(root, atoms[i], 3.5, 0, &neighbors, &neighbor_count, lattice);

        if (neighbor_count < MIN_PTS) {
            labels[i] = -1; // mark as noise
            free(neighbors);
            continue;
        }

        cluster_id++;
        expandCluster(atoms, n, labels, cluster_id, i, root, lattice);
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

void calculate_geometric_centers(Atom* atoms, int atom_count, int* cluster_labels, int cluster_count) {
    printf("Geometric centers of clusters:\n");

    double** centers = (double**)malloc(cluster_count * sizeof(double*));
    for (int i = 0; i < cluster_count; i++) {
        centers[i] = (double*)malloc(3 * sizeof(double));
        centers[i][0] = 0.0;
        centers[i][1] = 0.0;
        centers[i][2] = 0.0;
    }

    for (int cluster_id = 1; cluster_id <= cluster_count; cluster_id++) {
        int count = 0;

        for (int i = 0; i < atom_count; i++) {
            if (cluster_labels[i] == cluster_id) {
                centers[cluster_id - 1][0] += atoms[i].x;
                centers[cluster_id - 1][1] += atoms[i].y;
                centers[cluster_id - 1][2] += atoms[i].z;
                count++;
            }
        }

        if (count > 0) {
            centers[cluster_id - 1][0] /= count;
            centers[cluster_id - 1][1] /= count;
            centers[cluster_id - 1][2] /= count;
        }

        printf("Cluster %d: Center (%.3f, %.3f, %.3f)\n", cluster_id, centers[cluster_id - 1][0], centers[cluster_id - 1][1], centers[cluster_id - 1][2]);
    }

    for (int i = 0; i < cluster_count; i++) {
        free(centers[i]);
    }
    free(centers);
}

void calculate_cluster_covariance_matrix(Atom* atoms, int atom_count, int* cluster_labels, double lattice[3][3], double covariance_matrix[3][3]) {
    double mean[3] = {0.0, 0.0, 0.0};
    int total_atoms = 0;

    for (int i = 0; i < atom_count; i++) {
        if (cluster_labels[i] > 0) {
            double wrapped_x = fmod(atoms[i].x, lattice[0][0]);
            double wrapped_y = fmod(atoms[i].y, lattice[1][1]);
            double wrapped_z = fmod(atoms[i].z, lattice[2][2]);

            mean[0] += wrapped_x;
            mean[1] += wrapped_y;
            mean[2] += wrapped_z;
            total_atoms++;
        }
    }

    for (int i = 0; i < 3; i++) {
        mean[i] /= total_atoms;
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            covariance_matrix[i][j] = 0.0;
            for (int k = 0; k < atom_count; k++) {
                if (cluster_labels[k] > 0) {
                    double wrapped_x = fmod(atoms[k].x, lattice[0][0]);
                    double wrapped_y = fmod(atoms[k].y, lattice[1][1]);
                    double wrapped_z = fmod(atoms[k].z, lattice[2][2]);

                    double diff_i = (i == 0) ? wrapped_x - mean[0] : (i == 1) ? wrapped_y - mean[1] : wrapped_z - mean[2];
                    double diff_j = (j == 0) ? wrapped_x - mean[0] : (j == 1) ? wrapped_y - mean[1] : wrapped_z - mean[2];

                    covariance_matrix[i][j] += diff_i * diff_j;
                }
            }
            covariance_matrix[i][j] /= total_atoms;
        }
    }
}

void power_iteration_matrix(double matrix[3][3], Vector3 *eigenvector, int max_iter, double tol) {
    double v[3] = {1.0, 1.0, 1.0};
    double norm, new_v[3];

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < 3; i++) {
            new_v[i] = 0.0;
            for (int j = 0; j < 3; j++) {
                new_v[i] += matrix[i][j] * v[j];
            }
        }
        norm = sqrt(new_v[0] * new_v[0] + new_v[1] * new_v[1] + new_v[2] * new_v[2]);
        for (int i = 0; i < 3; i++) {
            new_v[i] /= norm;
        }
        double diff = fabs(new_v[0] - v[0]) + fabs(new_v[1] - v[1]) + fabs(new_v[2] - v[2]);
        if (diff < tol) break;
        memcpy(v, new_v, sizeof(v));
    }
    eigenvector->x = v[0];
    eigenvector->y = v[1];
    eigenvector->z = v[2];
}

void perform_pca(Atom* atoms, int atom_count, int* cluster_labels, Vector3* principal_axis, MillerIndex* hkl) {
    double covariance_matrix[3][3] = {0};
    calculate_cluster_covariance_matrix(atoms, atom_count, cluster_labels, structure.lattice, covariance_matrix);

    power_iteration_matrix(covariance_matrix, principal_axis, 100, 1e-6);

    printf("Principal Axis: (%.3f, %.3f, %.3f)\n", principal_axis->x, principal_axis->y, principal_axis->z);

    // calculate the slip plane (hkl)
    double norm = sqrt(principal_axis->x * principal_axis->x + principal_axis->y * principal_axis->y + principal_axis->z * principal_axis->z);
    hkl->h = round(principal_axis->x / norm * 10);
    hkl->k = round(principal_axis->y / norm * 10);
    hkl->l = round(principal_axis->z / norm * 10);

    printf("Slip Plane (hkl): (%d, %d, %d)\n", hkl->h, hkl->k, hkl->l);
}

void integrate_pca_with_structure(Structure *structure, int *cluster_labels) {
    Vector3 principal_axis;
    MillerIndex hkl;
    perform_pca(structure->atoms, structure->atom_count, cluster_labels, &principal_axis, &hkl);
}

bool validate_layer(int* cluster_labels, int num_clusters, Structure* structure, Vector3 principal_axis) {
    if (num_clusters < 3) {
        printf("Validation failed: Less than 3 clusters in the layer.\n");
        return false;
    }

    // calculate geometric centers of clusters
    Vector3* centers = (Vector3*)malloc(num_clusters * sizeof(Vector3));
    int* counts = (int*)calloc(num_clusters, sizeof(int));

    for (int i = 0; i < num_clusters; i++) {
        centers[i].x = centers[i].y = centers[i].z = 0.0;
    }

    for (int i = 0; i < structure->atom_count; i++) {
        int cluster_id = cluster_labels[i] - 1;
        if (cluster_id >= 0) {
            centers[cluster_id].x += structure->atoms[i].x;
            centers[cluster_id].y += structure->atoms[i].y;
            centers[cluster_id].z += structure->atoms[i].z;
            counts[cluster_id]++;
        }
    }

    for (int i = 0; i < num_clusters; i++) {
        if (counts[i] > 0) {
            centers[i].x /= counts[i];
            centers[i].y /= counts[i];
            centers[i].z /= counts[i];
        }
    }

    // calculate intra-layer and inter-layer distances
    double intra_layer_dist = 0.0, inter_layer_dist = 1e6;
    int intra_count = 0, inter_count = 0;

    for (int i = 0; i < num_clusters; i++) {
        for (int j = i + 1; j < num_clusters; j++) {
            double dx = centers[i].x - centers[j].x;
            double dy = centers[i].y - centers[j].y;
            double dz = centers[i].z - centers[j].z;

            // project the distance vector onto the principal axis
            double projected_dist = fabs(dx * principal_axis.x + dy * principal_axis.y + dz * principal_axis.z);

            if (projected_dist < 5.0) { // threshold for intra-layer distance
                intra_layer_dist += sqrt(dx * dx + dy * dy + dz * dz);
                intra_count++;
            } else {
                inter_layer_dist = fmin(inter_layer_dist, sqrt(dx * dx + dy * dy + dz * dz));
                inter_count++;
            }
        }
    }

    if (intra_count > 0 && inter_count > 0 && intra_layer_dist / intra_count >= inter_layer_dist) {
        printf("Validation failed: Intra-layer distances are not smaller than inter-layer distances.\n");
        free(centers);
        free(counts);
        return false;
    }

    // check R^2 of fitted centers in a line
    double mean_x = 0, mean_y = 0, mean_z = 0;
    for (int i = 0; i < num_clusters; i++) {
        mean_x += centers[i].x;
        mean_y += centers[i].y;
        mean_z += centers[i].z;
    }
    mean_x /= num_clusters;
    mean_y /= num_clusters;
    mean_z /= num_clusters;

    double ss_total = 0, ss_residual = 0;
    for (int i = 0; i < num_clusters; i++) {
        double proj = (centers[i].x - mean_x) * principal_axis.x +
                      (centers[i].y - mean_y) * principal_axis.y +
                      (centers[i].z - mean_z) * principal_axis.z;
        double predicted_x = mean_x + proj * principal_axis.x;
        double predicted_y = mean_y + proj * principal_axis.y;
        double predicted_z = mean_z + proj * principal_axis.z;

        double error_x = centers[i].x - predicted_x;
        double error_y = centers[i].y - predicted_y;
        double error_z = centers[i].z - predicted_z;
        ss_residual += error_x * error_x + error_y * error_y + error_z * error_z;

        double deviation_x = centers[i].x - mean_x;
        double deviation_y = centers[i].y - mean_y;
        double deviation_z = centers[i].z - mean_z;
        ss_total += deviation_x * deviation_x + deviation_y * deviation_y + deviation_z * deviation_z;
    }

    double r_squared = 1 - (ss_residual / ss_total);
    if (r_squared < 0.8) {
        printf("Validation failed: R^2 is below 0.8 (%.3f).\n", r_squared);
        free(centers);
        free(counts);
        return false;
    }

    free(centers);
    free(counts);
    return true;
}

bool try_alternative_layer_directions(int* cluster_labels, int num_clusters, Structure* structure) {
    double covariance_matrix[DIM][DIM];
    calculate_cluster_covariance_matrix(structure->atoms, structure->atom_count, cluster_labels, structure->lattice, covariance_matrix);

    Vector3 eigenvectors[DIM];
    for (int i = 0; i < DIM; i++) {
        power_iteration_matrix(covariance_matrix, &eigenvectors[i], 100, 1e-6);

        // orthogonalize the eigenvector to the previous ones
        for (int j = 0; j < i; j++) {
            double dot_product = eigenvectors[i].x * eigenvectors[j].x +
                                  eigenvectors[i].y * eigenvectors[j].y +
                                  eigenvectors[i].z * eigenvectors[j].z;
            eigenvectors[i].x -= dot_product * eigenvectors[j].x;
            eigenvectors[i].y -= dot_product * eigenvectors[j].y;
            eigenvectors[i].z -= dot_product * eigenvectors[j].z;
        }

        // normalize the eigenvector
        double norm = sqrt(eigenvectors[i].x * eigenvectors[i].x +
                            eigenvectors[i].y * eigenvectors[i].y +
                            eigenvectors[i].z * eigenvectors[i].z);
        eigenvectors[i].x /= norm;
        eigenvectors[i].y /= norm;
        eigenvectors[i].z /= norm;
    }

    for (int i = 0; i < DIM; i++) {
        Vector3 principal_axis = eigenvectors[i];
        printf("Attempt %d: Using eigenvector %d as the principal axis\n", i + 1, i + 1);

        if (validate_layer(cluster_labels, num_clusters, structure, principal_axis)) {
            return true;
        }
    }
    return false;
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

    if (freopen("output.txt", "w", stdout) == NULL) {
        fprintf(stderr, "Error redirecting stdout to file!\n");
        fclose(output_file);
        return 1;
    }

    if (parse_cif(filename, pizdec) != 0) {
        return !pizdec;
    }

    convert_fractional_to_direct();
    create_supercell(&structure, 5, 5, 5);

    int* labels = (int*)malloc(structure.atom_count * sizeof(int));
    if (labels == NULL) {
        fprintf(stderr, "Memory allocation failed for labels\n");
        return 1;
    }

    for (int i = 0; i < structure.atom_count; i++) {
        labels[i] = 0;
    }

    dbscan(structure.atoms, structure.atom_count, labels, structure.lattice);

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

    Vector3 principal_axis;
    MillerIndex hkl;
    perform_pca(structure.atoms, structure.atom_count, labels, &principal_axis, &hkl);

    if (!validate_layer(labels, num_clusters, &structure, principal_axis)) {
        printf("Layer validation failed. Trying alternative directions...\n");
        if (!try_alternative_layer_directions(labels, num_clusters, &structure)) {
            printf("All alternative layer directions failed. Exiting.\n");
            free(labels);
            free(structure.atoms);
            fclose(stdout);
            return 1;
        }
    }

    printf("Layer validation passed. Proceeding with analysis.\n");
    calculate_geometric_centers(structure.atoms, structure.atom_count, labels, num_clusters);
    integrate_pca_with_structure(&structure, labels);

    free(labels);
    free(structure.atoms);
    fclose(stdout);

    return 0;
}
