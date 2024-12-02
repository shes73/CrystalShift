#ifndef STRUCTURES_H
#define STRUCTURES_H

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


typedef struct {
    char symbol[3];
    int atomic_number;
} Element;

extern Structure structure;
extern Element periodic_table[118];

double normalize_coordinate(double coord);
Structure create_normalized_structure(const Structure* original);
int check_for_duplicates(const Structure* normalized);
void remove_duplicate_atoms(Structure* original, const Structure* normalized);

#endif
