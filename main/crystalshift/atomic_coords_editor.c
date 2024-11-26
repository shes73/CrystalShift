#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "atomic_coords_editor.h"

#define M_PI 3.14159265358979323846

void add_vectors(Atom *atoms, int atom_count, double *vector) {
    for (int i = 0; i < atom_count; i++) {
        atoms[i].x += vector[0];
        atoms[i].y += vector[1];
        atoms[i].z += vector[2];
    }
}