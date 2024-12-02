#ifndef POSCAR_PARSER_H
#define POSCAR_PARSER_H

#include "structures.h"

void calculate_lengths_and_angles(Structure *structure);

int parse_poscar(const char *filename);

#endif
