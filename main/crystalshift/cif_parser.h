#ifndef CIF_PARSER_H
#define CIF_PARSER_H

char *strdup(const char *s); // i've got warning bc compiler doesn't recognise strdup from string.h...

void calculate_lattice_matrix();

void remove_errors(char *line);

void parse_symmetry_operations_from_file(FILE *file);

int parse_cif(const char *filename);

#endif
