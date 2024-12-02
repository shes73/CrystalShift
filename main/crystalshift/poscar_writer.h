#ifndef POSCAR_WRITER_H
#define POSCAR_WRITER_H

int get_atomic_number(const char *symbol);
int compare_atoms(const void *a, const void *b, int ascending);
int compare_atoms_asc(const void *a, const void *b);
int compare_atoms_desc(const void *a, const void *b);
void sort_atoms_custom(char **custom_order, int custom_count);

void write_poscar(const char *filename, int sort_option, char **custom_order, int custom_count);

#endif
