#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "cell_filler.h"

#define M_PI 3.14159265358979323846

/// Парсит компоненту типа "-x+1/2", "y", "-z-1/4" и т.д.
double parse_operation_component(const char *component, double x, double y, double z) {
    double value = 0.0;
    int sign = 1;
    const char *ptr = component;

    // Пропускаем пробелы
    while (*ptr == ' ') ptr++;

    // Первый знак
    if (*ptr == '-') {
        sign = -1;
        ptr++;
    } else if (*ptr == '+') {
        ptr++;
    }

    // Основная переменная x/y/z
    if (*ptr == 'x') {
        value = x * sign;
        ptr++;
    } else if (*ptr == 'y') {
        value = y * sign;
        ptr++;
    } else if (*ptr == 'z') {
        value = z * sign;
        ptr++;
    } else {
        value = 0.0;  // не должно быть, но пусть будет
    }

    // Пропускаем пробелы
    while (*ptr == ' ') ptr++;

    // Дополнительный член: ±1/2 или ±0.25 и т.п.
    if (*ptr == '+' || *ptr == '-') {
        int second_sign = (*ptr == '+') ? 1 : -1;
        ptr++;

        while (*ptr == ' ') ptr++;

        int num = 0, denom = 1;
        if (sscanf(ptr, "%d/%d", &num, &denom) == 2) {
            value += second_sign * ((double)num / denom);
        } else {
            value += second_sign * atof(ptr);
        }
    }

    return value;
}

void fill_cell(Structure *structure) {
    if (structure == NULL || structure->symmetry_operations == NULL) {
        printf("No structure or symmetry operations provided.\n");
        return;
    }

    int op_count = structure->symmetry_operations_count;
    int original_atom_count = structure->atom_count;

    if (op_count == 0) {
        printf("No symmetry operations found.\n");
        return;
    }

    printf("Found %d symmetry operation%s:\n", op_count, op_count == 1 ? "" : "s");
    for (int i = 0; i < op_count; i++) {
        printf("  [%2d] %s\n", i + 1, structure->symmetry_operations[i]);
    }

    printf("Generating symmetrically equivalent atoms...\n");

    int total_expected = original_atom_count * op_count;
    structure->atoms = realloc(structure->atoms, sizeof(Atom) * total_expected);
    if (!structure->atoms) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    for (int i = 0; i < original_atom_count; i++) {
        Atom atom = structure->atoms[i];

        for (int j = 0; j < op_count; j++) {
            const char *op = structure->symmetry_operations[j];

            char op_x[64], op_y[64], op_z[64];
            if (sscanf(op, "%63[^,],%63[^,],%63s", op_x, op_y, op_z) != 3) {
                fprintf(stderr, "Failed to parse symmetry operation: %s\n", op);
                continue;
            }

            double new_x = parse_operation_component(op_x, atom.x, atom.y, atom.z);
            double new_y = parse_operation_component(op_y, atom.x, atom.y, atom.z);
            double new_z = parse_operation_component(op_z, atom.x, atom.y, atom.z);

            // Пропускаем j == 0 (тождественная операция), если не нужно дублировать
            if (j == 0) continue;

            Atom new_atom;
            new_atom.x = new_x;
            new_atom.y = new_y;
            new_atom.z = new_z;
            strncpy(new_atom.element, atom.element, sizeof(new_atom.element));
            new_atom.element[sizeof(new_atom.element) - 1] = '\0';

            structure->atoms[structure->atom_count++] = new_atom;
        }
    }

    printf("Final atom count after symmetry: %d\n", structure->atom_count);
    for (int i = 0; i < structure->atom_count; i++) {
        Atom *a = &structure->atoms[i];
        printf("%s: (%.6f, %.6f, %.6f)\n", a->element, a->x, a->y, a->z);
    }
}
