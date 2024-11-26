#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basis_changer.h"
#include "atomic_coords_editor.h"
#include "structures.h"
#include "cif_parser.h"
#include "poscar_parser.h"
#include "xyz_parser.h"
#include "cif_writer.h"
#include "poscar_writer.h"
#include "xyz_writer.h"

Structure structure;
Structure structure = {.atom_capacity = 100}; 

int main() {
    char filename[256];
    char poscar_filename[256];
    char xyz_filename[256];
    char output_filename[256];
    int format_type;
    int xyz_format;
    char order_input[1024];
    char *custom_order[256];
    int custom_count = 0;
    int sort_option = 1;

    int option;
    printf("Choose the option you need:\n"
    "1. Change basis of unit cell\n"
    "2. Edit atomic coordinates\n"
    "3. Converter CIF -> POSCAR\n"
    "4. Converter CIF -> xyz\n"
    "5. Converter POSCAR -> CIF\n"
    "6. Converter POSCAR -> xyz\n");
    scanf("%d", &option);

    switch (option) {
        case 1:
            printf("Change basis of unit cell\n");

            printf("Please, choose the format of the input file: \n");
            printf("1. CIF\n2. POSCAR\n");
            scanf("%d", &format_type);

            if (format_type == 1) {
                printf("Enter the path to the CIF file: ");
                scanf("%s", filename);
                if (parse_cif(filename) != 0) {
                    printf("Error parsing CIF file!\n");
                    return 1;
                }
            } else if (format_type == 2) {
                printf("Enter the path to the POSCAR file: ");
                scanf("%s", poscar_filename);
                if (parse_poscar(poscar_filename) != 0) {
                    printf("Error parsing POSCAR file!\n");
                    return 1;
                }
            } else {
                printf("Error! CrystalShift can't work with this file format!\n");
                return 1;
            }

            double transformation_matrix[3][3];
            double inversed_transformation_matrix[3][3];
            double new_basis_matrix[3][3];
            get_the_transformation_matrix(transformation_matrix);
            create_supercell_if_needed(transformation_matrix);
            inverse_transformation_matrix(transformation_matrix, inversed_transformation_matrix);
            double new_cell_lengths[3], new_cell_angles[3];
            calculate_new_basis_matrix(transformation_matrix, new_basis_matrix, new_cell_lengths, new_cell_angles);
            calculate_new_coords(inversed_transformation_matrix);
            
            if (format_type == 1) {
                printf("Enter the path to save the new CIF file: ");
                scanf("%s", output_filename);
                write_cif(output_filename);
            } else if (format_type == 2) {
                printf("Choose the order for writing atoms to POSCAR:\n");
                printf("1. From lightest to heaviest;\n");
                printf("2. From heaviest to lightest;\n");
                printf("3. Custom order (input symbols separated by spaces).\n");
                printf("Enter your choice: ");
                scanf("%d", &sort_option);
                getchar();

                if (sort_option == 3) {
                    printf("Enter custom order: ");
                    fgets(order_input, sizeof(order_input), stdin);

                    char *token = strtok(order_input, " \n");
                    while (token != NULL) {
                        custom_order[custom_count++] = token;
                        token = strtok(NULL, " \n");
                    }
                } else if (sort_option <= 0 || sort_option > 3) {
                    printf("Error! There's no such option!\n");
                    return 1;
                }
                printf("Enter the path to save the new POSCAR file: ");
                scanf("%s", output_filename);
                write_poscar(poscar_filename, sort_option, custom_order, custom_count);
            }

            break;

        case 2:
            printf("Edit atomic coordinates in the unit cell\n");

            double vector[3];

            printf("Please, choose the format of the input file: \n");
            printf("1. CIF\n2. POSCAR\n");
            scanf("%d", &format_type);

            if (format_type == 1) {
                printf("Enter the path to the CIF file: ");
                scanf("%s", filename);
                if (parse_cif(filename) != 0) {
                    printf("Error parsing CIF file!\n");
                    return 1;
                }
            } else if (format_type == 2) {
                printf("Enter the path to the POSCAR file: ");
                scanf("%s", poscar_filename);
                if (parse_poscar(poscar_filename) != 0) {
                    printf("Error parsing POSCAR file!\n");
                    return 1;
                }
            } else {
                printf("Error! CrystalShift can't work with this file format!\n");
                return 1;
            }

            printf("Enter the three float numbers separated by spaces that will be added to the coordinates of all atoms in the file (e.g. 0.25 0.0 0.0): ");
            scanf("%lf %lf %lf", &vector[0], &vector[1], &vector[2]);
            add_vectors(structure.atoms, structure.atom_count, vector);

            if (format_type == 1) {
                printf("Enter the path to save the new CIF file: ");
                scanf("%s", output_filename);
                write_cif(output_filename);
            } else if (format_type == 2) {
                printf("Choose the order for writing atoms to POSCAR:\n");
                printf("1. From lightest to heaviest;\n");
                printf("2. From heaviest to lightest;\n");
                printf("3. Custom order (input symbols separated by spaces).\n");
                printf("Enter your choice: ");
                scanf("%d", &sort_option);
                getchar();

                if (sort_option == 3) {
                    printf("Enter custom order: ");
                    fgets(order_input, sizeof(order_input), stdin);

                    char *token = strtok(order_input, " \n");
                    while (token != NULL) {
                        custom_order[custom_count++] = token;
                        token = strtok(NULL, " \n");
                    }
                } else if (sort_option <= 0 || sort_option > 3) {
                    printf("Error! There's no such option!\n");
                    return 1;
                }
                printf("Enter the path to save the new POSCAR file: ");
                scanf("%s", output_filename);
                write_poscar(poscar_filename, sort_option, custom_order, custom_count);
            }

            free(structure.atoms);
            break;

        case 3:
            printf("Converter CIF -> POSCAR\n");

            printf("Enter the path to the CIF file: ");
            scanf("%s", filename);

            if (parse_cif(filename) != 0) {
                printf("Error parsing CIF file!\n");
                return 1;
            }

            printf("Choose the order for writing atoms to POSCAR:\n");
            printf("1. From lightest to heaviest;\n");
            printf("2. From heaviest to lightest;\n");
            printf("3. Custom order (input symbols separated by spaces).\n");
            printf("Enter your choice: ");
            scanf("%d", &sort_option);
            getchar();

            if (sort_option == 3) {
                printf("Enter custom order: ");
                fgets(order_input, sizeof(order_input), stdin);

                char *token = strtok(order_input, " \n");
                while (token != NULL) {
                    custom_order[custom_count++] = token;
                    token = strtok(NULL, " \n");
                }
            } else if (sort_option <= 0 || sort_option > 3) {
                printf("Error! There's no such option!\n");
                return 1;
            }

            printf("Enter the path to the POSCAR file: ");
            scanf("%s", poscar_filename);
            write_poscar(poscar_filename, sort_option, custom_order, custom_count);

            free(structure.atoms);
            break;

        case 4:
            printf("Converter CIF -> xyz\n");

            printf("Enter the path to the CIF file: ");
            scanf("%s", filename);

            if (parse_cif(filename) != 0) {
                printf("Error parsing CIF file!\n");
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
            break;

        case 5:
            printf("Converter POSCAR -> CIF\n");

            printf("Enter the path to the POSCAR file: ");
            scanf("%s", poscar_filename);

            if (parse_poscar(poscar_filename) != 0) {
                printf("Error parsing POSCAR file!\n");
                return 1;
            }

            printf("Enter the path to the CIF file: ");
            scanf("%s", filename);
            write_cif(filename);

            free(structure.atoms);
            break;

        case 6:
            printf("Converter POSCAR -> xyz\n");

            printf("Enter the path to the POSCAR file: ");
            scanf("%s", poscar_filename);

            if (parse_poscar(poscar_filename) != 0) {
                printf("Error parsing POSCAR file!\n");
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
            break;

        default:
            printf("Sorry! Non-existent option!\n");
            break;
    }

    return 0;
}
