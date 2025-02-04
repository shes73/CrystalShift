#ifndef BASIS_CHANGER_H
#define BASIS_CHANGER_H

void get_the_transformation_matrix(double transformation_matrix[3][3]);
void swap_rows(double matrix[3][3], double inverse[3][3], int row1, int row2);
void inverse_transformation_matrix(double transformation_matrix[3][3], double inversed_transformation_matrix[3][3]);
void calculate_new_basis_matrix(double transformation_matrix[3][3], double new_basis_matrix[3][3], double new_cell_lengths[3], double new_cell_angles[3]);
void calculate_new_coords(double inversed_transformation_matrix[3][3]);
void create_supercell_if_needed(double transformation_matrix[3][3]);

#endif
