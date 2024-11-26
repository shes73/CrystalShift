#ifndef XYZ_PARSER_H
#define XYZ_PARSER_H

void inverse_matrix(double transformation_matrix[3][3], double inversed_transformation_matrix[3][3]);

int parse_xyz(const char *filename, int xyz_format);

#endif