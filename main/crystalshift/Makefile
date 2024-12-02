CC=gcc
CFLAGS=-Wall -Wextra -Werror

TARGET = crystalshift

SRC = main.c structures.c atomic_coords_editor.c basis_changer.c cif_parser.c cif_writer.c poscar_parser.c poscar_writer.c xyz_parser.c xyz_writer.c
HEADERS = structures.h atomic_coords_editor.h basis_changer.h cif_parser.h cif_writer.h poscar_parser.h poscar_writer.h xyz_parser.h xyz_writer.h
OBJ = $(SRC:.c=.o)

.PHONY: all clean rebuild

all: gcc $(FLAGS) -std=c99 $(SRC) -lm -o $(TARGET)

clean: all
	rm -f $(OBJ) $(TARGET)

rebuild: clean all
