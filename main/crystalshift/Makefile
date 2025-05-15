CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c99
LDFLAGS=-lm

TARGET = crystalshift

SRC = main.c structures.c atomic_coords_editor.c basis_changer.c cif_parser.c cif_writer.c poscar_parser.c poscar_writer.c xyz_parser.c xyz_writer.c cell_filler.c
HEADERS = structures.h atomic_coords_editor.h basis_changer.h cif_parser.h cif_writer.h poscar_parser.h poscar_writer.h xyz_parser.h xyz_writer.h cell_filler.h
OBJ = $(SRC:.c=.o)

.PHONY: all clean rebuild

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(TARGET)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

rebuild: clean all

