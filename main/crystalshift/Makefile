CC=gcc
CFLAGS=-Wall -Wextra -Werror

TARGET = CrystalEngineer

SRC = main.c layers_mover.c basis_changer.c format_converter.c
HEADERS = layers_mover.h basis_changer.h format_converter.h
OBJ = $(SRC:.c=.o)

.PHONY: all clean rebuild

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean: all
	rm -f $(OBJ) $(TARGET)

rebuild: clean all
