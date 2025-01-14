CC=gcc
CFLAGS=-Wall -Wextra -Werror -std=c99
LDFLAGS=-lm

TARGET = crystalshift_layers

SRC = crystalshift_layers.c
OBJ = $(SRC:.c=.o)

.PHONY: all clean rebuild

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)

rebuild: clean all

