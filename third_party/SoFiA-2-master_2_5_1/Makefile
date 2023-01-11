# Makefile for the compilation of SoFiA 2
#
# Usage examples:
#   make                               for GCC or Clang without OpenMP
#   make OMP=-fopenmp                  for GCC or Clang with OpenMP
#   make CC=icc OPT=-O3 OMP=-openmp    for Intel C Compiler with OpenMP (not tested)
#   make clean                         remove object files after compilation
#   make DEBUG=1                       for debug mode (no compiler optimisations)


SRC = src/Array_dbl.c \
      src/Array_siz.c \
      src/Catalog.c \
      src/common.c \
      src/DataCube.c \
      src/Flagger.c \
      src/Header.c \
      src/LinkerPar.c \
      src/Map.c \
      src/Matrix.c \
      src/Parameter.c \
      src/Path.c \
      src/Source.c \
      src/Stack.c \
      src/statistics_dbl.c \
      src/statistics_flt.c \
      src/String.c \
      src/Table.c \
      src/WCS.c

OBJ = $(SRC:.c=.o)

TEST = tests/test_LinkerPar.c

TEST_OBJ = $(TEST:.c=.o)

# OPENMP = -fopenmp
OMP     =
OPT     = --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3
LIBS    = -lm -lwcs
CC      = gcc
CFLAGS += $(OPT) $(OMP)

ifdef DEBUG
OPT     = -g -O0
endif

all:	sofia

sofia:	$(OBJ)
	$(CC) $(CFLAGS) -o sofia sofia.c $(OBJ) $(LIBS)

unittest:	$(OBJ) $(TEST_OBJ)
	$(CC) $(CFLAGS) -o unittest tests/unittest.c $(TEST_OBJ) $(OBJ) $(LIBS) `pkg-config --cflags --libs check`

clean:
	rm -rf $(OBJ) $(TEST_OBJ)
