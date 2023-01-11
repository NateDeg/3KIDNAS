#!/bin/sh
# ____________________________________________________________________ #
#                                                                      #
# SoFiA 2.5.1 (compile.sh) - Source Finding Application                #
# Copyright (C) 2022 The SoFiA 2 Authors                               #
# ____________________________________________________________________ #
#                                                                      #
# Address:  Tobias Westmeier                                           #
#           ICRAR M468                                                 #
#           The University of Western Australia                        #
#           35 Stirling Highway                                        #
#           Crawley WA 6009                                            #
#           Australia                                                  #
#                                                                      #
# E-mail:   tobias.westmeier [at] uwa.edu.au                           #
# ____________________________________________________________________ #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see http://www.gnu.org/licenses/.   #
# ____________________________________________________________________ #
#                                                                      #

# Usage: ./compile.sh [-fopenmp]
#
# The optional argument -fopenmp can be supplied to enable multi-threading.
# By default multi-threading will be disabled to allow compilation on Mac
# with Clang.

echo "_______________________________________________________________________"
echo
echo " Installing SoFiA"
echo "_______________________________________________________________________"
echo

# Compile source files
echo "  Compiling src/common.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/common.o -c src/common.c
echo "  Compiling src/statistics_flt.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/statistics_flt.o -c src/statistics_flt.c
echo "  Compiling src/statistics_dbl.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/statistics_dbl.o -c src/statistics_dbl.c
echo "  Compiling src/Table.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Table.o -c src/Table.c
echo "  Compiling src/String.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/String.o -c src/String.c
echo "  Compiling src/Stack.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Stack.o -c src/Stack.c
echo "  Compiling src/Path.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Path.o -c src/Path.c
echo "  Compiling src/Array_dbl.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Array_dbl.o -c src/Array_dbl.c
echo "  Compiling src/Array_siz.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Array_siz.o -c src/Array_siz.c
echo "  Compiling src/Map.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Map.o -c src/Map.c
echo "  Compiling src/Matrix.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Matrix.o -c src/Matrix.c
echo "  Compiling src/LinkerPar.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/LinkerPar.o -c src/LinkerPar.c $1
echo "  Compiling src/Parameter.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Parameter.o -c src/Parameter.c
echo "  Compiling src/Source.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Source.o -c src/Source.c
echo "  Compiling src/Catalog.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Catalog.o -c src/Catalog.c
echo "  Compiling src/Flagger.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Flagger.o -c src/Flagger.c
echo "  Compiling src/WCS.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/WCS.o -c src/WCS.c
echo "  Compiling src/Header.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/Header.o -c src/Header.c
echo "  Compiling src/DataCube.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o src/DataCube.o -c src/DataCube.c $1
echo "  Compiling sofia.c"
gcc --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O3 -o sofia src/common.o src/statistics_flt.o src/statistics_dbl.o src/Table.o src/String.o src/Stack.o src/Path.o src/Array_dbl.o src/Array_siz.o src/Map.o src/Matrix.o src/LinkerPar.o src/Parameter.o src/Flagger.o src/WCS.o src/Header.o src/DataCube.o src/Source.o src/Catalog.o sofia.c -lm -lwcs $1

# Remove object files
#rm -rf src/*.o

# Print instructions
echo "_______________________________________________________________________"
echo
echo " Installation complete"
echo "_______________________________________________________________________"
echo
echo "  Please check above for any error messages  produced by the compiler"
echo "  before proceeding with the instructions below. If no error messages"
echo "  have occurred, you can choose to create a symbolic link in /usr/bin"
echo "  to make SoFiA available across your system:"
echo
echo "    sudo ln -s ${PWD}/sofia /usr/bin/sofia"
echo
echo "  Alternatively, if you don't have root privileges on your system, an"
echo "  alias can be set in your .bashrc or .cshrc configuration file. Once"
echo "  complete, SoFiA can then be invoked via 'sofia <parameter_file>'."
echo
