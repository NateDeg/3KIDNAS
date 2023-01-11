// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Table.h) - Source Finding Application                   //
// Copyright (C) 2022 The SoFiA 2 Authors                               //
// ____________________________________________________________________ //
//                                                                      //
// Address:  Tobias Westmeier                                           //
//           ICRAR M468                                                 //
//           The University of Western Australia                        //
//           35 Stirling Highway                                        //
//           Crawley WA 6009                                            //
//           Australia                                                  //
//                                                                      //
// E-mail:   tobias.westmeier [at] uwa.edu.au                           //
// ____________________________________________________________________ //
//                                                                      //
// This program is free software: you can redistribute it and/or modify //
// it under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or    //
// (at your option) any later version.                                  //
//                                                                      //
// This program is distributed in the hope that it will be useful,      //
// but WITHOUT ANY WARRANTY; without even the implied warranty of       //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         //
// GNU General Public License for more details.                         //
//                                                                      //
// You should have received a copy of the GNU General Public License    //
// along with this program. If not, see http://www.gnu.org/licenses/.   //
// ____________________________________________________________________ //
//                                                                      //

/// @file   Table.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Container class for reading tabulated data from text files (header).


#ifndef TABLE_H
#define TABLE_H

#include "common.h"

#define TABLE_INDEX(r,c) ((c) + self->cols * (r))  ///< Returns the array index of the table entry in row `r` and column `c`.
#define TABLE_MAX_LINE_SIZE 65536   ///< Maximum supported length of a single line in the input file from which the table data are read.

typedef CLASS Table Table;


// ----------------------------------------------------------------- //
// Class 'Table'                                                     //
// ----------------------------------------------------------------- //
// This class provides a convenient container for reading tabulated  //
// data from an ASCII file. The container size will be automatically //
// adjusted to match the number of rows and columns in the file. The //
// data values will be stored as double-precision floating-point.    //
// The standard constructor has been set to private to disallow the  //
// creation of empty Table objects.                                  //
// ----------------------------------------------------------------- //

// Constructor and destructor
// NOTE: The standard constructor is private, as Table objects should
//       only be created using the Table_from_file constructor.
PRIVATE Table  *Table_new       (void);
PUBLIC  Table  *Table_from_file (const char *filename, const char *delimiters);
PUBLIC  void    Table_delete    (Table * self);

// Public methods
PUBLIC  size_t  Table_rows      (const Table *self);
PUBLIC  size_t  Table_cols      (const Table *self);
PUBLIC  double  Table_get       (const Table *self, const size_t row, const size_t col);
PUBLIC  void    Table_set       (Table *self, const size_t row, const size_t col, const double value);

#endif
