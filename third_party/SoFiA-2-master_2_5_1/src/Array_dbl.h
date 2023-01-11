// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Array_dbl.h) - Source Finding Application               //
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

/// @file   Array_dbl.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Container class template for storing a dynamic array of values of type `double` (header).


// WARNING: This is a template that needs to be instantiated before use.
//          Do not edit template instances, as they are auto-generated
//          and will be overwritten during instantiation!


#ifndef ARRAY_dbl_H
#define ARRAY_dbl_H

#include "common.h"


// ----------------------------------------------------------------- //
// Class 'Array_dbl'                                                 //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a convenient way to store //
// multiple values of a specific type in an array-like structure. A  //
// new array can either be of a given size and empty (using the      //
// standard constructor) or provided with a list of comma-separated  //
// values that will be stored in the array and used to determine its //
// size (using the alternative constructor).                         //
// ----------------------------------------------------------------- //

typedef CLASS Array_dbl Array_dbl;

// Constructor and destructor
PUBLIC Array_dbl    *Array_dbl_new      (const size_t size);
PUBLIC Array_dbl    *Array_dbl_new_str  (const char *string);
PUBLIC Array_dbl    *Array_dbl_copy     (const Array_dbl *source);
PUBLIC void          Array_dbl_delete   (Array_dbl *self);

// Public methods
PUBLIC const double *Array_dbl_get_ptr  (const Array_dbl *self);
PUBLIC size_t        Array_dbl_get_size (const Array_dbl *self);
PUBLIC Array_dbl    *Array_dbl_push     (Array_dbl *self, const double value);
PUBLIC double        Array_dbl_get      (const Array_dbl *self, const size_t index);
PUBLIC Array_dbl    *Array_dbl_set      (Array_dbl *self, const size_t index, const double value);
PUBLIC Array_dbl    *Array_dbl_add      (Array_dbl *self, const size_t index, const double value);
PUBLIC Array_dbl    *Array_dbl_mul      (Array_dbl *self, const size_t index, const double value);
PUBLIC Array_dbl    *Array_dbl_cat      (Array_dbl *self, const Array_dbl *source);
PUBLIC Array_dbl    *Array_dbl_sort     (Array_dbl *self);

#endif
