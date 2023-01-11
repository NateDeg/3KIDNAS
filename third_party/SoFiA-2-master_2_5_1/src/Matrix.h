// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Matrix.h) - Source Finding Application                  //
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

/// @file   Matrix.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for handling matrices and enabling basic matrix manipulation (header).


#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include "common.h"


// ----------------------------------------------------------------- //
// Class 'Matrix'                                                    //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a way of storing and      //
// handling matrices.                                                //
// ----------------------------------------------------------------- //

typedef CLASS Matrix Matrix;

// Constructor and destructor
PUBLIC  Matrix       *Matrix_new        (const size_t rows, const size_t cols);  // Standard constructor
PUBLIC  Matrix       *Matrix_copy       (const Matrix *source);                  // Copy constructor
PUBLIC  Matrix       *Matrix_identity   (const size_t size);                     // Constructor for square identity matrix
PUBLIC  Matrix       *Matrix_covar      (const size_t size, const size_t samples, const double *values); // Constructor for covariance matrix
PUBLIC  void          Matrix_delete     (Matrix *self);

// Public methods
PUBLIC  size_t        Matrix_rows       (const Matrix *self);
PUBLIC  size_t        Matrix_cols       (const Matrix *self);
PUBLIC  void          Matrix_set_value  (Matrix *self, const size_t row, const size_t col, const double value);
PUBLIC  void          Matrix_set_value_nocheck(Matrix *self, const size_t row, const size_t col, const double value);
PUBLIC  double        Matrix_get_value  (const Matrix *self, const size_t row, const size_t col);
PUBLIC  double        Matrix_get_value_nocheck(const Matrix *self, const size_t row, const size_t col);
PUBLIC  void          Matrix_add_value  (Matrix *self, const size_t row, const size_t col, const double value);
PUBLIC  void          Matrix_mul_value  (Matrix *self, const size_t row, const size_t col, const double value);
PUBLIC  void          Matrix_mul_scalar (Matrix *self, const double scalar);
PUBLIC  Matrix       *Matrix_mul_matrix (const Matrix *self, const Matrix *matrix);
PUBLIC  void          Matrix_add_matrix (Matrix *self, const Matrix *matrix);
PUBLIC  double        Matrix_vMv        (const Matrix *self, const Matrix *vector);
PUBLIC  double        Matrix_vMv_nocheck(const Matrix *self, const Matrix *vector);
PUBLIC  Matrix       *Matrix_transpose  (const Matrix *self);
PUBLIC  Matrix       *Matrix_invert     (const Matrix *self);
PUBLIC  void          Matrix_print      (const Matrix *self, const unsigned int width, const unsigned int decimals);
PUBLIC  double        Matrix_det        (const Matrix *self, const double scale_factor);
PUBLIC  double        Matrix_prob_dens  (const Matrix *covar_inv, const Matrix *vector, const double scal_fact);
PUBLIC  double        Matrix_prob_dens_nocheck(const Matrix *covar_inv, const Matrix *vector, const double scal_fact);
PUBLIC  void          Matrix_err_ellipse(const Matrix *covar, const size_t par1, const size_t par2, double *radius_maj, double *radius_min, double *pa);
//PUBLIC  void          Matrix_covariance (Matrix *self, const double values[], const size_t dim, const size_t length);

// Private methods
PRIVATE inline size_t Matrix_get_index  (const Matrix *self, const size_t row, const size_t col);
PRIVATE void          Matrix_swap_rows  (Matrix *self, const size_t row1, const size_t row2);
PRIVATE void          Matrix_add_row    (Matrix *self, const size_t row1, const size_t row2, const double factor);
PRIVATE void          Matrix_mul_row    (Matrix *self, const size_t row, const double factor);

#endif
