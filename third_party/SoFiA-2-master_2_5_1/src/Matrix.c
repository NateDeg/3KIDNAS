// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Matrix.c) - Source Finding Application                  //
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

/// @file   Matrix.c
/// @author Tobias Westmeier
/// @date   23/11/2021
/// @brief  Class for handling matrices and enabling basic matrix manipulation.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Array_dbl.h"



/// @brief Class for handling matrices
///
/// The purpose of this class is to provide a way of storing and
/// handling matrices.

CLASS Matrix
{
	size_t  rows;    ///< Number of matrix rows.
	size_t  cols;    ///< Number of matrix columns.
	double *values;  ///< Array containing the values of the matrix.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new, empty Matrix object and
/// return a pointer to the newly created object. The number of rows
/// and columns of the matrix needs to be specified. The matrix will
/// be initialised with 0. Note that the destructor will need to be
/// called explicitly once the object is no longer required to release
/// any memory allocated during the lifetime of the object.
///
/// @param rows  Number of matrix rows.
/// @param cols  Number of matrix columns.
///
/// @return Pointer to newly created Matrix object.

PUBLIC Matrix *Matrix_new(const size_t rows, const size_t cols)
{
	// Sanity checks
	ensure(rows && cols, ERR_USER_INPUT, "Number of matrix rows and cols must be > 0.");
	
	Matrix *self = (Matrix *)memory(MALLOC, 1, sizeof(Matrix));
	
	self->rows = rows;
	self->cols = cols;
	
	self->values = (double *)memory(CALLOC, rows * cols, sizeof(double));
	
	return self;
}



/// @brief Copy constructor
///
/// Copy constructor. Will create a new matrix with the same dimensions
/// as the specified source matrix and copy all entries from `source`. A
/// pointer to the newly created copy will be returned. Note that the
/// destructor will need to be called on the copy if it is no longer
/// required to release its memory.
///
/// @param source  Original matrix to be copied.
///
/// @return Pointer to newly created copy of matrix.

PUBLIC Matrix *Matrix_copy(const Matrix *source)
{
	// Sanity checks
	check_null(source);
	
	// Call standard constructor
	Matrix *self = Matrix_new(source->rows, source->cols);
	
	// Copy values
	memcpy(self->values, source->values, source->rows * source->cols * sizeof(double));
	
	return self;
}



/// @brief Constructor for square identity matrix
///
/// Variant of the standard constructor that will create a square
/// identity matrix where all diagonal elements are set to 1. A
/// pointer to the newly created matrix will be returned. Note that
/// the destructor will need to be called on the matrix if it is no
/// longer required to release its memory.
///
/// @param size  Size of the new matrix.
///
/// @return Pointer to newly created identity matrix.

PUBLIC Matrix *Matrix_identity(const size_t size)
{
	// Sanity checks
	ensure(size, ERR_USER_INPUT, "Matrix size must be > 0.");
	
	// Call standard constructor
	Matrix *self = Matrix_new(size, size);
	
	// Set diagonal elements to 1
	for(size_t i = 0; i < size; ++i) self->values[Matrix_get_index(self, i, i)] = 1.0;
	
	return self;
}


/// @brief Constructor for covariance matrix
///
/// Variant of the standard constructor that will create the covariance
/// matrix from a flattened array of values specified by the user. A
/// pointer to the newly created covariance will be returned. Note that
/// the destructor will need to be called on the matrix if it is no
/// longer required to release its memory.
///
/// The specified array of values must be a flattened array of `size`
/// columns representing the different variables or parameters for which
/// the covariance matrix is to be calculated and `samples` rows which
/// represent the list of samples supplied to the covariance calculation.
/// Hence, the array of `values` must be of length `size * samples`.
///
/// @param size      Size of covariance matrix. Must be equal to the number
///                  of columns in the array prior to flattening (i.e.
///                  number of variables or dimensions).
/// @param samples   Number of rows in the array prior to flattening, i.e.
///                  the number of samples.
/// @param values    Flattened array of values. Must be of length `size *
///                  samples`.
///
/// @return Pointer to newly created covariance matrix.

PUBLIC Matrix *Matrix_covar(const size_t size, const size_t samples, const double *values)
{
	// Create empty square matrix
	Matrix *covar = Matrix_new(size, size);
	
	// Calculate mean values
	Array_dbl *mean = Array_dbl_new(size);
	for(size_t i = size; i--;)
	{
		for(size_t j = 0; j < samples; ++j) Array_dbl_add(mean, i, values[size * j + i]);
		Array_dbl_mul(mean, i, 1.0 / samples);
	}
	
	// Then calculate the covariance matrix
	for(size_t i = size; i--;)
	{
		for(size_t j = size; j--;)
		{
			for(size_t k = 0; k < size * samples; k += size) Matrix_add_value(covar, i, j, (values[k + i] - Array_dbl_get(mean, i)) * (values[k + j] - Array_dbl_get(mean, j)));
			Matrix_mul_value(covar, i, j, 1.0 / samples);
		}
	}
	
	// Clean up
	Array_dbl_delete(mean);
	
	return covar;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Matrix_delete(Matrix *self)
{
	if(self != NULL) free(self->values);
	free(self);
	
	return;
}



/// @brief Returns number of rows
///
/// Public method for returning the number of rows of the specified
/// matrix.
///
/// @param self  Object self-reference.
///
/// @return Number of matrix rows.

PUBLIC size_t Matrix_rows(const Matrix *self)
{
	check_null(self);
	return self->rows;
}



/// @brief Return number of columns
///
/// @param self  Object self-reference.
///
/// @return Number of matrix columns.
///
/// Public method for returning the number of columns of the
/// specified matrix.

PUBLIC size_t Matrix_cols(const Matrix *self)
{
	check_null(self);
	return self->cols;
}



/// @brief Set matrix element to specified value
///
/// Public method for setting the element at position (`row`, `col`)
/// to the specified value. If `row` or `col` are out of range, the
/// process will be terminated. See Matrix_set_value_nocheck() for a
/// faster version without sanity checks.
///
/// @param self   Object self-reference.
/// @param row    Row of the element to be set.
/// @param col    Column of the element to be set.
/// @param value  Value to set the element to.

PUBLIC void Matrix_set_value(Matrix *self, const size_t row, const size_t col, const double value)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Matrix row or col out of range.");
	
	self->values[Matrix_get_index(self, row, col)] = value;
	
	return;
}



/// @brief Set matrix element to specified value (no sanity checks)
///
/// Public method for setting the element at position (`row`, `col`)
/// to the specified value. If `row` or `col` are out of range, the
/// process will be terminated. Identical to Matrix_set_value(), but
/// without sanity checks for better performance.
///
/// @param self   Object self-reference.
/// @param row    Row of the element to be set.
/// @param col    Column of the element to be set.
/// @param value  Value to set the element to.

PUBLIC void Matrix_set_value_nocheck(Matrix *self, const size_t row, const size_t col, const double value)
{
	self->values[Matrix_get_index(self, row, col)] = value;
	return;
}



/// @brief Get matrix element at specified position
///
/// Public method for retrieving the value of the matrix at the
/// position specified by `row` and `col`. If `row` or `col`
/// are out of range, the process will be terminated. Also see
/// Matrix_get_value_nocheck() for a version without sanity
/// checks for better performance.
///
/// @param self  Object self-reference.
/// @param row   Row of the element to be retrieved.
/// @param col   Column of the element to be retrieved.
///
/// @return Value of the matrix element at (row, col).

PUBLIC double Matrix_get_value(const Matrix *self, const size_t row, const size_t col)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Matrix row or col out of range.");
	
	return self->values[Matrix_get_index(self, row, col)];
}



/// @brief Get matrix element at specified position (no sanity checks)
///
/// Public method for retrieving the value of the matrix at the
/// position specified by `row` and `col`. If `row` or `col`
/// are out of range, the process will be terminated. Identical
/// to Matrix_get_value(), but without sanity checks for better
/// performance.
///
/// @param self  Object self-reference.
/// @param row   Row of the element to be retrieved.
/// @param col   Column of the element to be retrieved.
///
/// @return Value of the matrix element at (row, col).

PUBLIC double Matrix_get_value_nocheck(const Matrix *self, const size_t row, const size_t col)
{
	return self->values[Matrix_get_index(self, row, col)];
}



/// @brief Add specified value to matrix element
///
/// Public method for adding the specified value to the matrix element
/// specified by (`row`, `col`). If `row` or `col` are out of range,
/// the process will be terminated.
///
/// @param self   Object self-reference.
/// @param row    Row of the element to be set.
/// @param col    Column of the element to be set.
/// @param value  Value to add to the element.

PUBLIC void Matrix_add_value(Matrix *self, const size_t row, const size_t col, const double value)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Matrix row or col out of range.");
	
	self->values[Matrix_get_index(self, row, col)] += value;
	
	return;
}



/// @brief Multiply matrix element by specified value
///
/// Public method for multiplying the matrix element specified by
/// (`row`, `col`) by value. If `row` or `col` are out of range,
/// the process will be terminated.
///
/// @param self   Object self-reference.
/// @param row    Row of the element to be set.
/// @param col    Column of the element to be set.
/// @param value  Value to multiply element by.

PUBLIC void Matrix_mul_value(Matrix *self, const size_t row, const size_t col, const double value)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Matrix row or col out of range.");
	
	self->values[Matrix_get_index(self, row, col)] *= value;
	
	return;
}



/// @brief Multiply matrix by scalar
///
/// Public method for multiplying the matrix with the specified
/// scalar value. The matrix will be multiplied in situ.
///
/// @param self    Object self-reference.
/// @param scalar  Scalar to multiply the matrix by.

PUBLIC void Matrix_mul_scalar(Matrix *self, const double scalar)
{
	// Sanity checks
	check_null(self);
	
	for(double *ptr = self->values + self->rows * self->cols; ptr --> self->values;) *ptr *= scalar;
	
	return;
}



/// @brief Matrix multiplication
///
/// Public method for multiplying a matrix with another one. The
/// result of the multiplication will be written to a new matrix
/// that will be returned by the method. For this to work, the number
/// of columns in the left operand (self) must be equal to the number
/// of rows in the right operand (matrix). If not, the process will
/// be terminated. The user is responsible for calling the destructor
/// on the returned result matrix once it is no longer needed.
///
/// @param self    Object self-reference.
/// @param matrix  Matrix by which to multiply.
///
/// @return Result of the matrix multiplication.

PUBLIC Matrix *Matrix_mul_matrix(const Matrix *self, const Matrix *matrix)
{
	// Sanity checks
	check_null(self);
	check_null(matrix);
	ensure(self->cols == matrix->rows, ERR_USER_INPUT, "Incompatible row and column numbers in matrix multiplication.");
	
	// Create result matrix
	Matrix *result = Matrix_new(self->rows, matrix->cols);
	
	// Calculate entries
	for(size_t i = 0; i < result->rows; ++i)
	{
		for(size_t k = 0; k < result->cols; ++k)
		{
			double value = 0.0;
			for(size_t j = 0; j < self->cols; ++j) value += self->values[Matrix_get_index(self, i, j)] * matrix->values[Matrix_get_index(matrix, j, k)];
			result->values[Matrix_get_index(result, i, k)] += value;
		}
	}
	
	return result;
}



/// @brief Add two matrices
///
/// Public method for adding one matrix to another. Both matrices
/// must have the same size, otherwise the process will be terminated.
/// The input matrix will be modified in situ.
///
/// @param self    Object self-reference.
/// @param matrix  Matrix to be added.

PUBLIC void Matrix_add_matrix(Matrix *self, const Matrix *matrix)
{
	// Sanity checks
	check_null(self);
	check_null(matrix);
	ensure(self->rows == matrix->rows && self->cols == matrix->cols, ERR_USER_INPUT, "Incompatible row and column numbers in matrix addition.");
	
	// Calculate entries
	for(size_t i = 0; i < self->rows; ++i)
	{
		for(size_t j = 0; j < self->cols; ++j)
		{
			self->values[Matrix_get_index(self, i, j)] += matrix->values[Matrix_get_index(matrix, i, j)];
		}
	}
	
	return;
}


/// @brief Calculate v^T M v
///
/// Public method for calculating the product of v^T M v between
/// the specified matrix, M, and vector, v, where v^T denotes the
/// transpose of the vector. The matrix must be a square matrix
/// with a size equal to that of the vector. The vector is expected
/// to be a column vector. Also see Matrix_vMv_nocheck() for a
/// version without sanity checks for better performance.
///
/// @param self    Object self-reference.
/// @param vector  Vector to be multiplied.
///
/// @return Result of v^T M v (which is a scalar).

PUBLIC double Matrix_vMv(const Matrix *self, const Matrix *vector)
{
	// Sanity checks
	check_null(self);
	check_null(vector);
	ensure(self->rows == self->cols, ERR_USER_INPUT, "Matrix is not square.");
	ensure(vector->cols == 1, ERR_USER_INPUT, "Vector has more than one column.");
	ensure(self->rows == vector->rows, ERR_USER_INPUT, "Vector size (%zu) does not match matrix (%zu x %zu).", vector->rows, self->rows, self->cols);
	
	double *array = (double *)memory(CALLOC, self->rows, sizeof(double));
	
	for(size_t col = self->cols; col--;)
	{
		for(size_t row = self->rows; row--;)
		{
			*(array + col) += Matrix_get_value_nocheck(vector, row, 0) * Matrix_get_value_nocheck(self, row, col);
		}
	}
	
	double result = 0.0;
	for(size_t i = self->rows; i--;) result += *(array + i) * Matrix_get_value_nocheck(vector, i, 0);
	
	free(array);
	
	return result;
}



/// @brief Calculate v^T M v (no sanity checks)
///
/// Public method for calculating the product of v^T M v between
/// the specified matrix, M, and vector, v, where v^T denotes the
/// transpose of the vector. The matrix must be a square matrix
/// with a size equal to that of the vector. The vector is expected
/// to be a column vector. Identical to Matrix_vMv(), but without
/// sanity checks for better performance.
///
/// @param self    Object self-reference.
/// @param vector  Vector to be multiplied.
///
/// @return Result of v^T M v (which is a scalar).

PUBLIC double Matrix_vMv_nocheck(const Matrix *self, const Matrix *vector)
{
	double *array = (double *)memory(CALLOC, self->rows, sizeof(double));
	
	for(size_t col = self->cols; col--;)
	{
		for(size_t row = self->rows; row--;)
		{
			*(array + col) += Matrix_get_value_nocheck(vector, row, 0) * Matrix_get_value_nocheck(self, row, col);
		}
	}
	
	double result = 0.0;
	for(size_t i = self->rows; i--;) result += *(array + i) * Matrix_get_value_nocheck(vector, i, 0);
	
	free(array);
	
	return result;
}



/// @brief Transpose matrix
///
/// Public method for transposing the specified matrix. The transposed
/// matrix will be returned. The user will be responsible for calling
///  the destructor on the transposed matrix once it is no longer needed.
///
/// @param self  Object self-reference.
///
/// @return Transpose of the input matrix.

PUBLIC Matrix *Matrix_transpose(const Matrix *self)
{
	// Sanity checks
	check_null(self);
	
	// Create result matrix with rows and cols swapped
	Matrix *result = Matrix_new(self->cols, self->rows);
	
	// Copy values across
	for(size_t i = 0; i < self->rows; ++i)
	{
		for(size_t j = 0; j < self->cols; ++j)
		{
			result->values[Matrix_get_index(result, j, i)] = self->values[Matrix_get_index(self, i, j)];
		}
	}
	
	return result;
}



/// @brief Invert matrix
///
/// Public method for inverting the specified matrix. The method
/// makes use of the Gauss-Jordan elimination algorithm for this
/// purpose, unless the matrix size is <= 3, in which case the
/// solution is calculated analytically. The inverted matrix will
/// be returned. If the matrix is not invertible, a `NULL` pointer
/// will instead be returned. The user is responsible for calling
/// the destructor once the matrix is no longer needed.
///
/// @param self  Object self-reference.
///
/// @return Inverse of the input matrix, or NULL if not invertible.
///
/// @note The Gauss-Jordan elimination algorithm can be numerically
/// unstable, and integer numbers may get represented as non-integer
/// values.

PUBLIC Matrix *Matrix_invert(const Matrix *self)
{
	// Sanity checks
	check_null(self);
	ensure(self->rows == self->cols, ERR_USER_INPUT, "Cannot invert non-square matrix.");
	
	const size_t size = self->rows;
	
	
	// Special case: matrix size up to 3 x 3 --> use analytic solution
	// ---------------------------------------------------------------
	
	if(size <= 3)
	{
		// Calculate and check determinant
		const double det = Matrix_det(self, 1.0);
		if(det == 0.0)
		{	
			// TODO(austin): Why a warning here and not an ensure?
			warning("Matrix is not invertible.");
			return NULL;
		}
		
		// Create empty inverse matrix
		Matrix *result = Matrix_new(size, size);
		
		// Set inverse matrix entries
		if(size == 1)
		{
			// 1 x 1 matrix
			Matrix_set_value(result, 0, 0, 1.0 / det);
		}
		else if(size == 2)
		{
			// 2 x 2 matrix
			Matrix_set_value(result, 0, 0,  Matrix_get_value(self, 1, 1) / det);
			Matrix_set_value(result, 0, 1, -Matrix_get_value(self, 0, 1) / det);
			Matrix_set_value(result, 1, 0, -Matrix_get_value(self, 1, 0) / det);
			Matrix_set_value(result, 1, 1,  Matrix_get_value(self, 0, 0) / det);
		}
		else if(size == 3)
		{
			// 3 x 3 matrix
			const double a = Matrix_get_value(self, 0, 0);
			const double b = Matrix_get_value(self, 0, 1);
			const double c = Matrix_get_value(self, 0, 2);
			const double d = Matrix_get_value(self, 1, 0);
			const double e = Matrix_get_value(self, 1, 1);
			const double f = Matrix_get_value(self, 1, 2);
			const double g = Matrix_get_value(self, 2, 0);
			const double h = Matrix_get_value(self, 2, 1);
			const double i = Matrix_get_value(self, 2, 2);
			
			Matrix_set_value(result, 0, 0,  (e * i - f * h) / det);
			Matrix_set_value(result, 0, 1,  (c * h - b * i) / det);
			Matrix_set_value(result, 0, 2,  (b * f - c * e) / det);
			Matrix_set_value(result, 1, 0,  (f * g - d * i) / det);
			Matrix_set_value(result, 1, 1,  (a * i - c * g) / det);
			Matrix_set_value(result, 1, 2,  (c * d - a * f) / det);
			Matrix_set_value(result, 2, 0,  (d * h - e * g) / det);
			Matrix_set_value(result, 2, 1,  (b * g - a * h) / det);
			Matrix_set_value(result, 2, 2,  (a * e - b * d) / det);
		}
		
		return result;
	}
	
	
	// General case: N x N matrix --> use Gauss-Jordan elimination
	// -----------------------------------------------------------
	
	// Create initial left and right square matrix
	Matrix *L = Matrix_copy(self);
	Matrix *R = Matrix_identity(size);
	
	// For each column (i) in the matrix
	for(size_t i = 0; i < size; ++i)
	{
		// Find maximum pivot
		double pivot_max = fabs(L->values[Matrix_get_index(L, i, i)]);
		size_t pivot_max_row = i;
		
		for(size_t j = i + 1; j < size; ++j)
		{
			double value = Matrix_get_value(L, j, i);
			if(pivot_max < fabs(value))
			{
				pivot_max = fabs(value);
				pivot_max_row = j;
			}
		}
		
		// Ensure that the maximum pivot is > 0
		// as otherwise the matrix is not invertible
		if(pivot_max == 0.0)
		{
			warning("Matrix is not invertible.");
			Matrix_delete(L);
			Matrix_delete(R);
			return NULL;
		}
		
		// Swap rows if necessary
		if(pivot_max_row != i)
		{
			Matrix_swap_rows(L, i, pivot_max_row);
			Matrix_swap_rows(R, i, pivot_max_row);
		}
		
		// Extract diagonal element (i, i) as pivot
		double pivot = L->values[Matrix_get_index(L, i, i)];
		
		// Divide pivot row (i) by pivot to get a diagonal of 1
		Matrix_mul_row(L, i, 1.0 / pivot);
		Matrix_mul_row(R, i, 1.0 / pivot);
		
		// Subtract multiple of pivot row (i) from every other row (j) to get non-diagonals of 0
		for(size_t j = 0; j < size; ++j)
		{
			if(j != i)
			{
				const double factor = -Matrix_get_value(L, j, i);
				Matrix_add_row(L, j, i, factor);
				Matrix_add_row(R, j, i, factor);
			}
		}
	}
	
	// Clean up
	Matrix_delete(L);
	
	return R;
}



/// @brief Print matrix
///
/// Public method for printing the matrix to the standard output.
/// The width parameter specifies how many characters each printed
/// column is wide. The number of decimals after the decimal point
/// is specified with the decimals argument.
///
/// @param self      Object self-reference.
/// @param width     Width of each column in characters.
/// @param decimals  Number of decimals after the decimal point.

PUBLIC void Matrix_print(const Matrix *self, const unsigned int width, const unsigned int decimals)
{
	// Sanity checks
	check_null(self);
	
	for(size_t row = 0; row < self->rows; ++row)
	{
		for(size_t col = 0; col < self->cols; ++col) printf("%*.*f ", width, decimals, Matrix_get_value(self, row, col));
		printf("\n");
	}
	
	return;
}



/// @brief Calculate determinant
///
/// Public method for returning the determinant of the specified
/// matrix. The determinant for matrices of size <= 3 is calculated
/// analytically, while the determinant of larger matrices will be
/// calculated numerically.
///
/// The purpose of the scale factor is to be able to make use of
/// the fact that det(f * M) = f^n * det(M). By setting the scale
/// factor to a value != 1, one can efficiently calculate the product
/// f^n * det(M), which is needed in the calculation of the PDF of
/// the multivariate normal distribution below.
///
/// @param self          Object self-reference.
/// @param scale_factor  A scale factor by which each element of the
///                      matrix is multiplied before calculating the
///                      determinant. Set to 1 for no scaling.
///
/// @return Determinant of specified matrix.

PUBLIC double Matrix_det(const Matrix *self, const double scale_factor)
{
	// Sanity checks
	check_null(self);
	ensure(self->rows == self->cols, ERR_USER_INPUT, "Cannot calculate determinant of non-square matrix.");
	
	if(self->rows == 1)
	{
		return scale_factor * *(self->values);
	}
	
	if(self->rows == 2)
	{
		return scale_factor * scale_factor * (Matrix_get_value(self, 0, 0) * Matrix_get_value(self, 1, 1) - Matrix_get_value(self, 0, 1) * Matrix_get_value(self, 1, 0));
	}
	
	if(self->rows == 3)
	{
		return scale_factor * scale_factor * scale_factor
			*  (Matrix_get_value(self, 0, 0) * Matrix_get_value(self, 1, 1) * Matrix_get_value(self, 2, 2)
			+   Matrix_get_value(self, 0, 1) * Matrix_get_value(self, 1, 2) * Matrix_get_value(self, 2, 0)
			+   Matrix_get_value(self, 0, 2) * Matrix_get_value(self, 1, 0) * Matrix_get_value(self, 2, 1)
			-   Matrix_get_value(self, 0, 2) * Matrix_get_value(self, 1, 1) * Matrix_get_value(self, 2, 0)
			-   Matrix_get_value(self, 0, 1) * Matrix_get_value(self, 1, 0) * Matrix_get_value(self, 2, 2)
			-   Matrix_get_value(self, 0, 0) * Matrix_get_value(self, 1, 2) * Matrix_get_value(self, 2, 1));
	}
	
	//warning("Determinant calculation for N > 3 not yet implemented.");
	//return 0.0;
	
	// General case for N x N matrix
	// ALERT: THIS ALGORITHM HAS NOT BEEN TESTED YET!
	double det = 1.0;
	Matrix *matrix = Matrix_copy(self);  // Create a copy to work on
	
	for(size_t i = 0; i < matrix->rows; ++i)
	{
		double pivot = Matrix_get_value(matrix, i, i);
		size_t pivot_row = i;
		
		for(size_t row = i + 1; row < matrix->rows; ++row)
		{
			if(fabs(Matrix_get_value(matrix, row, i)) > fabs(pivot))
			{
				pivot = Matrix_get_value(matrix, row, i);
				pivot_row = row;
			}
		}
		
		if(pivot == 0.0) return 0.0;
		
		if(pivot_row != i)
		{
			Matrix_swap_rows(matrix, i, pivot_row);
			det *= -1.0;
		}
		
		det *= pivot;
		
		for(size_t row = i + 1; row < matrix->rows; ++row)
		{
			for(size_t col = i + 1; col < matrix->cols; ++col)
			{
				Matrix_add_value(matrix, row, col, -Matrix_get_value(matrix, row, i) * Matrix_get_value(matrix, i, col) / pivot);
			}
		}
	}
	
	Matrix_delete(matrix);
	
	return det;
}



/// @brief Probability density of a multivariate normal distribution
///
/// Public method for calculating the probability density of a
/// multivariate normal distribution described by the specified
/// covariance matrix at the location specified by `vector`
/// (relative to the mean). For this to work, the covariance
/// matrix must be square and invertible, and the vector must be
/// in column form with a size equal to that of the covariance
/// matrix. The vector entries must be relative to the centroid
/// (x - x_cen) such that 0 corresponds to the peak of the PDF.
/// The PDF will be correctly normalised such that the integral
/// over the entire n-dimensional space is 1 if the correct scale
/// factor is supplied. Also see Matrix_prob_dens_nocheck() for a
/// version without sanity checks for better performance.
///
/// @param covar_inv  Inverse of the covariance matrix. Can be
///                   calculated with Matrix_invert().
/// @param vector     Coordinate (or parameter) vector relative to
///                   the mean for which the probability density is
///                   to be returned.
/// @param scal_fact  Scale factor of 1 divided by the square root of
///                   the determinant of 2 pi times the covariance
///                   matrix, 1 / SQRT(|2 pi COV|).
///
/// @return Probability density at the position of vector.
///
/// @note det(f * M) = f^n * det(M), hence the scale factor is
/// such that 1 / SQRT(|2 pi COV|) = 1 / SQRT((2 pi)^n |COV|),
/// which is the correct normalisation factor of the multivariate
/// normal distribution in n dimensions.

PUBLIC double Matrix_prob_dens(const Matrix *covar_inv, const Matrix *vector, const double scal_fact)
{
	// Sanity checks
	check_null(covar_inv);
	check_null(vector);
	ensure(covar_inv->rows == covar_inv->cols, ERR_USER_INPUT, "Covariance matrix must be square.");
	ensure(covar_inv->rows == vector->rows && vector->cols == 1, ERR_USER_INPUT, "Vector size does not match covariance matrix size.");
	
	// Return PDF = exp(-0.5 v^T C^-1 v) / SQRT((2 pi)^n |C|) of multivariate normal distribution
	return scal_fact * exp(-0.5 * Matrix_vMv_nocheck(covar_inv, vector));
}



/// @brief Probability density of a multivariate normal distribution (no sanity checks)
///
/// Public method for calculating the probability density of a
/// multivariate normal distribution described by the specified
/// covariance matrix at the location specified by `vector`
/// (relative to the mean). For this to work, the covariance
/// matrix must be square and invertible, and the vector must be
/// in column form with a size equal to that of the covariance
/// matrix. The vector entries must be relative to the centroid
/// (x - x_cen) such that 0 corresponds to the peak of the PDF.
/// The PDF will be correctly normalised such that the integral
/// over the entire n-dimensional space is 1 if the correct scale
/// factor is supplied. Identical to Matrix_prob_dens(), but
/// without sanity checks for better performance.
///
/// @param covar_inv  Inverse of the covariance matrix. Can be
///                   calculated with Matrix_invert().
/// @param vector     Coordinate (or parameter) vector relative to
///                   the mean for which the probability density is
///                   to be returned.
/// @param scal_fact  Scale factor of 1 divided by the square root of
///                   the determinant of 2 pi times the covariance
///                   matrix, 1 / SQRT(|2 pi COV|).
///
/// @return Probability density at the position of vector.
///
/// @note det(f * M) = f^n * det(M), hence the scale factor is
/// such that 1 / SQRT(|2 pi COV|) = 1 / SQRT((2 pi)^n |COV|),
/// which is the correct normalisation factor of the multivariate
/// normal distribution in n dimensions.

PUBLIC double Matrix_prob_dens_nocheck(const Matrix *covar_inv, const Matrix *vector, const double scal_fact)
{
	// Return PDF = exp(-0.5 v^T C^-1 v) / SQRT((2 pi)^n |C|) of multivariate normal distribution
	return scal_fact * exp(-0.5 * Matrix_vMv_nocheck(covar_inv, vector));
}



/// @brief Create error ellipse from covariance matrix
///
/// Public method for determining the radii and position angle of
/// the error ellipse corresponding to the specified covariance
/// matrix. The radii correspond to the standard deviation (sigma).
/// The results will be written into the parameters `radius_maj`,
/// `radius_min` and `pa`.
///
/// @param covar       Covariance matrix.
/// @param par1        First parameter to use.
/// @param par2        Second parameter to use.
/// @param radius_maj  Radius of the ellipse along major axis.
/// @param radius_min  Radius of the ellipse along minor axis.
/// @param pa          Position angle of the ellipse.
///
/// @note The covariance matrix must be a square matrix for this
/// to make sense. The parameters `par1` and `par2` are used to
/// select two of of the N parameters of the covariance matrix, and
/// the ellipse parameters will then be determined for this specific
/// 2-D projection.

PUBLIC void Matrix_err_ellipse(const Matrix *covar, const size_t par1, const size_t par2, double *radius_maj, double *radius_min, double *pa)
{
	// Sanity checks
	check_null(covar);
	ensure(covar->rows == covar->cols, ERR_USER_INPUT, "Error ellipse can only be derived for square matrix.");
	ensure(covar->rows >= 2, ERR_USER_INPUT, "Covariance matrix must have size 2 or greater.");
	ensure(par1 < covar->rows && par2 < covar->rows, ERR_INDEX_RANGE, "Covariance matrix row index out of range.");
	ensure(par1 != par2, ERR_USER_INPUT, "Please specify two different rows.");
	
	// Extract values from covariance matrix
	const double v1 = Matrix_get_value(covar, par1, par1);
	const double v2 = Matrix_get_value(covar, par2, par2);
	const double c  = Matrix_get_value(covar, par2, par1);
	
	// Some settings
	//const double scale_factor = -2.0 * log(1.0 - confidence);  // NOTE: Can be used to draw confidence level rather than sigma.
	const double scale_factor = 1.0;
	const double tmp = sqrt((v1 - v2) * (v1 - v2) + 4 * c * c);
	
	// Define eigenvalues
	Matrix *eigenvalues = Matrix_new(2, 1);
	Matrix_set_value(eigenvalues, 0, 0, 0.5 * (v1 + v2 - tmp));
	Matrix_set_value(eigenvalues, 1, 0, 0.5 * (v1 + v2 + tmp));
	
	// Define eigenvectors (not yet normalised)
	Matrix *eigenvectors = Matrix_new(2, 2);
	Matrix_set_value(eigenvectors, 0, 0, 1.0);
	Matrix_set_value(eigenvectors, 1, 0, (Matrix_get_value(eigenvalues, 0, 0) - v1) / c);
	Matrix_set_value(eigenvectors, 0, 1, 1.0);
	Matrix_set_value(eigenvectors, 1, 1, (Matrix_get_value(eigenvalues, 1, 0) - v1) / c);
	
	// Normalise eigenvectors to length 1
	const double norm1 = 1.0 / sqrt(Matrix_get_value(eigenvectors, 1, 0) * Matrix_get_value(eigenvectors, 1, 0) + 1.0);
	const double norm2 = 1.0 / sqrt(Matrix_get_value(eigenvectors, 1, 1) * Matrix_get_value(eigenvectors, 1, 1) + 1.0);
	Matrix_mul_value(eigenvectors, 0, 0, norm1);
	Matrix_mul_value(eigenvectors, 1, 0, norm1);
	Matrix_mul_value(eigenvectors, 0, 1, norm2);
	Matrix_mul_value(eigenvectors, 1, 1, norm2);
	
	// Check which eigenvalue is largest
	const size_t idx1 = Matrix_get_value(eigenvalues, 0, 0) > Matrix_get_value(eigenvalues, 1, 0) ? 0 : 1;
	const size_t idx2 = 1 - idx1;
	
	// Calculate ellipse sizes
	const double x1 = Matrix_get_value(eigenvectors, 0, idx1) * sqrt(Matrix_get_value(eigenvalues, idx1, 0) * scale_factor);
	const double y1 = Matrix_get_value(eigenvectors, 1, idx1) * sqrt(Matrix_get_value(eigenvalues, idx1, 0) * scale_factor);
	
	const double x2 = Matrix_get_value(eigenvectors, 0, idx2) * sqrt(Matrix_get_value(eigenvalues, idx2, 0) * scale_factor);
	const double y2 = Matrix_get_value(eigenvectors, 1, idx2) * sqrt(Matrix_get_value(eigenvalues, idx2, 0) * scale_factor);
	
	*radius_maj = sqrt(x1 * x1 + y1 * y1);
	*radius_min = sqrt(x2 * x2 + y2 * y2);
	*pa = atan2(y1, x1);
	
	// Clean up again
	Matrix_delete(eigenvalues);
	Matrix_delete(eigenvectors);
	
	return;
}



/// @brief Get array index from row and column
///
/// Private method for returning the array index corresponding to
/// the specified row and column number. This can be used to
/// directly access the data array from within public and private
/// methods without having to use the Matrix_get_value() and
/// Matrix_set_value() methods.
///
/// @param self  Object self-reference.
/// @param row   Row number.
/// @param col   Column number.
///
/// @return Array index corresponding to the row and column number.

PRIVATE inline size_t Matrix_get_index(const Matrix *self, const size_t row, const size_t col)
{
	return row + self->rows * col;
}



/// @brief Swap two matrix rows
///
/// Private method for swapping the two rows specified by `row1` and
/// `row2` of the matrix. This method is needed for the Gauss-Jordan
/// elimination algorithm.
///
/// @param self  Object self-reference.
/// @param row1  First row to be swapped.
/// @param row2  Second row to be swapped.

PRIVATE void Matrix_swap_rows(Matrix *self, const size_t row1, const size_t row2)
{
	for(size_t i = self->cols; i--;) swap(self->values + Matrix_get_index(self, row1, i), self->values + Matrix_get_index(self, row2, i));
	return;
}



/// @brief Add `factor * row2` to `row1`
///
/// Private method for multiplying `row2` by `factor` and adding
/// the result to `row1` of the specified matrix. This method is
/// needed for the Gauss-Jordan elimination algorithm.
///
/// @param self    Object self-reference.
/// @param row1    Row to be manipulated.
/// @param row2    Row to be scaled and added to row1.
/// @param factor  Factor to multiply row2 by.

PRIVATE void Matrix_add_row(Matrix *self, const size_t row1, const size_t row2, const double factor)
{
	for(size_t i = self->cols; i--;) self->values[Matrix_get_index(self, row1, i)] += factor * self->values[Matrix_get_index(self, row2, i)];
	return;
}



/// @brief Multiply row by factor
///
/// Private method for multiplying the specified row of the
/// matrix by the specified factor. This method is needed for
/// the Gauss-Jordan elimination algorithm.
///
/// @param self    Object self-reference.
/// @param row     Row to be multipled.
/// @param factor  Factor to multiply row by.

PRIVATE void Matrix_mul_row(Matrix *self, const size_t row, const double factor)
{
	for(size_t i = self->cols; i--;) self->values[Matrix_get_index(self, row, i)] *= factor;
	return;
}
