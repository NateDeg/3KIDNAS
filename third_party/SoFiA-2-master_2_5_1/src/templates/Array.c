// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Array_SFX.c) - Source Finding Application               //
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

/// @file   Array_SFX.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Container class template for storing a dynamic array of values of type `DATA_T`.


// WARNING: This is a template that needs to be instantiated before use.
//          Do not edit template instances, as they are auto-generated
//          and will be overwritten during instantiation!


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Array_SFX.h"



/// @brief Container class template for storing a dynamic array of values of type `DATA_T`.
///
/// The purpose of this class is to provide a convenient way to store
/// multiple values of a specific type in an array-like structure. A
/// new array can either be of a given size and empty (using the
/// standard constructor) or provided with a list of comma-separated
/// values that will be stored in the array and used to determine its
/// size (using the alternative constructor).

CLASS Array_SFX
{
	size_t size;     ///< Size of the array.
	DATA_T *values;  ///< Pointer to array of values.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new Array_SFX object of
/// given size and type and return a pointer to the newly created
/// object. Sufficient memory will be allocated to store the array
/// values of the specified type. Note that the destructor will
/// need to be called explicitly once the object is no longer
/// required to release any memory allocated during the lifetime
/// of the object. Note that the array will be initialised to 0.
///
/// @param size  Size of the array to be created.
///
/// @return Pointer to newly created Array_SFX object.

PUBLIC Array_SFX *Array_SFX_new(const size_t size)
{
	Array_SFX *self = (Array_SFX *)memory(MALLOC, 1, sizeof(Array_SFX));
	
	self->size = size;
	
	if(self->size) self->values = (DATA_T *)memory(CALLOC, size, sizeof(DATA_T));
	else self->values = NULL;
	
	return self;
}



/// @brief Alternative constructor
///
/// Alternative constructor. Will create a new Array_SFX object,
/// the size of which is determined by the number of comma-separated
/// values specified in `string`. A pointer to the newly created
/// object will be returned. Sufficient memory will be allocated
/// to store the array values. Note that the destructor will need
/// to be called explicitly once the object is no longer required to
/// release any memory allocated during the lifetime of the object.
///
/// @param string  String containing the values to be stored in the
///                array, separated by commas.
///
/// @return Pointer to newly created Array_SFX object.

PUBLIC Array_SFX *Array_SFX_new_str(const char *string)
{
	// Sanity checks
	check_null(string);
	
	// Return empty array if string is empty
	if(!strlen(string)) return Array_SFX_new(0);
	
	// Create a copy of the string
	char *copy = (char *)memory(MALLOC, strlen(string) + 1, sizeof(char));
	strcpy(copy, string);
	
	// Count number of commas
	size_t size = 1;
	size_t i = strlen(copy);
	while(i--) if(copy[i] == ',') ++size;
	
	// Create array of given size
	Array_SFX *self = Array_SFX_new(size);
	
	// Fill array with values
	char *token = strtok(copy, ",");
	ensure(token != NULL, ERR_USER_INPUT, "Failed to parse string as array.");
	
	self->values[0] = (DATA_T)strtod(token, NULL);
	
	for(i = 1; i < size; ++i)
	{
		token = strtok(NULL, ",");
		ensure(token != NULL, ERR_USER_INPUT, "Failed to parse string as array.");
		
		self->values[i] = (DATA_T)strtod(token, NULL);
	}
	
	// Delete string copy again
	free(copy);
	
	return self;
}



/// @brief Copy constructor
///
/// Copy constructor. Will create a new array of the same size as
/// the source array and then copy all elements from source. A
/// pointer to the newly created array copy will be returned. Note
/// that the destructor will need to be called explicitly once the
/// object is no longer required to release any memory allocated
/// during the lifetime the object.
///
/// @param source  Array to be copied.
///
/// @return Pointer to newly created copy of array object.

PUBLIC Array_SFX *Array_SFX_copy(const Array_SFX *source)
{
	// Sanity checks
	check_null(source);
	
	// Create new array of same size as source
	Array_SFX *self = Array_SFX_new(source->size);
	
	// Copy all elements
	for(size_t i = self->size; i--;) self->values[i] = source->values[i];
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Array_SFX_delete(Array_SFX *self)
{
	if(self != NULL) free(self->values);
	free(self);
	return;
}



/// @brief Get size of array
///
/// Public method for returning the size of the specified array.
///
/// @param self  Object self-reference.
///
/// @return Size of the array.

PUBLIC size_t Array_SFX_get_size(const Array_SFX *self)
{
	check_null(self);
	return self->size;
}



/// @brief Get pointer to data array
///
/// Public method for returning a pointer to the first element of
/// the array.
///
/// @param self  Object self-reference.
///
/// @return Pointer to the first element of the array. If the array has
///         size 0, a `NULL` pointer will be returned instead.

PUBLIC const DATA_T *Array_SFX_get_ptr(const Array_SFX *self)
{
	check_null(self);
	return self->size ? self->values : NULL;
}



/// @brief Push new element
///
/// Public method for adding a new element at the end of the array.
///
/// @param self   Object self-reference.
/// @param value  Value to be added.
///
/// @return Pointer to modified array.

PUBLIC Array_SFX *Array_SFX_push(Array_SFX *self, const DATA_T value)
{
	check_null(self);
	self->values = (DATA_T *)memory_realloc(self->values, ++self->size, sizeof(DATA_T));
	self->values[self->size - 1] = value;
	return self;
}



/// @brief Get array element
///
/// Public method for retrieving the array value at the specified
/// index.
///
/// @param self   Object self-reference.
/// @param index  Index of the element to be returned.
///
/// @return Value of the requested element.

PUBLIC DATA_T Array_SFX_get(const Array_SFX *self, const size_t index)
{
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Array index out of range.");
	return self->values[index];
}



/// @brief Set array element
///
/// Public method for setting the value of the array element at the
/// specified index.
///
/// @param self   Object self-reference.
/// @param index  Index of the element to be set.
/// @param value  Value of the element to be set.
///
/// @return Pointer to modified array.

PUBLIC Array_SFX *Array_SFX_set(Array_SFX *self, const size_t index, const DATA_T value)
{
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Array index out of range.");
	self->values[index] = value;
	return self;
}



/// @brief Add value to array element
///
/// Public method for adding the specified value to the array
/// element at the specified index.
///
/// @param self   Object self-reference.
/// @param index  Index of the element to be set.
/// @param value  Value to be added to the element.
///
/// @return Pointer to modified array.

PUBLIC Array_SFX *Array_SFX_add(Array_SFX *self, const size_t index, const DATA_T value)
{
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Array index out of range.");
	self->values[index] += value;
	return self;
}



/// @brief Multiply array element by value
///
/// Public method for multiplying the array element at the
/// specified index by the specified value.
///
/// @param self   Object self-reference.
/// @param index  Index of the element to be set.
/// @param value  Value to be multiplied by.
///
/// @return Pointer to modified array.

PUBLIC Array_SFX *Array_SFX_mul(Array_SFX *self, const size_t index, const DATA_T value)
{
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Array index out of range.");
	self->values[index] *= value;
	return self;
}



/// @brief Concatenate two arrays
///
/// Public method for concatenating two arrays by adding the
/// elements of `source` at the end of `self`.
///
/// @param self    Object self-reference.
/// @param source  Array to be added.
///
/// @return Pointer to modified array.

PUBLIC Array_SFX *Array_SFX_cat(Array_SFX *self, const Array_SFX *source)
{
	check_null(self);
	if(source == NULL) return self;
	
	self->values = (DATA_T *)memory_realloc(self->values, self->size + source->size, sizeof(DATA_T));
	for(size_t i = 0; i < source->size; ++i) self->values[self->size + i] = source->values[i];
	self->size += source->size;
	
	return self;
}



/// @brief Sort array elements
///
/// Public method for sorting the array in ascending order. A pointer
/// to the sorted array will be returned for convenience to allow
/// chaining of methods.
///
/// @param self  Object self-reference.
///
/// @return Pointer to sorted array.

PUBLIC Array_SFX *Array_SFX_sort(Array_SFX *self)
{
	check_null(self);
	if(self->size < 2) return self;
	
	// Use insertion sort algorithm to sort array
	for(size_t i = 1; i < self->size; ++i)
	{
		size_t j = i;
		while(j > 0 && self->values[j - 1] > self->values[j])
		{
			DATA_T tmp = self->values[j];
			self->values[j] = self->values[j - 1];
			self->values[j - 1] = tmp;
			--j;
		}
	}
	
	return self;
}
