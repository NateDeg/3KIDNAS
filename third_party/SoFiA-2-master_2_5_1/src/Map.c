// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Map.c) - Source Finding Application                     //
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

/// @file   Map.c
/// @author Tobias Westmeier
/// @date   23/11/2021
/// @brief  Container class for storing key-value pairs of type `size_t`.


#include <stdlib.h>
#include "Map.h"



/// @brief Container class for storing key-value pairs of type `size_t`.
///
/// The purpose of this class is to provide a structure for storing
/// and updating source parameters handled by the linker implemented
/// in the class DataCube.

CLASS Map
{
	size_t  size;    ///< Number of key-value pairs stored.
	size_t *keys;    ///< Pointer to array of keys.
	size_t *values;  ///< Pointer to array of values.
};



/// Standard constructor
///
/// Standard constructor. Will create a new, empty Map object and
/// return a pointer to the newly created object. Note that the
/// destructor will need to be called explicitly once the object is
/// no longer required to release any memory allocated during the
/// lifetime of the object.
///
/// @return Pointer to newly created Map object.

PUBLIC Map *Map_new(void)
{
	Map *self = (Map *)memory(MALLOC, 1, sizeof(Map));
	
	self->size = 0;
	self->keys = NULL;
	self->values = NULL;
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Map_delete(Map *self)
{
	if(self != NULL)
	{
		free(self->keys);
		free(self->values);
		free(self);
	}
	
	return;
}



/// @brief Push new key-value pair onto map
///
/// Public method for pushing a new key-value pair onto the specified
/// map. Note that there will be no check as to whether the key
/// already exists, and it is therefore possible to create more than
/// one entry with the same key.
///
/// @param self   Object self-reference.
/// @param key    Key to be created.
/// @param value  Value to be added.

PUBLIC void Map_push(Map *self, const size_t key, const size_t value)
{
	// Sanity checks
	check_null(self);
	
	++self->size;
	
	self->keys   = (size_t *)memory_realloc(self->keys, self->size, sizeof(size_t));
	self->values = (size_t *)memory_realloc(self->values, self->size, sizeof(size_t));
	
	self->keys[self->size - 1] = key;
	self->values[self->size - 1] = value;
	
	return;
}



/// @brief Retrieve value by key
///
/// Public method for retrieving the value associated with the
/// specified key. If the same key exists more than once, the
/// last occurrence will be retrieved. The process will be
/// terminated of the specified map is empty or the key is not
/// found.
///
/// @param self  Object self-reference.
/// @param key   Key to be retrieved.
///
/// @return Value belonging to `key`.

PUBLIC size_t Map_get_value(const Map *self, const size_t key)
{
	// Sanity checks
	check_null(self);
	
	// Search for key and return value
	for(size_t i = self->size; i--;) if(self->keys[i] == key) return self->values[i];
	
	// Key not found
	ensure(false, ERR_USER_INPUT, "Key \'%zu\' not found in map.", key);
	
	// While this return statement is never executed, it is required by the compiler:
	return 0;
}



/// @brief Return size of map
///
/// Public method for returning the current size of the specified
/// map, i.e. the number of key-value pairs currently stored. If
/// the same key exists multiple times, it will also be counted
/// repeatedly.
///
/// @param self  Object self-reference.
///
/// @return Current size of map, i.e. number of entries.

PUBLIC size_t Map_get_size(const Map *self)
{
	check_null(self);
	return self->size;
}



/// @brief Check if key exists
///
/// Public method for checking if the specified key exists. The
/// method will return `true` if the key is found and `false`
/// otherwise.
///
/// @param self  Object self-reference.
/// @param key   Key to be checked.
///
/// @return Returns `true` if key exists, `false` otherwise.

PUBLIC bool Map_key_exists(const Map *self, const size_t key)
{
	// Sanity checks
	check_null(self);
	
	// Search for key
	for(size_t i = 0; i < self->size; ++i) if(self->keys[i] == key) return true;
	return false;
}
