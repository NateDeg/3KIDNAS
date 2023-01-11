// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Stack.c) - Source Finding Application                   //
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

/// @file   Stack.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Class implementing a simple LIFO stack for the recursive linker algorithm.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Stack.h"



/// @brief Class implementing a simple LIFO stack for the recursive linker algorithm.
///
/// The purpose of this class is to provide a simple LIFO stack for a
/// recursive pixel linking algorithm. The stack is capable of storing
/// the indices of the pixels identified as part of a source so their
/// neighbours can be checked recursively.

CLASS Stack
{
	size_t  size;  ///< Current size of the stack, i.e. number of elements stored.
	size_t *data;  ///< Pointer to array of currently stored stack elements.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new and empty Stack object.
/// Note that the destructor will need to be called explicitly once
/// the object is no longer required to release its memory again.
///
/// @return Pointer to newly created Stack object.

PUBLIC Stack *Stack_new(void)
{
	Stack *self = (Stack *)memory(MALLOC, 1, sizeof(Stack));
	
	self->size = 0;
	self->data = NULL;
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Stack_delete(Stack *self)
{
	if(self != NULL) free(self->data);
	free(self);
	
	return;
}



/// @brief Push element onto stack
///
/// Public method for pushing a new element onto the stack. The
/// stack will automatically expand its memory to hold the extra
/// item and terminate with a stack overflow error if memory
/// allocation fails.
///
/// @param self   Object self-reference.
/// @param value  Value to be pushed.

PUBLIC void Stack_push(Stack *self, const size_t value)
{
	// Sanity checks
	check_null(self);
	
	self->size += 1;
	self->data = (size_t *)realloc(self->data, self->size * sizeof(size_t));
	ensure(self->data != NULL, ERR_MEM_ALLOC, "Stack overflow error at %.5f GB memory usage.", (double)(self->size * sizeof(size_t)) / GIGABYTE);
	self->data[self->size - 1] = value;
	
	return;
}



/// @brief Pop element from stack
///
/// Public method for popping the last element from the stack. The
/// stack will automatically adjust its memory allocation to the
/// new, smaller size after popping. A stack underflow error will
/// be raised and the process terminated if the method is called on
/// an empty stack.
///
/// @param self  Object self-reference.
///
/// @return Value of the element to be popped.

PUBLIC size_t Stack_pop(Stack *self)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_FAILURE, "Stack underflow error.");
	check_null(self->data);
	
	const size_t value = self->data[self->size - 1];
	self->size -= 1;
	
	if(self->size)
	{
		self->data = (size_t *)realloc(self->data, self->size * sizeof(size_t));
		ensure(self->data != NULL, ERR_MEM_ALLOC, "Stack overflow error at %.5f GB memory usage.", (double)(self->size * sizeof(size_t)) / GIGABYTE);
	}
	else
	{
		free(self->data);
		self->data = NULL;
	}
	
	return value;
}



/// @brief Retrieve current stack size
///
/// Public method for retrieving the current size of the stack,
/// i.e. the number of elements currently stored.
///
/// @param self  Object self-reference.
///
/// @return Size of the stack (i.e. number of elements).

PUBLIC size_t Stack_get_size(const Stack *self)
{
	check_null(self);
	return self->size;
}
