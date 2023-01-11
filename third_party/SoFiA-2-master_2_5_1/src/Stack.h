// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Stack.h) - Source Finding Application                   //
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

/// @file   Stack.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class implementing a simple LIFO stack for the recursive linker algorithm (header).


#ifndef STACK_H
#define STACK_H

#include "common.h"


// ----------------------------------------------------------------- //
// Class 'Stack'                                                     //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a simple LIFO stack for a //
// recursive pixel linking algorithm. The stack is capable of stor-  //
// ing the indices of the pixels identified as part of a source so   //
// their neighbours can be checked recursively.                      //
// ----------------------------------------------------------------- //

typedef CLASS Stack Stack;

// Constructor and destructor
PUBLIC Stack        *Stack_new      (void);
PUBLIC void          Stack_delete   (Stack *self);

// Public methods
PUBLIC void          Stack_push     (Stack *self, const size_t value);
PUBLIC size_t        Stack_pop      (Stack *self);
PUBLIC size_t        Stack_get_size (const Stack *self);

#endif
