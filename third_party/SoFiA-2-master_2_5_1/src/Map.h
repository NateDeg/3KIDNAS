// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Map.h) - Source Finding Application                     //
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

/// @file   Map.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Container class for storing key-value pairs of type `size_t` (header).


#ifndef MAP_H
#define MAP_H

#include "common.h"


// ----------------------------------------------------------------- //
// Class 'Map'                                                       //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a simple map structure    //
// that allows key-value pairs to be handled. Both key and value are //
// of type size_t (to facilitate handling of source mask labels).    //
// Multiple entries with the same key are allowed.                   //
// ----------------------------------------------------------------- //

typedef CLASS Map Map;

// Constructor and destructor
PUBLIC Map          *Map_new        (void);
PUBLIC void          Map_delete     (Map *self);

// Methods
PUBLIC void          Map_push       (Map *self, const size_t key, const size_t value);
PUBLIC size_t        Map_get_value  (const Map *self, const size_t key);
PUBLIC size_t        Map_get_size   (const Map *self);
PUBLIC bool          Map_key_exists (const Map *self, const size_t key);

#endif
