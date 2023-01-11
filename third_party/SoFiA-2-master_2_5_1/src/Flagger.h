// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Flagger.h) - Source Finding Application                 //
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

/// @file   Flagger.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for storing data flagging information (header). **WORK IN PROGRESS!**


#ifndef FLAGGER_H
#define FLAGGER_H

#include "common.h"

// List of supported flagging region shapes
// NOTE: SHAPE_COUNT must always be the last item and keeps
//       track of the total number of supported shapes.
enum {PIXEL, CHANNEL, REGION, CIRCLE, SHAPE_COUNT};


// ----------------------------------------------------------------- //
// Class 'Flagger'                                                   //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a structure for storing   //
// flagging information. It supports several different shapes of     //
// regions, such as individual pixels or 3D rectangular regions. All //
// parameters are expected to be in units of integer pixels.         //
// ----------------------------------------------------------------- //

typedef CLASS Flagger Flagger;

// Constructor and destructor
PUBLIC Flagger *Flagger_new(void);
PUBLIC void     Flagger_delete(Flagger *self);

// Public methods
PUBLIC size_t   Flagger_size(const Flagger *self);
PUBLIC int      Flagger_npar(const Flagger *self, const int shape);
PUBLIC Flagger *Flagger_add(Flagger *self, const int shape, ...);
PUBLIC void     Flagger_get(const Flagger *self, const size_t index, int *shape, const long int **parameters);

#endif
