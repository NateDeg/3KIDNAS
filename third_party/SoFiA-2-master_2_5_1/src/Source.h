// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Source.h) - Source Finding Application                  //
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

/// @file   Source.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for storing and handling source parameters (header).


#ifndef SOURCE_H
#define SOURCE_H

#include "common.h"

#define SOURCE_TYPE_INT 0
#define SOURCE_TYPE_FLT 1


// ----------------------------------------------------------------- //
// Class 'Source'                                                    //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a structure for storing   //
// and handling the measured parameters of a source. Parameters are  //
// composed of a name, value, type and unit. Both long integer and   //
// double-precision floating-point values are supported. In addi-    //
// tion, a source can be assigned an identifier in the form of a     //
// string, e.g. a source name.                                       //
// ----------------------------------------------------------------- //

typedef CLASS Source Source;

// Constructor and destructor
PUBLIC  Source       *Source_new                 (const bool verbosity);
PUBLIC  void          Source_delete              (Source *self);

// Public methods
PUBLIC  void          Source_set_identifier      (Source *self, const char *name);
PUBLIC  const char   *Source_get_identifier      (const Source *self);

PUBLIC  size_t        Source_get_num_par         (const Source *self);

PUBLIC  void          Source_add_par_flt         (Source *self, const char *name, const double   value, const char *unit, const char *ucd);
PUBLIC  void          Source_add_par_int         (Source *self, const char *name, const long int value, const char *unit, const char *ucd);
PUBLIC  void          Source_set_par_flt         (Source *self, const char *name, const double   value, const char *unit, const char *ucd);
PUBLIC  void          Source_set_par_int         (Source *self, const char *name, const long int value, const char *unit, const char *ucd);
PUBLIC  double        Source_get_par_flt         (const Source *self, const size_t index);
PUBLIC  long int      Source_get_par_int         (const Source *self, const size_t index);
PUBLIC  double        Source_get_par_by_name_flt (const Source *self, const char *name);
PUBLIC  long int      Source_get_par_by_name_int (const Source *self, const char *name);

PUBLIC  void          Source_offset_xyz          (Source *self, const size_t dx, const size_t dy, const size_t dz);

PUBLIC  bool          Source_par_exists          (const Source *self, const char *name, size_t *index);

PUBLIC  const char   *Source_get_name            (const Source *self, const size_t index);
PUBLIC  const char   *Source_get_unit            (const Source *self, const size_t index);
PUBLIC  unsigned char Source_get_type            (const Source *self, const size_t index);
PUBLIC  const char   *Source_get_ucd             (const Source *self, const size_t index);

// Private methods
PRIVATE void          Source_append_memory       (Source *self);

#endif
