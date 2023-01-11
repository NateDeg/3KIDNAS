// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Catalog.h) - Source Finding Application                 //
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

/// @file   Catalog.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for storing source catalogues (header).


#ifndef CATALOG_H
#define CATALOG_H

#include <stdint.h>
#include "common.h"
#include "Source.h"
#include "Parameter.h"

#define CATALOG_COLUMN_WIDTH 14  ///< Defines the width of each column in the plain-text SoFiA source catalogue.

typedef enum {CATALOG_FORMAT_ASCII, CATALOG_FORMAT_XML, CATALOG_FORMAT_SQL} file_format;


// ----------------------------------------------------------------- //
// Class 'Catalog'                                                   //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a structure for storing   //
// and handling source catalogues as a simple list of objects of     //
// class 'Source'.                                                   //
// ----------------------------------------------------------------- //

typedef CLASS Catalog Catalog;

// Constructor and destructor
PUBLIC  Catalog *Catalog_new           (void);
PUBLIC  void     Catalog_delete        (Catalog *self);

// Public methods
PUBLIC  void     Catalog_add_source    (Catalog *self, Source *src);
PUBLIC  Source  *Catalog_get_source    (const Catalog *self, const size_t index);
PUBLIC  size_t   Catalog_get_index     (const Catalog *self, const Source *src);
PUBLIC  bool     Catalog_source_exists (const Catalog *self, const Source *src, size_t *index);

PUBLIC  size_t   Catalog_get_size      (const Catalog *self);

PUBLIC  void     Catalog_save          (const Catalog *self, const char *filename, const file_format format, const bool overwrite, const Parameter *par);

// Private methods
PRIVATE void     Catalog_append_memory (Catalog *self);

#endif
