// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Header.h) - Source Finding Application                  //
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

/// @file   Header.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for storing and managing FITS file headers (header).


#ifndef HEADER_H
#define HEADER_H

#include "common.h"
#include "String.h"

#define FITS_HEADER_BLOCK_SIZE   2880  ///< Size of a single FITS header block (in bytes).
#define FITS_HEADER_LINE_SIZE      80  ///< Size of a single FITS header line (in bytes).
#define FITS_HEADER_LINES          36  ///< Number of lines in a single FITS header block.
#define FITS_HEADER_KEYWORD_SIZE    8  ///< Maximum size of a FITS header keyword (in bytes).
#define FITS_HEADER_KEY_SIZE       10  ///< Size of a FITS header key, including '=' assignment (in bytes).
#define FITS_HEADER_VALUE_SIZE     70  ///< Maximum size of a FITS header value (in bytes).
#define FITS_HEADER_FIXED_WIDTH    20  ///< Size of a FITS header value of fixed width (in bytes).


// ----------------------------------------------------------------- //
// Class 'Header'                                                    //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a structure for storing   //
// header information from a FITS data file and implement ways of    //
// reading and manipulating individual header entries.               //
// ----------------------------------------------------------------- //

typedef CLASS Header Header;

// Constructor and destructor
PUBLIC  Header     *Header_new        (const char *header, const size_t size, const bool verbosity);
PUBLIC  Header     *Header_copy       (const Header *source);
PUBLIC  Header     *Header_blank      (const bool verbosity);
PUBLIC  void        Header_delete     (Header *self);

// Public methods
PUBLIC  const char *Header_get        (const Header *self);
PUBLIC  size_t      Header_get_size   (const Header *self);

// Extract header entries
PUBLIC  long int    Header_get_int    (const Header *self, const char *key);
PUBLIC  double      Header_get_flt    (const Header *self, const char *key);
PUBLIC  bool        Header_get_bool   (const Header *self, const char *key);
PUBLIC  int         Header_get_str    (const Header *self, const char *key, char *value);
PUBLIC  String     *Header_get_string (const Header *self, const char *key);

// Manipulate header entries
PUBLIC  int         Header_set_int    (Header *self, const char *key, const long int value);
PUBLIC  int         Header_set_flt    (Header *self, const char *key, const double value);
PUBLIC  int         Header_set_bool   (Header *self, const char *key, const bool value);
PUBLIC  int         Header_set_str    (Header *self, const char *key, const char *value);
PUBLIC  size_t      Header_comment    (Header *self, const char *value, const bool history);

// Miscellaneous header operations
PUBLIC  size_t      Header_check      (const Header *self, const char *key);
PUBLIC  bool        Header_compare    (const Header *self, const char *key, const char *value, const size_t n);
PUBLIC  int         Header_remove     (Header *self, const char *key);
PUBLIC  void        Header_copy_wcs   (const Header *source, Header *target);
PUBLIC  void        Header_copy_misc  (const Header *source, Header *target, const bool copy_bunit, const bool copy_beam);
PUBLIC  void        Header_adjust_wcs_to_subregion(Header *self, const size_t x_min, const size_t x_max, const size_t y_min, const size_t y_max, const size_t z_min, const size_t z_max);

// Private methods
PRIVATE int         Header_get_raw    (const Header *self, const char *key, char *buffer);
PRIVATE int         Header_set_raw    (Header *self, const char *key, const char *buffer);

#endif
