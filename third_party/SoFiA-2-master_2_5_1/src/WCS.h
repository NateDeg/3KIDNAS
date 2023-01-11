// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (WCS.h) - Source Finding Application                     //
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

/// @file   WCS.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for handling and converting world coordinates (header).


#ifndef WCS_H
#define WCS_H

#include "common.h"


// ----------------------------------------------------------------- //
// Class 'WCS'                                                       //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a support for World Coor- //
// dinate System (WCS) conversions in the form of a wrapper around   //
// the wcslib package. The class provides methods for setting up WCS //
// information from a FITS file header and converting between pixel  //
// coordinate (x, y, z) and world coordinates (lon, lat, spec) in    //
// three dimensions.                                                 //
// ----------------------------------------------------------------- //

typedef CLASS WCS WCS;

// Constructor and destructor
PUBLIC  WCS  *WCS_new            (const char *header, const int n_keys, const int n_axes, const int *dim_axes);
PUBLIC  void  WCS_delete         (WCS *self);

// Public methods
PUBLIC  bool  WCS_is_valid       (const WCS *self);
PUBLIC  void  WCS_convertToWorld (const WCS *self, const double x, const double y, const double z, double *longitude, double *latitude, double *spectral);
PUBLIC  void  WCS_convertToPixel (const WCS *self, const double longitude, const double latitude, const double spectral, double *x, double *y, double *z);

// Private methods
PRIVATE void  WCS_setup          (WCS *self, const char *header, const int n_keys, const int n_axes, const int *dim_axes);

#endif
