// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Path.h) - Source Finding Application                    //
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

/// @file   Path.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Class for storing and handling file paths (header).


#ifndef SoFiA_PATH_H
#define SoFiA_PATH_H

#include "common.h"


// ----------------------------------------------------------------- //
// Class 'Path'                                                      //
// ----------------------------------------------------------------- //
// The purpose of this class is to provide a structure for storing   //
// and handling file paths in a Linux/Unix directory system. Various //
// methods for setting and reading the different components of the   //
// path are available.                                               //
// ----------------------------------------------------------------- //

typedef CLASS Path Path;

// Constructor and destructor
PUBLIC Path       *Path_new                      (void);
PUBLIC void        Path_delete                   (Path *self);

// Public methods
PUBLIC void        Path_set                      (Path *self, const char *path);
PUBLIC void        Path_set_file                 (Path *self, const char *file);
PUBLIC void        Path_set_dir                  (Path *self, const char *dir);
PUBLIC void        Path_set_file_from_template   (Path *self, const char *basename, const char *suffix, const char *mimetype);
PUBLIC void        Path_append_dir_from_template (Path *self, const char *basename, const char *appendix);
PUBLIC void        Path_append_file              (Path *self, const char *appendix);
PUBLIC const char *Path_get                      (Path *self);
PUBLIC const char *Path_get_dir                  (const Path *self);
PUBLIC const char *Path_get_file                 (const Path *self);

PUBLIC bool        Path_file_is_readable         (Path *self);

#endif
