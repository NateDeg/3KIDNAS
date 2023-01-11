// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (common.h) - Source Finding Application                  //
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

/// @file   common.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  This file defines functionality commonly required by all classes (header).


#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#define SOFIA_VERSION "2.5.1"              ///< SoFiA version number.
#define SOFIA_VERSION_FULL "SoFiA 2.5.1"   ///< Full SoFiA version string.
#define SOFIA_CREATION_DATE "24-Jun-2022"  ///< Date of current SoFiA version.

#ifndef M_PI
#define M_PI 3.141592653589793  ///< Archimedes' constant (pi).
#endif

#ifndef MAD_TO_STD
#define MAD_TO_STD 1.482602218505602  ///< Conversion factor between MAD and standard deviation of normal distribution.
#endif
// NOTE: Calculated as 1.0 / scipy.stats.norm.ppf(3.0 / 4.0)

#ifndef INV_SQRT_TWO_PI
#define INV_SQRT_TWO_PI 0.3989422804014327  ///< Value of 1 / sqrt(2 * pi).
#endif

#ifndef IS_NAN
#define IS_NAN(x) ((x) != (x))  ///< Returns `true` if `x` is Not a Number (NaN).
#endif
#ifndef IS_NOT_NAN
#define IS_NOT_NAN(x) ((x) == (x))  ///< Returns `false` if `x` is Not a Number (NaN).
#endif
#ifndef FILTER_NAN
#define FILTER_NAN(x) ((x) == (x) ? (x) : 0)  ///< Returns zero if `x` is Not a Number (NaN); otherwise returns `x`.
#endif

#ifndef IS_ODD
#define IS_ODD(x) ((x) & 1)  ///< Returns `true` if `x` is an odd number.
#endif
#ifndef IS_EVEN
#define IS_EVEN(x) (!((x) & 1))  ///< Returns `true` if `x` is an even number.
#endif

enum {MALLOC, CALLOC};  // Define memory allocation modes

#define NOISE_SAMPLE_SIZE 999983  ///< Define maximum RMS measurement sample size.
// NOTE: This is chosen to be a prime number to reduce the risk of
//       obtaining a stride that is a multiple of the x-axis size.

#define KILOBYTE       1024  ///< Size of a kilobyte (in bytes).
#define MEGABYTE    1048576  ///< Size of a megabyte (in bytes).
#define GIGABYTE 1073741824  ///< Size of a gigabyte (in bytes).

// Define object-oriented terminology
#define CLASS struct    ///< Define `CLASS` as alias of `struct` (in analogy to C++ '`class`').
#define PUBLIC extern   ///< Define `PUBLIC` as alias of `extern` (in analogy to C++ '`public`').
#define PRIVATE static  ///< Define `PRIVATE` as alias of `static` (in analogy to C++ '`private`').

// Define error codes
#define ERR_SUCCESS      0  ///< Return code on success.
#define ERR_FAILURE      1  ///< Return code on failure; indicating that an unspecified error occurred.
#define ERR_NULL_PTR     2  ///< Return code on failure; indicating an attempt to dereference `NULL` pointer.
#define ERR_MEM_ALLOC    3  ///< Return code on failure; indicating a memory allocation error.
#define ERR_INDEX_RANGE  4  ///< Return code on failure; indicating that an array index is out of range.
#define ERR_FILE_ACCESS  5  ///< Return code on failure; indicating a file access error.
#define ERR_INT_OVERFLOW 6  ///< Return code on failure; indicating an integer overflow error.
#define ERR_USER_INPUT   7  ///< Return code on failure; indicating invalid user input.
#define ERR_NO_SRC_FOUND 8  ///< Return code on failure; indicating that no sources were found by SoFiA.

// Generic compile time check; should result in a compiler error if
// condition is false due to attempt to create array of negative size.
// NOTE: This does not actually create a physical array, but merely
//       defines a new type.
#define COMPILE_TIME_CHECK(condition, message) typedef char message[(condition) ? 1 : -1]  ///< Check `condition` at compile time and terminate compilation if `false`.

// Check condition and exit if not met
void ensure(const bool condition, const int errorCode, const char *format, ...);
void check_null(const void *ptr);

// Print info and warning messages
void message(const char *format, ...);
void message_verb(const bool verbosity, const char *format, ...);
void status(const char *format, ...);
void warning(const char *format, ...);
void warning_verb(const bool verbosity, const char *format, ...);

// Display progress bar and time stamp
void progress_bar(const char *text, const size_t progress, const size_t maximum);
void timestamp(const time_t start, const clock_t start_clock);

// Memory allocation
void *memory(const int mode, const size_t n_blocks, const size_t block_size);
void *memory_realloc(void *ptr, const size_t n_blocks, const size_t block_size);

// String functions
char *trim_string(char *str);
void int_to_str(char *str, const size_t size, const long int value);

// Swap two values
void swap(double *val1, double *val2);

// Plotting aids
double auto_tick(const double range, const size_t n);
void write_eps_header(FILE *fp, const char *title, const char *creator, const char *bbox);
void write_eps_footer(FILE *fp);

// Byte-order functions
bool is_little_endian(void);
void swap_byte_order(char *word, const size_t size);

#endif
