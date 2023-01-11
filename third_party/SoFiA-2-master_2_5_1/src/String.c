// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (String.c) - Source Finding Application                  //
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

/// @file   String.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Container class for handling strings.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "String.h"



/// @brief Container class for handling strings.
///
/// This class provides a convenient structure for storing and handling
/// strings. These strings are dynamic, meaning that they adjust
/// their memory allocation automatically and dynamically as needed,
/// removing the restrictions imposed by the static nature of native
/// C strings.

CLASS String
{
	size_t size;   ///< String size **without** the terminating null character.
	char *string;  ///< C string containing the string value, including '`\0`'.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new String object. Note
/// that the destructor will have to be called explicitly once the
/// object is no longer required to release its memory again.
///
/// @param string  C string to be assigned to the new String object.
///                Use `""` to create an empty string.
///
/// @return Pointer to newly created String object.

PUBLIC String *String_new(const char *string)
{
	// Sanity check
	if(string == NULL) return NULL;
	
	// Allocate memory for new String object
	String *self = (String *)memory(MALLOC, 1, sizeof(String));
	
	// Initialise properties
	self->size = 0;
	self->string = NULL;
	
	// Call method for setting string value
	String_set(self, string);
	
	// Return new String object
	return self;
}



/// @brief Copy constructor
///
/// Copy constructor. Will create a copy of the specified String
/// object. Note that this copy constructor has only been included
/// for semantic reasons, as its sole job is to call the standard
/// constructor.
///
/// @param string  String object to be copied.
///
/// @return Pointer to newly created String object.

PUBLIC String *String_copy(const String *string)
{
	return String_new(string == NULL ? NULL : string->string);
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release all
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void String_delete(String *self)
{
	if(self != NULL) free(self->string);
	free(self);
	return;
}



/// @brief Return string size
///
/// Public method for returning the size of the string excluding
/// the terminating null character.
///
/// @param self  Object self-reference.
///
/// @return Size of the string without the terminating null character.

PUBLIC size_t String_size(const String *self)
{
	return self == NULL ? 0 : self->size;
}



/// @brief Return pointer to C string
///
/// Public method for returning a pointer to the C string stored in
/// `self`. If `self` is `NULL`, then a `NULL` pointer will instead
/// be returned.
///
/// @param self  Object self-reference.
///
/// @return Pointer to C string of `self`, or `NULL` if `self` is `NULL`.

PUBLIC const char *String_get(const String *self)
{
	return self == NULL ? NULL : self->string;
}



/// @brief Extract character at specified index
///
/// Public method for returning the character at the specified index.
/// The process will be terminated if the index is found to be out
/// of range.
///
/// @param self   Object self-reference.
/// @param index  Index of the character to be extracted.
///
/// @return Character at the specified index.

PUBLIC char String_at(const String *self, const size_t index)
{
	// Sanity checks
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "String index out of range.");
	
	return self->string[index];
}



/// @brief Compare two strings
///
/// Public method for checking if two strings are equal, i.e. if 
/// they contain the same character sequence. The method will return
/// `true` if the strings are equal and `false` otherwise.
///
/// @param self    Object self-reference.
/// @param string  C string to compare with.
///
/// @return `true` if strings are equal, `false` otherwise.

PUBLIC bool String_compare(const String *self, const char *string)
{
	// Sanity checks
	check_null(self);
	check_null(string);
	
	return strcmp(self->string, string) == 0;
}



/// @brief Set string
///
/// Public method for setting a string to the specified value. The
/// string will be cleared if the specified value is an empty
/// string or a `NULL` pointer.
///
/// @param self    Object self-reference.
/// @param string  C string value to be set.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_set(String *self, const char *string)
{
	// Sanity checks
	check_null(self);
	
	// Empty string?
	if(string == NULL || *string == '\0')
	{
		String_clear(self);
	}
	else
	{
		self->size = strlen(string);
		self->string = (char *)memory_realloc(self->string, self->size + 1, sizeof(char));
		strcpy(self->string, string);
	}
	
	return self;
}



/// @brief Set character at specified index
///
/// Public method for setting the character at the specified index
/// position to the value of `c`. The method will terminate if the
/// index is out of range.
///
/// @param self   Object self-reference.
/// @param index  Index of the character to be set.
/// @param c      Character.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_set_char(String *self, const size_t index, const char c)
{
	// Sanity checks
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "String index out of range.");
	self->string[index] = c;
	return self;
}



/// @brief Set string from integer
///
/// Public method for setting a string to contain a textual
/// representation of the specified integer value. The size of
/// the resulting string is restricted to at most 32 characters
/// which should be sufficient to hold all possible integer values.
/// A warning would be issued if the buffer size was insufficient
/// and the string be truncated to 32 characters (including the
/// terminating null character) in this case.
///
/// @param self    Object self-reference.
/// @param format  Format specifier similar to `printf()`.
/// @param value   Integer value to set the string to.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_set_int(String *self, const char *format, const long int value)
{
	// Sanity checks
	check_null(self);
	
	// Create buffer
	char *buffer = (char *)memory(MALLOC, 32, sizeof(char));
	
	// Write number into buffer
	const int flag = snprintf(buffer, 32, format != NULL ? format : "%ld", value);
	if(flag >= 32) warning("Buffer overflow in int-to-string conversion.");
	if(flag <   0) warning("Encoding error in int-to-string conversion.");
	
	// Copy buffer into string
	String_set(self, buffer);
	
	// Clean up
	free(buffer);
	
	return self;
}



/// @brief Set string from or until delimiting character
///
/// Public method for setting a String object to the specified
/// string from or until the specified delimiting character. The
/// delimiting character can either be the first from the start or
/// the last before the end (argument `first`), and either the
/// sub-string until the delimiter or from the delimiter onward can
/// be copied (argument `until`). The delimiting character itself
/// will be excluded in all cases. This method can be used to copy
/// part of a string as defined by a specific delimiting character.
///
/// @param  self       Object self-reference.
/// @param  string     C string value to set the String object to.
/// @param  delimiter  Delimiting character from or up to which to
///                    copy the string value.
/// @param  first      If `true`, use the first occurrence of the
///                    delimiter (first from start), otherwise use
///                    the last occurrence (last before end).
/// @param  until      If `true`, copy the string until the delimiting
///                    character, otherwise copy the string from the
///                    delimiting character onward. The delimiter
///                    itself will be excluded in both cases.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_set_delim(String *self, const char *string, const char delimiter, const bool first, const bool until)
{
	// Sanity checks
	check_null(self);
	check_null(string);
	
	// Some settings
	const size_t size_old = strlen(string);
	size_t size_new = 0;
	
	// Find delimiter
	const char *pos_delim = first ? strchr(string, delimiter) : strrchr(string, delimiter);
	
	// Not found?
	if(pos_delim == NULL) return String_set(self, string);
	
	// Determine size of new sub-string
	if(until) size_new = pos_delim - string;
	else size_new = string + size_old - pos_delim - 1;
	
	// If delimiter at edge, clear string
	if(size_new == 0) return String_clear(self);
	
	// Otherwise copy sub-string
	char *tmp = (char *)memory(MALLOC, size_new + 1, sizeof(char));
	if(until) strncpy(tmp, string, size_new);
	else strncpy(tmp, pos_delim + 1, size_new);
	*(tmp + size_new) = '\0';
	
	String_set(self, tmp);
	free(tmp);
	
	return self;
}



/// @brief Append string
///
/// Public method for appending a String object with the specified
/// C string value. The size of the String object will be automatically
/// adjusted.
///
/// @param self    Object self-reference.
/// @param string  C string value to be appended.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_append(String *self, const char *string)
{
	// Sanity checks
	check_null(self);
	if(string == NULL || *string == '\0') return self;
	
	const size_t size = strlen(string);
	
	self->string = (char *)memory_realloc(self->string, self->size + size + 1, sizeof(char));
	strcpy(self->string + self->size, string);
	self->size += size;
	
	return self;
}



/// @brief Append integer value to string
///
/// Public method for appending the textual representation of the
/// specified integer value to a string. The size of the resulting
/// string is restricted to at most 32 characters which should be
/// sufficient to hold all possible integer values. A warning would
/// be issued if the buffer size was insufficient and the string be
/// truncated to 32 characters (including the terminating null
/// character) in this case.
///
/// @param self    Object self-reference.
/// @param format  Format specifier similar to `printf()`.
/// @param value   Integer value the textual representation of
///                which will be appended to the string.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_append_int(String *self, const char *format, const long int value)
{
	// Sanity checks
	check_null(self);
	
	// Create buffer
	char *buffer = (char *)memory(MALLOC, 32, sizeof(char));
	
	// Write number into buffer
	const int flag = snprintf(buffer, 32, format != NULL ? format : "%ld", value);
	if(flag >= 32) warning("Buffer overflow in int-to-string conversion.");
	if(flag <   0) warning("Encoding error in int-to-string conversion.");
	
	// Append buffer to string
	String_append(self, buffer);
	
	// Clean up
	free(buffer);
	
	return self;
}



/// @brief Append floating-point value to string
///
/// Public method for appending the textual representation of the
/// specified floating-point value to a string. The format can be
/// specified as in the `printf()` function. If no format is given
/// (i.e. format is a `NULL` pointer) then `%.5f` will be used by
/// default. The resulting string is restricted to 32 characters
/// (including the terminating null character) and will be truncated
/// with a warning message if exceeded.
///
/// @param self    Object self-reference.
/// @param format  Format specifier as used by `printf()`. If set
///                to `NULL`, the default format will be applied.
/// @param value   Floating-point value the textual representation
///                of which will be appended to the string.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_append_flt(String *self, const char *format, const double value)
{
	// Sanity checks
	check_null(self);
	
	// Create buffer
	char *buffer = (char *)memory(MALLOC, 32, sizeof(char));
	
	// Write number into buffer
	const int flag = snprintf(buffer, 32, format != NULL ? format : "%.5f", value);
	if(flag >= 32) warning("Buffer overflow in float-to-string conversion.");
	if(flag <   0) warning("Encoding error in float-to-string conversion.");
	
	// Append buffer to string
	String_append(self, buffer);
	
	// Clean up
	free(buffer);
	
	return self;
}



/// @brief Prepend string
///
/// Public method for prepending a string with the specified string
/// value. The size of the string will be automatically adjusted.
///
/// @param self    Object self-reference.
/// @param string  C string to be prepended.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_prepend(String *self, const char *string)
{
	// Sanity checks
	check_null(self);
	if(string == NULL || *string == '\0') return self;
	
	const size_t size = strlen(string);
	
	self->string = (char *)memory_realloc(self->string, self->size + size + 1, sizeof(char));
	memmove(self->string + size, self->string, self->size + 1);
	memcpy(self->string, string, size);
	self->size += size;
	
	return self;
}



/// @brief Clear string
///
/// Public method for clearing the specified String object. This
/// will set the string to an empty C string ("`\0`") and the
/// size of the string to zero.
///
/// @param self  Object self-reference.
///
/// @return Pointer to `self` after assignment.

PUBLIC String *String_clear(String *self)
{
	// Sanity checks
	check_null(self);
	
	self->size = 0;
	self->string = (char *)memory_realloc(self->string, 1, sizeof(char));
	*self->string = '\0';
	
	return self;
}



/// @brief Trim string
///
/// Public method for trimming a string by removing any contiguous
/// sequence of whitespace characters from the beginning and end of
/// the string. Whitespace will be anything considered as white-
/// space by the standard library function `isspace()`.
///
/// @param self  Object self-reference.
///
/// @return Pointer to `self` after trimming.

PUBLIC String *String_trim(String *self)
{
	// Sanity checks
	if(self == NULL || self->size == 0) return self;
	
	// Find first non-whitespace character
	char *start = self->string;
	while(isspace((unsigned char)*start)) ++start;
	
	// All space?
	if(*start == '\0')
	{
		String_clear(self);
		return self;
	}
	
	// Find last non-whitespace character
	char *end = self->string + self->size - 1;
	while(end > start && isspace((unsigned char)*end)) --end;
	
	// Shift sub-string to beginning of string buffer
	self->size = end - start + 1;
	memmove(self->string, start, self->size);
	*(self->string + self->size) = '\0';
	
	// Adjust memory allocation
	self->string = (char *)memory_realloc(self->string, self->size + 1, sizeof(char));
	
	return self;
}



/// @brief Convert string to lower case
///
/// Public method for converting a string to lower case. This makes
/// use of the `tolower()` function from the C standard library. A
/// pointer to `self` after conversion will be returned to allow
/// chaining of methods.
///
/// @param self  Object self-reference.
///
/// @return Pointer to `self` after conversion.

PUBLIC String *String_to_lower(String *self)
{
	// Sanity checks
	if(self == NULL || self->size == 0) return self;
	
	// Loop over string to convert to lower-case
	for(char *ptr = self->string; *ptr; ++ptr) *ptr = tolower(*ptr);
	
	return self;
}



/// @brief Convert string to upper case
///
/// Public method for converting a string to upper case. This makes
/// use of the `toupper()` function from the C standard library. A
/// pointer to `self` after conversion will be returned to allow
/// chaining of methods.
///
/// @param self  Object self-reference.
///
/// @return Pointer to `self` after conversion.

PUBLIC String *String_to_upper(String *self)
{
	// Sanity checks
	if(self == NULL || self->size == 0) return self;
	
	// Loop over string to convert to upper-case
	for(char *ptr = self->string; *ptr; ++ptr) *ptr = toupper(*ptr);
	
	return self;
}
