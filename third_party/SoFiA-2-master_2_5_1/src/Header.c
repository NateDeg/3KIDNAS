// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Header.c) - Source Finding Application                  //
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

/// @file   Header.c
/// @author Tobias Westmeier
/// @date   26/11/2021
/// @brief  Class for storing and managing FITS file headers.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Header.h"


/// @brief Class for storing and managing FITS file headers.
///
/// The purpose of this class is to provide a structure for storing
/// header information from a FITS data file and implement ways of
/// reading and manipulating individual header entries.

CLASS Header
{
	char   *header;     ///< Pointer to `char` array containing the header block.
	size_t  size;       ///< Size of header block in bytes.
	bool    verbosity;  ///< Verbosity level (0 or 1).
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new Header object occupied
/// with the information provided through the `header` argument and
/// return a pointer to the newly created object. Note that the
/// destructor will need to be called explicitly once the object is
/// no longer required to release any memory allocated during the
/// lifetime of the object.
///
/// @param header     Character array containing the header
///                   information; must not be null-terminated.
/// @param size       Size of the header; must be a multiple of
///                   the native FITS header block size.
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created Header object.

PUBLIC Header *Header_new(const char *header, const size_t size, const bool verbosity)
{
	// Sanity checks
	check_null(header);
	ensure(size, ERR_USER_INPUT, "Received empty header array.");
	
	// Allocate memory for new Header object
	Header *self = (Header *)memory(MALLOC, 1, sizeof(Header));
	
	// Make copy of header data
	char *header_copy = (char *)memory(MALLOC, size, sizeof(char));
	memcpy(header_copy, header, size);
	
	// Initialise properties
	self->size      = size;
	self->header    = header_copy;
	self->verbosity = verbosity;
	
	return self;
}



/// @brief Copy constructor
///
/// Copy constructor. Will create a new Header object that is a
/// physical copy of the object pointed to by `source`. A pointer
/// to the newly created object will be returned. Note that the
/// destructor will need to be called explicitly once the object
/// is no longer required to release any memory allocated to the
/// object.
///
/// @param source  Header object to be copied.
///
/// @return Pointer to newly created Header object.

PUBLIC Header *Header_copy(const Header *source)
{
	check_null(source);
	return Header_new(source->header, source->size, source->verbosity);
}



/// @brief Variant of standard constructor
///
/// Alternative standard constructor. Will create a new Header object
/// that contains an empty header, i.e. a single FITS header block
/// containing just the `END` keyword. A pointer to the newly
/// created object will be returned. Note that the destructor will
/// need to be called explicitly once the object is no longer
/// required to release any memory allocated to the object.
///
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created Header object.

PUBLIC Header *Header_blank(const bool verbosity)
{
	// Allocate memory for new Header object
	Header *self = (Header *)memory(MALLOC, 1, sizeof(Header));
	
	// Create basic header
	self->verbosity = verbosity;
	self->size      = FITS_HEADER_BLOCK_SIZE;
	self->header    = (char *)memory(MALLOC, self->size, sizeof(char));
	
	// Fill entire header with spaces
	memset(self->header, ' ', self->size);
	
	// Insert END keyword
	memcpy(self->header, "END", 3);
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Header_delete(Header *self)
{
	if(self != NULL) free(self->header);
	free(self);
	return;
}



/// @brief Return pointer to header data
///
/// Public method for returning a pointer to the header
/// data of the Header object.
///
/// @param self  Object self-reference.
///
/// @return Pointer to header data.

PUBLIC const char *Header_get(const Header *self)
{
	check_null(self);
	return self->header;
}



/// @brief Return size of header
///
/// Public method for returning the size of the Header object, i.e.
/// the number of characters occupied by the header data.
///
/// @param self  Object self-reference.
///
/// @return Size of header.

PUBLIC size_t Header_get_size(const Header *self)
{
	check_null(self);
	return self->size;
}



/// @brief Retrieve header element as raw string buffer
///
/// Private method for retrieving the specified header element as a
/// raw string buffer. The resulting value will be written to the
/// `char` array pointed to by `buffer`, which needs to be large enough
/// to hold the maximum permissible FITS header value size. If the
/// header keyword is not found, the buffer will remain unchanged
/// and a value of 1 will be returned by the function.
///
/// @param self    Object self-reference.
/// @param key     Name of the header element to be retrieved.
/// @param buffer  Pointer to `char` buffer for holding result.
///
/// @return Returns 0 on success or 1 if the header keyword was not found.

PRIVATE int Header_get_raw(const Header *self, const char *key, char *buffer)
{
	// Sanity checks (done here, as the respective public methods all call this function)
	check_null(self);
	check_null(self->header);
	check_null(buffer);
	check_null(key);
	
	const char *ptr = self->header;
	
	while(ptr < self->header + self->size)
	{
		if(strncmp(ptr, key, strlen(key)) == 0)
		{
			memcpy(buffer, ptr + FITS_HEADER_KEY_SIZE, FITS_HEADER_VALUE_SIZE);
			buffer[FITS_HEADER_VALUE_SIZE] = '\0';
			return 0;
		}
		
		ptr += FITS_HEADER_LINE_SIZE;
	}
	
	warning_verb(self->verbosity, "Header keyword \'%s\' not found.", key);
	return 1;
}



/// @brief Retrieve header element as `int`
///
/// Public method for retrieving the specified header element as an
/// Integer value. The method will call Header_get_raw(); see there
/// for more information.
///
/// @param self  Object self-reference.
/// @param key   Name of the header element to be retrieved.
///
/// @return Returns the requested header value as an Integer value.
///         If the header keyword was not found, then 0 will be
///         returned.

PUBLIC long int Header_get_int(const Header *self, const char *key)
{
	char buffer[FITS_HEADER_VALUE_SIZE + 1] = ""; // Note that "" initialises entire array with 0
	const int flag = Header_get_raw(self, key, buffer);
	return flag ? 0 : strtol(buffer, NULL, 10);
}



/// @brief Retrieve header element as `double`
///
/// Public method for retrieving the specified header element as
/// a double-precision floating-point value. The method will call
/// Header_get_raw(); see there for more information.
///
/// @param self  Object self-reference.
/// @param key   Name of the header element to be retrieved.
///
/// @return Returns the requested header value as a double-precision
///         floating-point value. If the header keyword was not found,
///         then `NaN` will be returned.

PUBLIC double Header_get_flt(const Header *self, const char *key)
{
	char buffer[FITS_HEADER_VALUE_SIZE + 1] = "";
	const int flag = Header_get_raw(self, key, buffer);
	return flag ? NAN : strtod(buffer, NULL);
}



/// @brief Retrieve header element as `bool`
///
/// Public method for retrieving the specified header element as a
/// Boolean value. The method will call Header_get_raw(); see there
/// for more information.
///
/// @param self  Object self-reference.
/// @param key   Name of the header element to be retrieved.
///
/// @return Returns the requested header value as a Boolean value.
///         If the header keyword was not found, then `false` will
///         be returned.

PUBLIC bool Header_get_bool(const Header *self, const char *key)
{
	char buffer[FITS_HEADER_VALUE_SIZE + 1] = "";
	const int flag = Header_get_raw(self, key, buffer);
	
	if(!flag)
	{
		const char *ptr = buffer;
		while(ptr <= buffer + FITS_HEADER_VALUE_SIZE)
		{
			if(*ptr == ' ') ++ptr;
			else return *ptr == 'T';
		}
	}
	
	return false;
}



/// @brief Retrieve header element as string
///
/// Public method for retrieving the specified header element as a
/// string. The string value will be written to the `char` array
/// pointed to by `value`, which will need to be large enough to hold
/// the maximum permissible FITS header value size. This function
/// will call Header_get_raw(); see there for more information.
/// If the header keyword is not found, `value` will be set to an
/// empty string and a value of 1 will be returned.
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be retrieved.
/// @param value  Pointer to `char` buffer for holding result.
///
/// @return Returns 0 on success and 1 if the header keyword was not found.

PUBLIC int Header_get_str(const Header *self, const char *key, char *value)
{
	// WARNING: This function will fail if there are quotation marks inside a comment!
	char buffer[FITS_HEADER_VALUE_SIZE + 1] =  "";
	if(Header_get_raw(self, key, buffer))
	{
		value[0] = '\0';
		return 1;
	}
	
	const char *left = strchr(buffer, '\'');
	ensure(left != NULL, ERR_USER_INPUT, "FITS header entry is not a string.");
	
	const char *right = strchr(left + 1, '\'');
	while(right != NULL && *(right + 1) == '\'') right = strchr(right + 2, '\'');
	ensure(right != NULL, ERR_USER_INPUT, "Unbalanced quotation marks in FITS header entry.");
	
	memcpy(value, left + 1, right - left - 1);
	value[right - left - 1] = '\0';
	return 0;
}



/// @brief Retrieve header element as String object
///
/// Public method for retrieving the specified header element as a
/// String object. This function will call Header_get_raw(); see there
/// for more information. If the header keyword is not found, the
/// returned String object will contain an empty string.
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be retrieved.
///
/// @return String object containing the requested header element.

PUBLIC String *Header_get_string(const Header *self, const char *key)
{
	// WARNING: This function will fail if there are quotation marks inside a comment!
	char buffer[FITS_HEADER_VALUE_SIZE + 1] =  "";
	
	if(Header_get_raw(self, key, buffer))
	{
		// Keyword doesn't exist --> return empty string
		buffer[0] = '\0';
	}
	else
	{
		// Keyword does exist --> remove quotation marks
		const char *left = strchr(buffer, '\'');
		ensure(left != NULL, ERR_USER_INPUT, "FITS header entry is not a string.");
		
		const char *right = strchr(left + 1, '\'');
		while(right != NULL && *(right + 1) == '\'') right = strchr(right + 2, '\'');
		ensure(right != NULL, ERR_USER_INPUT, "Unbalanced quotation marks in FITS header entry.");
		
		memmove(buffer, left + 1, right - left - 1);
		buffer[right - left - 1] = '\0';
	}
	
	// Create and return String object
	return String_new(buffer);
}



/// @brief Write raw C string to header
///
/// Private method for writing a raw string buffer into the header.
/// The buffer needs to of size `FITS_HEADER_VALUE_SIZE` and padded
/// with spaces (ASCII 32) at the end. If the specified keyword
/// already exists, its first occurrence will be overwritten with the
/// new buffer. If the keyword does not exist, a new entry will be
/// inserted at the end of the header just before the `END` keyword.
/// If necessary, the header size will be automatically adjusted to
/// be able to accommodate the new entry.
///
/// @param self    Object self-reference.
/// @param key     Name of the header element to be written.
/// @param buffer  Character buffer to be written to header.
///
/// @return Returns 0 if a new header entry was created and 1 if an
///         existing header entry was overwritten.
///
/// @note `COMMENT` and `HISTORY` items will always be appended at
///       the end, and existing ones will never be overwritten.

PRIVATE int Header_set_raw(Header *self, const char *key, const char *buffer)
{
	// Sanity checks
	check_null(self);
	check_null(self->header);
	check_null(key);
	check_null(buffer);
	ensure(strlen(key) > 0 && strlen(key) <= FITS_HEADER_KEYWORD_SIZE, ERR_USER_INPUT, "Illegal length of header keyword.");
	
	const bool is_not_comment = (strcmp(key, "COMMENT") != 0);
	const bool is_not_history = (strcmp(key, "HISTORY") != 0);
	
	char *ptr = self->header;
	size_t line = 0;
	
	if(is_not_comment && is_not_history)
	{
		// Check if keyword already exists
		line = Header_check(self, key);
		
		// Overwrite header entry if already present
		if(line > 0)
		{
			memcpy(ptr + (line - 1) * FITS_HEADER_LINE_SIZE + FITS_HEADER_KEY_SIZE, buffer, FITS_HEADER_VALUE_SIZE);
			return 0;
		}
	}
	
	// Create a new entry
	warning_verb(self->verbosity, "Header keyword \'%s\' not found. Creating new entry.", key);
	
	// Check current length
	line = Header_check(self, "END");
	ensure(line > 0, ERR_USER_INPUT, "No END keyword found in header of Header object.");
	
	// Expand header if necessary
	if(line % FITS_HEADER_LINES == 0)
	{
		warning_verb(self->verbosity, "Expanding header to fit new entry.");
		self->size += FITS_HEADER_BLOCK_SIZE;
		self->header = (char *)memory_realloc(self->header, self->size, sizeof(char));
		memset(self->header + self->size - FITS_HEADER_BLOCK_SIZE, ' ', FITS_HEADER_BLOCK_SIZE); // fill with space
	}
	
	ptr = self->header;
	
	// Add new header keyword at end
	memcpy(ptr + (line - 1) * FITS_HEADER_LINE_SIZE, key, strlen(key)); // key
	if(is_not_comment && is_not_history) memcpy(ptr + (line - 1) * FITS_HEADER_LINE_SIZE + FITS_HEADER_KEYWORD_SIZE, "=", 1); // =
	memcpy(ptr + (line - 1) * FITS_HEADER_LINE_SIZE + FITS_HEADER_KEY_SIZE, buffer, FITS_HEADER_VALUE_SIZE); // value
	memcpy(ptr + line * FITS_HEADER_LINE_SIZE, "END", 3); // new end
	
	return 1;
}



/// @brief Write `int` value to header
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be written.
/// @param value  Value to be written to header.
///
/// @return Returns 0 if a new header entry was created and 1 if an
///         existing header entry was overwritten.
///
/// Public method for writing an `int` value into the header. The
/// method will call Header_set_raw(); see there for more information.

PUBLIC int Header_set_int(Header *self, const char *key, const long int value)
{
	char buffer[FITS_HEADER_VALUE_SIZE];
	memset(buffer, ' ', FITS_HEADER_VALUE_SIZE);
	int size = snprintf(buffer, FITS_HEADER_FIXED_WIDTH + 1, "%20ld", value);
	ensure(size > 0 && size <= FITS_HEADER_FIXED_WIDTH, ERR_FAILURE, "Creation of new header entry failed for unknown reasons.");
	buffer[size] = ' '; // get rid of NUL character
	
	return Header_set_raw(self, key, buffer);
}



/// @brief Write `double` value to header
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be written.
/// @param value  Value to be written to header.
///
/// @return Returns 0 if a new header entry was created and 1 if an
///         existing header entry was overwritten.
///
/// Public method for writing a `double` value into the header. The
/// method will call Header_set_raw(); see there for more information.

PUBLIC int Header_set_flt(Header *self, const char *key, const double value)
{
	char buffer[FITS_HEADER_VALUE_SIZE];
	memset(buffer, ' ', FITS_HEADER_VALUE_SIZE);
	int size = snprintf(buffer, FITS_HEADER_FIXED_WIDTH + 1, "%20.11E", value);
	ensure(size > 0 && size <= FITS_HEADER_FIXED_WIDTH, ERR_FAILURE, "Creation of new header entry failed for unknown reasons.");
	buffer[size] = ' '; // get rid of NUL character
	
	return Header_set_raw(self, key, buffer);
}



/// @brief Write `bool` value to header
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be written.
/// @param value  Value to be written to header.
///
/// @return Returns 0 if a new header entry was created and 1 if an
///         existing header entry was overwritten.
///
/// Public method for writing a `bool` value into the header. The
/// method will call Header_set_raw(); see there for more information.

PUBLIC int Header_set_bool(Header *self, const char *key, const bool value)
{
	char buffer[FITS_HEADER_VALUE_SIZE];
	memset(buffer, ' ', FITS_HEADER_VALUE_SIZE);
	buffer[FITS_HEADER_FIXED_WIDTH - 1] = value ? 'T' : 'F';
	
	return Header_set_raw(self, key, buffer);
}



/// @brief Write string value to header
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to be written.
/// @param value  Value to be written to header.
///
/// @return Returns 0 if a new header entry was created and 1 if an
///         existing header entry was overwritten.
///
/// Public method for writing a string value into the header. The
/// method will call Header_set_raw(); see there for more information.

PUBLIC int Header_set_str(Header *self, const char *key, const char *value)
{
	const size_t size = strlen(value);
	ensure(size <= FITS_HEADER_VALUE_SIZE - 2, ERR_USER_INPUT, "String too long for FITS header line.");
	char buffer[FITS_HEADER_VALUE_SIZE];
	memset(buffer, ' ', FITS_HEADER_VALUE_SIZE);
	memcpy(buffer, "\'", 1); // opening quotation mark
	memcpy(buffer + 1, value, size); // string
	memcpy(buffer + 1 + size, "\'", 1); // closing quotation mark
	
	return Header_set_raw(self, key, buffer);
}



/// @brief Write `COMMENT` or `HISTORY` item to header
///
/// Public method for writing a `COMMENT` or `HISTORY` item into the
/// header. If the `value` string exceeds the maximum permissible
/// length, then the comment will automatically be broken into multiple
/// lines. The number of lines written will be returned. If `history`
/// is set to `true`, then a `HISTORY` keyword will be created.
/// Otherwise, a standard `COMMENT` keyword will be written.
///
/// @param self     Object self-reference.
/// @param value    Comment or history string to be written.
/// @param history  If `true`, then a `HISTORY` item will be created,
///                 otherwise a standard `COMMENT` will be written.
///
/// @return Returns the number of lines written.

PUBLIC size_t Header_comment(Header *self, const char *value, const bool history)
{
	size_t size = strlen(value);
	const size_t lines = (size == 0) ? 1 : (size - 1) / FITS_HEADER_VALUE_SIZE + 1;
	char buffer[FITS_HEADER_VALUE_SIZE];
	
	for(size_t i = 0; i < lines; ++i)
	{
		memset(buffer, ' ', FITS_HEADER_VALUE_SIZE);
		if(size > FITS_HEADER_VALUE_SIZE)
		{
			memcpy(buffer, value + i * FITS_HEADER_VALUE_SIZE, FITS_HEADER_VALUE_SIZE);
			size -= FITS_HEADER_VALUE_SIZE;
		}
		else if(size > 0) memcpy(buffer, value + i * FITS_HEADER_VALUE_SIZE, size);
		Header_set_raw(self, history ? "HISTORY" : "COMMENT", buffer);
	}
	
	return lines;
}



/// @brief Check for existence of header keyword
///
/// Searches for the first occurrence of the specified header keyword
/// and returns the corresponding line number. If the header keyword
/// is not found, the function will return 0.
///
/// @param self  Object self-reference.
/// @param key   Name of the header element to be checked for.
///
/// @return Line number of the first occurrence of the specified key
///         in the header. If the key is not found, 0 is returned.

PUBLIC size_t Header_check(const Header *self, const char *key)
{
	// Sanity checks
	check_null(self);
	check_null(self->header);
	check_null(key);
	const size_t size = strlen(key);
	ensure(size > 0 && size <= FITS_HEADER_KEYWORD_SIZE, ERR_USER_INPUT, "Illegal FITS header keyword: %s.", key);
	
	char *ptr = self->header;
	size_t line = 1;
	
	while(ptr < self->header + self->size)
	{
		if(strncmp(ptr, key, size) == 0 && (*(ptr + size) == ' ' || *(ptr + size) == '=')) return line;
		ptr += FITS_HEADER_LINE_SIZE;
		++line;
	}
	
	warning_verb(self->verbosity, "Header keyword \'%s\' not found.", key);
	return 0;
}



/// @brief Check if header value equal to specified string
///
/// Public method for comparing the first `n` characters of the value
/// of the header entry specified by `key` to the specified string
/// value. Returns `true` if the two are the same, `false` otherwise.
///
/// @param self   Object self-reference.
/// @param key    Name of the header element to check.
/// @param value  String value to compare against.
/// @param n      Number of characters to compare. Set to 0 to
///               compare all characters.
///
/// @return `true` if the values are equal, `false` otherwise.
///
/// @note This will only work for strings, but not for any other
///       header data type.

PUBLIC bool Header_compare(const Header *self, const char *key, const char *value, const size_t n)
{
	// Sanity checks
	check_null(self);
	check_null(key);
	check_null(value);
	
	char buffer[FITS_HEADER_VALUE_SIZE + 1];
	Header_get_str(self, key, buffer);
	
	if(n) return strncmp(buffer, value, n) == 0;
	return strcmp(buffer, value) == 0;
}



/// @brief Delete header keyword
///
/// Deletes all occurrences of the specified header keyword. Any
/// empty blocks at the end of the new header will be removed.
///
/// @param self  Object self-reference.
/// @param key   Name of the header element to be deleted.
///
/// @return Returns 1 if the header keyword was not found, 0 otherwise.

PUBLIC int Header_remove(Header *self, const char *key)
{
	size_t line = Header_check(self, key);
	if(line == 0) return 1;
	
	// Header keyword found; shift all subsequent lines up by 1
	// and fill last line with spaces. Do this repeatedly until
	// the keyword is no longer found.
	while(line)
	{
		memmove(self->header + (line - 1) * FITS_HEADER_LINE_SIZE, self->header + line * FITS_HEADER_LINE_SIZE, self->size - line * FITS_HEADER_LINE_SIZE);
		memset(self->header + self->size - FITS_HEADER_LINE_SIZE, ' ', FITS_HEADER_LINE_SIZE);
		line = Header_check(self, key);
	}
	
	// Check if the header block can be shortened.
	line = Header_check(self, "END");
	ensure(line, ERR_USER_INPUT, "END keyword missing from FITS header.");
	const size_t last_line = self->size / FITS_HEADER_LINE_SIZE;
	const size_t empty_blocks = (last_line - line) / FITS_HEADER_LINES;
	
	if(empty_blocks)
	{
		warning_verb(self->verbosity, "Reducing size of header to remove empty block(s).");
		self->size -= empty_blocks * FITS_HEADER_BLOCK_SIZE;
		self->header = (char *)memory_realloc(self->header, self->size, sizeof(char));
	}
	
	return 0;
}



/// @brief Copy WCS information from one header to another
///
/// Public method for copying WCS header information from one data
/// cube to another. This method is intended to be used when a
/// blank data cube is created with the same dimensions and region
/// on the sky as an existing cube; the WCS information of the
/// existing cube can then simply be copied across to populate the
/// header of the blank cube with the appropriate axis descriptors
/// and WCS keywords.
///
/// Note that this will also work if the blank cube has a reduced
/// dimensionality compared to the original cube, e.g. when a moment
/// map is to be created from a 3-D data cube. Only the relevant
/// axes will be copied in this case based on the `NAXIS` keyword
/// of the target cube.
///
/// @param source  Data cube from which to copy WCS information.
/// @param target  Data cube to which to copy WCS information.

PUBLIC void Header_copy_wcs(const Header *source, Header *target)
{
	// Sanity checks
	check_null(source);
	check_null(target);
	check_null(source->header);
	check_null(target->header);
	
	char value[FITS_HEADER_VALUE_SIZE + 1];
	const size_t dimensions = Header_get_int(target, "NAXIS");
	ensure(dimensions, ERR_USER_INPUT, "\'NAXIS\' keyword is missing from header.");
	
	// First axis
	if(dimensions >= 1)
	{
		if(Header_check(source, "CTYPE1"))
		{
			Header_get_str(source, "CTYPE1", value);
			Header_set_str(target, "CTYPE1", value);
		}
		if(Header_check(source, "CRVAL1")) Header_set_flt(target, "CRVAL1", Header_get_flt(source, "CRVAL1"));
		if(Header_check(source, "CRPIX1")) Header_set_flt(target, "CRPIX1", Header_get_flt(source, "CRPIX1"));
		if(Header_check(source, "CDELT1")) Header_set_flt(target, "CDELT1", Header_get_flt(source, "CDELT1"));
		if(Header_check(source, "CUNIT1"))
		{
			Header_get_str(source, "CUNIT1", value);
			Header_set_str(target, "CUNIT1", value);
		}
		if(Header_check(source, "CROTA1")) Header_set_flt(target, "CROTA1", Header_get_flt(source, "CROTA1"));
	}
	
	// Second axis
	if(dimensions >= 2)
	{
		if(Header_check(source, "CTYPE2"))
		{
			Header_get_str(source, "CTYPE2", value);
			Header_set_str(target, "CTYPE2", value);
		}
		if(Header_check(source, "CRVAL2")) Header_set_flt(target, "CRVAL2", Header_get_flt(source, "CRVAL2"));
		if(Header_check(source, "CRPIX2")) Header_set_flt(target, "CRPIX2", Header_get_flt(source, "CRPIX2"));
		if(Header_check(source, "CDELT2")) Header_set_flt(target, "CDELT2", Header_get_flt(source, "CDELT2"));
		if(Header_check(source, "CUNIT2"))
		{
			Header_get_str(source, "CUNIT2", value);
			Header_set_str(target, "CUNIT2", value);
		}
		if(Header_check(source, "CROTA2")) Header_set_flt(target, "CROTA2", Header_get_flt(source, "CROTA2"));
	}
	
	// Third axis
	if(dimensions >= 3)
	{
		if(Header_check(source, "CTYPE3"))
		{
			Header_get_str(source, "CTYPE3", value);
			Header_set_str(target, "CTYPE3", value);
		}
		if(Header_check(source, "CRVAL3")) Header_set_flt(target, "CRVAL3", Header_get_flt(source, "CRVAL3"));
		if(Header_check(source, "CRPIX3")) Header_set_flt(target, "CRPIX3", Header_get_flt(source, "CRPIX3"));
		if(Header_check(source, "CDELT3")) Header_set_flt(target, "CDELT3", Header_get_flt(source, "CDELT3"));
		if(Header_check(source, "CUNIT3"))
		{
			Header_get_str(source, "CUNIT3", value);
			Header_set_str(target, "CUNIT3", value);
		}
		if(Header_check(source, "CROTA3")) Header_set_flt(target, "CROTA3", Header_get_flt(source, "CROTA3"));
		if(Header_check(source, "CELLSCAL"))
		{
			// NOTE: CELLSCAL keyword will only be copied for 3-D data!
			Header_get_str(source, "CELLSCAL", value);
			Header_set_str(target, "CELLSCAL", value);
		}
	}
	
	// PCi_j, PVi_j and CDi_j keywords
	if(Header_check(source, "PC1_1") && dimensions >= 2) Header_set_flt(target, "PC1_1", Header_get_flt(source, "PC1_1"));
	if(Header_check(source, "PC2_1") && dimensions >= 2) Header_set_flt(target, "PC2_1", Header_get_flt(source, "PC2_1"));
	if(Header_check(source, "PC3_1") && dimensions >= 3) Header_set_flt(target, "PC3_1", Header_get_flt(source, "PC3_1"));
	if(Header_check(source, "PC1_2") && dimensions >= 2) Header_set_flt(target, "PC1_2", Header_get_flt(source, "PC1_2"));
	if(Header_check(source, "PC2_2") && dimensions >= 2) Header_set_flt(target, "PC2_2", Header_get_flt(source, "PC2_2"));
	if(Header_check(source, "PC3_2") && dimensions >= 3) Header_set_flt(target, "PC3_2", Header_get_flt(source, "PC3_2"));
	if(Header_check(source, "PC1_3") && dimensions >= 3) Header_set_flt(target, "PC1_3", Header_get_flt(source, "PC1_3"));
	if(Header_check(source, "PC2_3") && dimensions >= 3) Header_set_flt(target, "PC2_3", Header_get_flt(source, "PC2_3"));
	if(Header_check(source, "PC3_3") && dimensions >= 3) Header_set_flt(target, "PC3_3", Header_get_flt(source, "PC3_3"));
	if(Header_check(source, "PC01_01") && dimensions >= 2) Header_set_flt(target, "PC01_01", Header_get_flt(source, "PC01_01"));
	if(Header_check(source, "PC02_01") && dimensions >= 2) Header_set_flt(target, "PC02_01", Header_get_flt(source, "PC02_01"));
	if(Header_check(source, "PC03_01") && dimensions >= 3) Header_set_flt(target, "PC03_01", Header_get_flt(source, "PC03_01"));
	if(Header_check(source, "PC01_02") && dimensions >= 2) Header_set_flt(target, "PC01_02", Header_get_flt(source, "PC01_02"));
	if(Header_check(source, "PC02_02") && dimensions >= 2) Header_set_flt(target, "PC02_02", Header_get_flt(source, "PC02_02"));
	if(Header_check(source, "PC03_02") && dimensions >= 3) Header_set_flt(target, "PC03_02", Header_get_flt(source, "PC03_02"));
	if(Header_check(source, "PC01_03") && dimensions >= 3) Header_set_flt(target, "PC01_03", Header_get_flt(source, "PC01_03"));
	if(Header_check(source, "PC02_03") && dimensions >= 3) Header_set_flt(target, "PC02_03", Header_get_flt(source, "PC02_03"));
	if(Header_check(source, "PC03_03") && dimensions >= 3) Header_set_flt(target, "PC03_03", Header_get_flt(source, "PC03_03"));
	
	if(Header_check(source, "PV1_1") && dimensions >= 2) Header_set_flt(target, "PV1_1", Header_get_flt(source, "PV1_1"));
	if(Header_check(source, "PV2_1") && dimensions >= 2) Header_set_flt(target, "PV2_1", Header_get_flt(source, "PV2_1"));
	if(Header_check(source, "PV3_1") && dimensions >= 3) Header_set_flt(target, "PV3_1", Header_get_flt(source, "PV3_1"));
	if(Header_check(source, "PV1_2") && dimensions >= 2) Header_set_flt(target, "PV1_2", Header_get_flt(source, "PV1_2"));
	if(Header_check(source, "PV2_2") && dimensions >= 2) Header_set_flt(target, "PV2_2", Header_get_flt(source, "PV2_2"));
	if(Header_check(source, "PV3_2") && dimensions >= 3) Header_set_flt(target, "PV3_2", Header_get_flt(source, "PV3_2"));
	if(Header_check(source, "PV1_3") && dimensions >= 3) Header_set_flt(target, "PV1_3", Header_get_flt(source, "PV1_3"));
	if(Header_check(source, "PV2_3") && dimensions >= 3) Header_set_flt(target, "PV2_3", Header_get_flt(source, "PV2_3"));
	if(Header_check(source, "PV3_3") && dimensions >= 3) Header_set_flt(target, "PV3_3", Header_get_flt(source, "PV3_3"));
	if(Header_check(source, "PV01_01") && dimensions >= 2) Header_set_flt(target, "PV01_01", Header_get_flt(source, "PV01_01"));
	if(Header_check(source, "PV02_01") && dimensions >= 2) Header_set_flt(target, "PV02_01", Header_get_flt(source, "PV02_01"));
	if(Header_check(source, "PV03_01") && dimensions >= 3) Header_set_flt(target, "PV03_01", Header_get_flt(source, "PV03_01"));
	if(Header_check(source, "PV01_02") && dimensions >= 2) Header_set_flt(target, "PV01_02", Header_get_flt(source, "PV01_02"));
	if(Header_check(source, "PV02_02") && dimensions >= 2) Header_set_flt(target, "PV02_02", Header_get_flt(source, "PV02_02"));
	if(Header_check(source, "PV03_02") && dimensions >= 3) Header_set_flt(target, "PV03_02", Header_get_flt(source, "PV03_02"));
	if(Header_check(source, "PV01_03") && dimensions >= 3) Header_set_flt(target, "PV01_03", Header_get_flt(source, "PV01_03"));
	if(Header_check(source, "PV02_03") && dimensions >= 3) Header_set_flt(target, "PV02_03", Header_get_flt(source, "PV02_03"));
	if(Header_check(source, "PV03_03") && dimensions >= 3) Header_set_flt(target, "PV03_03", Header_get_flt(source, "PV03_03"));
	
	if(Header_check(source, "CD1_1") && dimensions >= 2) Header_set_flt(target, "CD1_1", Header_get_flt(source, "CD1_1"));
	if(Header_check(source, "CD2_1") && dimensions >= 2) Header_set_flt(target, "CD2_1", Header_get_flt(source, "CD2_1"));
	if(Header_check(source, "CD3_1") && dimensions >= 3) Header_set_flt(target, "CD3_1", Header_get_flt(source, "CD3_1"));
	if(Header_check(source, "CD1_2") && dimensions >= 2) Header_set_flt(target, "CD1_2", Header_get_flt(source, "CD1_2"));
	if(Header_check(source, "CD2_2") && dimensions >= 2) Header_set_flt(target, "CD2_2", Header_get_flt(source, "CD2_2"));
	if(Header_check(source, "CD3_2") && dimensions >= 3) Header_set_flt(target, "CD3_2", Header_get_flt(source, "CD3_2"));
	if(Header_check(source, "CD1_3") && dimensions >= 3) Header_set_flt(target, "CD1_3", Header_get_flt(source, "CD1_3"));
	if(Header_check(source, "CD2_3") && dimensions >= 3) Header_set_flt(target, "CD2_3", Header_get_flt(source, "CD2_3"));
	if(Header_check(source, "CD3_3") && dimensions >= 3) Header_set_flt(target, "CD3_3", Header_get_flt(source, "CD3_3"));
	if(Header_check(source, "CD01_01") && dimensions >= 2) Header_set_flt(target, "CD01_01", Header_get_flt(source, "CD01_01"));
	if(Header_check(source, "CD02_01") && dimensions >= 2) Header_set_flt(target, "CD02_01", Header_get_flt(source, "CD02_01"));
	if(Header_check(source, "CD03_01") && dimensions >= 3) Header_set_flt(target, "CD03_01", Header_get_flt(source, "CD03_01"));
	if(Header_check(source, "CD01_02") && dimensions >= 2) Header_set_flt(target, "CD01_02", Header_get_flt(source, "CD01_02"));
	if(Header_check(source, "CD02_02") && dimensions >= 2) Header_set_flt(target, "CD02_02", Header_get_flt(source, "CD02_02"));
	if(Header_check(source, "CD03_02") && dimensions >= 3) Header_set_flt(target, "CD03_02", Header_get_flt(source, "CD03_02"));
	if(Header_check(source, "CD01_03") && dimensions >= 3) Header_set_flt(target, "CD01_03", Header_get_flt(source, "CD01_03"));
	if(Header_check(source, "CD02_03") && dimensions >= 3) Header_set_flt(target, "CD02_03", Header_get_flt(source, "CD02_03"));
	if(Header_check(source, "CD03_03") && dimensions >= 3) Header_set_flt(target, "CD03_03", Header_get_flt(source, "CD03_03"));
	
	// Rest frequency and velocity
	if(Header_check(source, "RESTFREQ")) Header_set_flt(target, "RESTFREQ", Header_get_flt(source, "RESTFREQ"));
	if(Header_check(source, "RESTFRQ"))  Header_set_flt(target, "RESTFRQ",  Header_get_flt(source, "RESTFRQ"));
	if(Header_check(source, "VELREF"))  Header_set_int(target, "VELREF",  Header_get_int(source, "VELREF"));
	if(Header_check(source, "SPECSYS"))
	{
		Header_get_str(source, "SPECSYS", value);
		Header_set_str(target, "SPECSYS", value);
	}
	
	// Equinox and coordinate system
	if(Header_check(source, "EQUINOX"))  Header_set_flt(target, "EQUINOX",  Header_get_flt(source, "EQUINOX"));
	if(Header_check(source, "EPOCH"))    Header_set_flt(target, "EPOCH",    Header_get_flt(source, "EPOCH"));
	if(Header_check(source, "RADESYS"))
	{
		Header_get_str(source, "RADESYS", value);
		Header_set_str(target, "RADESYS", value);
	}
	
	// Longitude and latitude of celestial pole
	if(Header_check(source, "LONPOLE")) Header_set_flt(target, "LONPOLE", Header_get_flt(source, "LONPOLE"));
	if(Header_check(source, "LATPOLE")) Header_set_flt(target, "LATPOLE", Header_get_flt(source, "LATPOLE"));
	
	return;
}



/// @brief Copy miscellaneous information from one header to another
///
/// Public method for copying miscellaneous header information from
/// one data cube to another. This method is intended to be used if
/// information about the flux unit (keyword: `BUNIT`) or the beam
/// (`BMAJ`, `BMIN`, `BPA`) needs to be copied from one cube to another.
///
/// @param source      Data cube from which to copy WCS information.
/// @param target      Data cube to which to copy WCS information.
/// @param copy_bunit  Should the `BUNIT` keyword be copied?
/// @param copy_beam   Should beam information be copied?

PUBLIC void Header_copy_misc(const Header *source, Header *target, const bool copy_bunit, const bool copy_beam)
{
	// Sanity checks
	check_null(source);
	check_null(target);
	check_null(source->header);
	check_null(target->header);
	
	if(copy_bunit && Header_check(source, "BUNIT"))
	{
		char value[FITS_HEADER_VALUE_SIZE + 1];
		Header_get_str(source, "BUNIT", value);
		Header_set_str(target, "BUNIT", value);
	}
	
	if(copy_beam)
	{
		if(Header_check(source, "BMAJ")) Header_set_flt(target, "BMAJ", Header_get_flt(source, "BMAJ"));
		if(Header_check(source, "BMIN")) Header_set_flt(target, "BMIN", Header_get_flt(source, "BMIN"));
		if(Header_check(source, "BPA"))  Header_set_flt(target, "BPA",  Header_get_flt(source, "BPA"));
	}
	
	return;
}



/// @brief Adjust WCS information to subregion
///
/// Public method for adjusting the `NAXIS`i and `CRPIX`i keywords in
/// the header to account for a subregion of the specified dimensions.
/// This is useful if a cut-out of a larger cube has been created, but
/// the header is still that of the original cube.
///
/// @param  self   Object self-reference.
/// @param  x_min  Lower limit of new region in `x`.
/// @param  x_max  Upper limit of new region in `x`.
/// @param  y_min  Lower limit of new region in `y`.
/// @param  y_max  Upper limit of new region in `y`.
/// @param  z_min  Lower limit of new region in `z`.
/// @param  z_max  Upper limit of new region in `z`.

PUBLIC void Header_adjust_wcs_to_subregion(Header *self, const size_t x_min, const size_t x_max, const size_t y_min, const size_t y_max, const size_t z_min, const size_t z_max)
{
	const size_t nx = x_max - x_min + 1;
	const size_t ny = y_max - y_min + 1;
	const size_t nz = z_max - z_min + 1;
	
	if(Header_check(self, "NAXIS1")) Header_set_int(self, "NAXIS1", nx);
	if(Header_check(self, "NAXIS2")) Header_set_int(self, "NAXIS2", ny);
	if(Header_check(self, "NAXIS3")) Header_set_int(self, "NAXIS3", nz);
	if(Header_check(self, "CRPIX1")) Header_set_flt(self, "CRPIX1", Header_get_flt(self, "CRPIX1") - x_min);
	if(Header_check(self, "CRPIX2")) Header_set_flt(self, "CRPIX2", Header_get_flt(self, "CRPIX2") - y_min);
	if(Header_check(self, "CRPIX3")) Header_set_flt(self, "CRPIX3", Header_get_flt(self, "CRPIX3") - z_min);
	
	return;
}
