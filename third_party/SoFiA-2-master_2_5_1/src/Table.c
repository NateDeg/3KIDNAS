// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Table.c) - Source Finding Application                   //
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

/// @file   Table.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Container class for reading tabulated data from text files.


#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Table.h"



/// @brief Container class for reading tabulated data from text files.
///
/// This class provides a convenient container for reading tabulated
/// data from a text file. The container size will be automatically
/// adjusted to match the number of rows and columns in the file. The
/// data values will be stored as double-precision floating-point.
/// The standard constructor has been set to private to disallow the
/// creation of empty Table objects.

CLASS Table
{
	size_t cols;   ///< Number of table columns.
	size_t rows;   ///< Number of table rows.
	double *data;  ///< Pointer to array of data values.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new and empty Table object.
/// Note that the destructor will have to be called explicitly once
/// the object is no longer required to release its memory again.
///
/// @return Pointer to newly created Table object.
///
/// @note The standard constructor has been made private to disable
///       the creation of empty Table objects, as there are no methods
///       to change the size of a Table object. The Table_from_file()
///       constructor must instead be used to directly generate Tables
///       from tabulated data stored in a text file.

PRIVATE Table *Table_new(void)
{
	Table *self = (Table *)memory(MALLOC, 1, sizeof(Table));
	
	self->cols = 0;
	self->rows = 0;
	self->data = NULL;
	
	return self;
}



/// @brief Alternative constructor for creating tables from text file
///
/// Alternative constructor. Will create a new Table object from
/// tabulated data read from the specified text file. The list of
/// delimiters to be used to separate data columns can be specified
/// using the `delimiters` string (e.g. "` \t`" will use space and
/// tab characters as delimiters). Consecutive delimiting characters
/// will always be merged even when of mixed type. A new object
/// filled with the data from the file will be returned.
///
/// If no valid data are found in the specified file, then an empty
/// Table object with 0 rows and columns will be returned. The size
/// of the table will be defined by the first data row in the file.
/// All subsequent rows must have the same number of columns. If
/// fewer columns are encountered, the method will terminate with
/// an error message. If more columns are encountered, then any
/// excess columns will be silently ignored. If an entry is
/// encountered that does not constitute a floating-point number,
/// then the resulting table entry will be set to 0.
///
/// The data file can contain empty lines and lines beginning with
/// a comment character (`#`); these will be ignored.
///
/// @param filename    The path to the file from which to read data.
/// @param delimiters  List of delimiting characters to be used to
///                    separate columns in the data file.
///
/// @return Pointer to newly created Table object.

PUBLIC Table *Table_from_file(const char *filename, const char *delimiters)
{
	// Sanity checks
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "Empty file name provided.");
	check_null(delimiters);
	
	// Allocate memory for a single line
	char *line = (char *)memory(MALLOC, TABLE_MAX_LINE_SIZE, sizeof(char));
	size_t counter = 0;
	
	// Try to open file
	FILE *fp = fopen(filename, "r");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open input file: %s.", filename);
	
	// Read beginning of file to establish number of columns
	while(fgets(line, TABLE_MAX_LINE_SIZE, fp))
	{
		// Trim line and check for comments and empty lines
		char *trimmed = trim_string(line);
		if(strlen(trimmed) == 0 || !isalnum(trimmed[0])) continue;
		
		char *entry = strtok(line, delimiters);
		while(entry != NULL)
		{
			++counter;
			entry = strtok(NULL, delimiters);
		}
		
		if(counter > 0) break;
	}
	
	// Return empty table if no valid data found
	if(counter == 0)
	{
		warning("No valid data found in file %s.\n         Returning empty table.", filename);
		fclose(fp);
		return Table_new();
	}
	
	// Otherwise reset file pointer
	rewind(fp);
	
	// Create a new table
	Table *self = Table_new();
	self->cols  = counter;
	
	while(fgets(line, TABLE_MAX_LINE_SIZE, fp))
	{
		// Trim line and check for comments and empty lines
		char *trimmed = trim_string(line);
		if(strlen(trimmed) == 0 || !isalnum(trimmed[0])) continue;
		
		// Add row
		self->rows += 1;
		self->data = memory_realloc(self->data, self->rows * self->cols, sizeof(double));
		
		char *entry = strtok(trimmed, delimiters);
		for(size_t i = 0; i < self->cols; ++i)
		{
			ensure(entry != NULL, ERR_USER_INPUT, "Inconsistent number of data columns in file %s.\n       %zu columns expected, but only %zu columns found in data row %zu.", filename, self->cols, i, self->rows);
			self->data[TABLE_INDEX(self->rows - 1, i)] = strtod(entry, NULL);
			entry = strtok(NULL, delimiters);
		}
	}
	
	// Close file again and clean up
	fclose(fp);
	free(line);
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release all
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Table_delete(Table *self)
{
	if(self != NULL) free(self->data);
	free(self);
	return;
}



/// @brief Return the number of table rows
///
/// Public methods for retrieving the number of rows in
/// the specified Table object.
///
/// @param self  Object self-reference.
///
/// @return Number of table rows.

PUBLIC size_t Table_rows(const Table *self)
{
	return self == NULL ? 0 : self->rows;
}



/// @brief Return the number of table columns
///
/// Public methods for retrieving the number of columns in
/// the specified Table object.
///
/// @param self  Object self-reference.
///
/// @return Number of table columns.

PUBLIC size_t Table_cols(const Table *self)
{
	return self == NULL ? 0 : self->cols;
}



/// @brief Get table value at specified row and column
///
/// Public method for retrieving the table entry at the specified
/// row and column of the Table object. The method will terminate
/// with an error message if the requested row or column is out of
/// range.
///
/// @param self  Object self-reference.
/// @param row   Requested row.
/// @param col   Requested column.
///
/// @return Table value at specified row and column.

PUBLIC double Table_get(const Table *self, const size_t row, const size_t col)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Requested Table column or row out of range.");
	
	return self->data[TABLE_INDEX(row, col)];
}



/// @brief Set table value at specified row and column
///
/// Public method for writing the specified value into the specified
/// row and column of the Table object, thus overwriting any
/// existing value. The method will terminate with an error message
/// if the requested row or column is out of range.
///
/// @param self   Object self-reference.
/// @param row    Requested row.
/// @param col    Requested column.
/// @param value  Value to write into the table.

PUBLIC void Table_set(Table *self, const size_t row, const size_t col, const double value)
{
	// Sanity checks
	check_null(self);
	ensure(row < self->rows && col < self->cols, ERR_INDEX_RANGE, "Requested Table column or row out of range.");
	
	self->data[TABLE_INDEX(row, col)] = value;
	return;
}
