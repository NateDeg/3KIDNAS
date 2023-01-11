// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Source.c) - Source Finding Application                  //
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

/// @file   Source.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Class for storing and handling source parameters.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Source.h"
#include "String.h"

/// @brief Make SourceValue union a type.
typedef union SourceValue SourceValue;

/// @brief Union for storing either integer or floating-point parameter value.
union SourceValue
{
	double value_flt;  ///< Floating-point parameter value.
	size_t value_int;  ///< Integer parameter value.
};



/// @brief Class for storing and handling source parameters
///
/// The purpose of this class is to provide a structure for storing
/// and handling the measured parameters of a source. Parameters are
/// composed of a name, value, type and unit. Both long integer and
/// double-precision floating-point values are supported. In addition,
/// a source can be assigned an identifier in the form of a string,
/// e.g. a source name.

CLASS Source
{
	String         *identifier;  ///< Source name.
	size_t          n_par;       ///< Total number of source parameters stored.
	SourceValue    *values;      ///< Array of source parameter values.
	unsigned char  *types;       ///< Array of source parameter data types. Can be SOURCE_TYPE_INT or SOURCE_TYPE_FLT.
	String        **names;       ///< Array of source parameter names.
	String        **units;       ///< Array of source parameter units.
	String        **ucds;        ///< Array of Unified Content Descriptors (UCDs).
	int             verbosity;   ///< Verbosity level (0 or 1).
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new and empty Source object
/// and return a pointer to the newly created object. No memory
/// will be allocated other than for the object itself. Note that
/// the destructor will need to be called explicitly once the object
/// is no longer required to release any memory allocated during the
/// the lifetime of the object.
///
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created Source object.

PUBLIC Source *Source_new(const bool verbosity)
{
	// Allocate memory for new source
	Source *self = (Source *)memory(MALLOC, 1, sizeof(Source));
	
	// Initialise properties
	self->identifier = String_new("");
	self->n_par      = 0;
	self->values     = NULL;
	self->types      = NULL;
	self->names      = NULL;
	self->units      = NULL;
	self->ucds       = NULL;
	
	self->verbosity  = verbosity;
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the me-
/// mory occupied by the object.
///
/// @param self  Object self-reference.
///
/// @note If the source is added to a Catalog object, the
///       destructor must not be called manually as the Catalog
///       object will take ownership of the source and call its
///       destructor at the end of its lifetime.

PUBLIC void Source_delete(Source *self)
{
	if(self != NULL)
	{
		for(size_t i = self->n_par; i--;)
		{
			if(self->names != NULL) String_delete(self->names[i]);
			if(self->units != NULL) String_delete(self->units[i]);
			if(self->ucds  != NULL) String_delete(self->ucds[i]);
		}
		
		String_delete(self->identifier);
		
		free(self->values);
		free(self->types);
		free(self->names);
		free(self->units);
		free(self->ucds);
		
		free(self);
	}
	
	return;
}



/// @brief Set source identifier
///
/// Public method for setting the identifier of a source. Note that
/// names are case-sensitive.
///
/// @param self  Object self-reference.
/// @param name  Name to be used as identifier.

PUBLIC void Source_set_identifier(Source *self, const char *name)
{
	check_null(self);
	String_set(self->identifier, name);
	return;
}



/// @brief Append a new parameter of floating-point type
///
/// Public method for appending a new parameter of floating-point
/// type to an existing source. The value will be stored as a
/// double-precision floating-point number. Note that the function
/// does not check if the specified parameter name already exists
/// in the current parameter list; if it did, a new parameter with
/// the same name would be appended at the end. Note that names are
/// case-sensitive.
///
/// @param self   Object self-reference.
/// @param name   Name of the new parameter.
/// @param value  Value of the new parameter.
/// @param unit   Unit of the new parameter.
/// @param ucd    Unified Content Descriptor for the new parameter.

PUBLIC void Source_add_par_flt(Source *self, const char *name, const double value, const char *unit, const char *ucd)
{
	// Sanity checks
	check_null(self);
	check_null(name);
	check_null(unit);
	check_null(ucd);
	
	// Reserve memory for one additional parameter
	Source_append_memory(self);
	
	// Copy new parameter information
	self->values[self->n_par - 1].value_flt = value;
	self->types[self->n_par - 1] = SOURCE_TYPE_FLT;
	String_set(self->names[self->n_par - 1], name);
	String_set(self->units[self->n_par - 1], unit);
	String_set(self->ucds [self->n_par - 1],  ucd);
	
	return;
}



/// @brief Append a new parameter of integer type
///
/// Public method for appending a new parameter of integer type to
/// an existing source. The value will be stored as a long integer
/// number. Note that the function does not check if the specified
/// parameter name already exists in the current parameter list;
/// if it did, a new parameter with the same name would be
/// appended at the end. Note that names are case-sensitive.
///
/// @param self   Object self-reference.
/// @param name   Name of the new parameter.
/// @param value  Value of the new parameter.
/// @param unit   Unit of the new parameter.
/// @param ucd    Unified Content Descriptor for the new parameter.

PUBLIC void Source_add_par_int(Source *self, const char *name, const long int value, const char *unit, const char *ucd)
{
	// Sanity checks
	check_null(self);
	check_null(name);
	check_null(unit);
	check_null(ucd);
	
	// Reserve memory for one additional parameter
	Source_append_memory(self);
	
	// Copy new parameter information
	self->values[self->n_par - 1].value_int = value;
	self->types[self->n_par - 1] = SOURCE_TYPE_INT;
	String_set(self->names[self->n_par - 1], name);
	String_set(self->units[self->n_par - 1], unit);
	String_set(self->ucds [self->n_par - 1],  ucd);
	
	return;
}



/// @brief Set source parameter as floating-point value
///
/// Public method for setting the specified parameter of the specified
/// source as a double-precision floating-point number. If a parameter
/// of the same name already exists, its value will be replaced;
/// otherwise, a new parameter will be created and added to the
/// existing parameter list. Note that names are case-sensitive.
///
/// If the parameter already exists and only its value should be
/// updated, then `unit` and `ucd` can be set to `NULL`. This would
/// cause an error, however, if the parameter did not yet exist.
///
/// @param self   Object self-reference.
/// @param name   Name of the parameter to be set.
/// @param value  Value of the parameter to be set.
/// @param unit   Unit of the parameter to be set.
/// @param ucd    Unified Content Descriptor for the parameter.

PUBLIC void Source_set_par_flt(Source *self, const char *name, const double value, const char *unit, const char *ucd)
{
	// Sanity checks
	check_null(self);
	check_null(name);
	
	// Check if parameter of same name already exists
	size_t index = 0;
	
	if(Source_par_exists(self, name, &index))
	{
		// If so, overwrite with new parameter information
		self->values[index].value_flt = value;
		self->types[index] = SOURCE_TYPE_FLT;
		String_set(self->names[index], name);
		if(unit != NULL) String_set(self->units[index], unit);
		if(ucd  != NULL) String_set(self->ucds [index],  ucd);
	}
	else
	{
		// Otherwise add as new parameter
		check_null(unit);
		check_null(ucd);
		Source_add_par_flt(self, name, value, unit, ucd);
	}
	
	return;
}



/// @brief Set source parameter as integer value
///
/// Public method for setting the specified parameter of the specified
/// source as a long integer number. If a parameter of the same name
/// already exists, its value will be replaced; otherwise, a new
/// parameter will be created and added to the existing parameter list.
/// Note that names are case-sensitive.
///
/// If the parameter already exists and only its value should be
/// updated, then `unit` and `ucd` can be set to `NULL`. This would
/// cause an error, however, if the parameter did not yet exist.
///
/// @param self   Object self-reference.
/// @param name   Name of the parameter to be set.
/// @param value  Value of the parameter to be set.
/// @param unit   Unit of the parameter to be set.
/// @param ucd    Unified Content Descriptor for the parameter.

PUBLIC void Source_set_par_int(Source *self, const char *name, const long int value, const char *unit, const char *ucd)
{
	// Sanity checks
	check_null(self);
	check_null(name);
	
	// Check if parameter already exists
	size_t index = 0;
	
	if(Source_par_exists(self, name, &index))
	{
		// If so, overwrite with new parameter information
		self->values[index].value_int = value;
		self->types[index] = SOURCE_TYPE_INT;
		String_set(self->names[index], name);
		if(unit != NULL) String_set(self->units[index], unit);
		if(ucd  != NULL) String_set(self->ucds [index], ucd);
	}
	else
	{
		// Otherwise add as new parameter
		check_null(unit);
		check_null(ucd);
		Source_add_par_int(self, name, value, unit, ucd);
	}
	
	return;
}



/// @brief Extract a parameter as floating-point value
///
/// Public method for extracting the value of the specified parameter
/// from the specified source as a double-precision floating-point
/// number.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter to be extracted.
///
/// @return Requested parameter value as type `double`.

PUBLIC double Source_get_par_flt(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return self->values[index].value_flt;
}



/// @brief Extract a parameter as integer value
///
/// Public method for extracting the value of the specified parameter
/// from the specified source as a long integer number.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter to be extracted.
///
/// @return Requested parameter value as type `long int`.

PUBLIC long int Source_get_par_int(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return self->values[index].value_int;
}



/// @brief Extract a parameter by name as floating-point value
///
/// Public method for extracting the value of the specified parameter
/// from the specified source as a double-precision floating-point
/// number. If a parameter of the same name does not exist, a value
/// of `NaN` will instead be returned.
///
/// @param self  Object self-reference.
/// @param name  Name of the parameter to be extracted.
///
/// @return Requested parameter value as type `double`.

PUBLIC double Source_get_par_by_name_flt(const Source *self, const char *name)
{
	check_null(self);
	check_null(name);
	for(size_t i = self->n_par; i--;) if(String_compare(self->names[i], name)) return self->values[i].value_flt;
	warning_verb(self->verbosity, "Parameter \'%s\' not found.", name);
	return NAN;
}



/// @brief Extract a parameter by name as integer value
///
/// Public method for extracting the value of the specified parameter
/// from the specified source as a signed long integer value. If a
/// parameter of the same name does not exist, a value of 0 will
/// instead be returned.
///
/// @param self  Object self-reference.
/// @param name  Name of the parameter to be extracted.
///
/// @return Requested parameter value as type `long int`.

PUBLIC long int Source_get_par_by_name_int(const Source *self, const char *name)
{
	check_null(self);
	check_null(name);
	for(size_t i = self->n_par; i--;) if(String_compare(self->names[i], name)) return self->values[i].value_int;
	warning_verb(self->verbosity, "Parameter \'%s\' not found.", name);
	return 0;
}



/// @brief Add a position offset to `x`, `y`, `z` parameters
///
/// Public method for adding a position offset to parameters named
/// `x`, `y`, `z`, `x_min`, `x_max`, `y_min`, `y_max`,
/// `z_min`, and `z_max`. Only existing parameters will be shifted
/// and non-existing ones ignored. Offsets can only be positive, as
/// negative pixel coordinates are not possible.
///
/// @param self  Object self-reference.
/// @param dx    Position offset in x.
/// @param dy    Position offset in y.
/// @param dz    Position offset in z.

PUBLIC void Source_offset_xyz(Source *self, const size_t dx, const size_t dy, const size_t dz)
{
	check_null(self);
	
	size_t index = 0;
	
	if(Source_par_exists(self, "x", &index)) self->values[index].value_flt += (double)dx;
	if(Source_par_exists(self, "y", &index)) self->values[index].value_flt += (double)dy;
	if(Source_par_exists(self, "z", &index)) self->values[index].value_flt += (double)dz;
	
	if(Source_par_exists(self, "x_min", &index)) self->values[index].value_int += dx;
	if(Source_par_exists(self, "x_max", &index)) self->values[index].value_int += dx;
	if(Source_par_exists(self, "y_min", &index)) self->values[index].value_int += dy;
	if(Source_par_exists(self, "y_max", &index)) self->values[index].value_int += dy;
	if(Source_par_exists(self, "z_min", &index)) self->values[index].value_int += dz;
	if(Source_par_exists(self, "z_max", &index)) self->values[index].value_int += dz;
	
	return;
}



/// @brief Check if source parameter exists
///
/// Public method for checking if a parameter of the specified name
/// already exists in the specified source. The function will return
/// `true` if the parameter exists and `false` otherwise. Note
/// that name is case-sensitive. The variable `index` will be set
/// to the index of the parameter if found. Otherwise it will be
/// left untouched. If no index is required, a `NULL` pointer can
/// instead be provided.
///
/// @param self   Object self-reference.
/// @param name   Name of the parameter to be checked.
/// @param index  Pointer to a variable that will hold the index
///               of the source parameter if found.
///
/// @return Returns `true` if the parameter exists, `false` otherwise.

PUBLIC bool Source_par_exists(const Source *self, const char *name, size_t *index)
{
	// Sanity checks
	check_null(self);
	check_null(name);
	
	for(size_t i = self->n_par; --i;)
	{
		if(String_compare(self->names[i], name))
		{
			if(index != NULL) *index = i;
			return true;
		}
	}
	
	return false;
}



/// @brief Extract name of parameter by index
///
/// Public method for returning a pointer to the name string of the
/// specified parameter.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter the unit of which is to
///               be returned.
///
/// @return Pointer to the name string of the specified parameter.

PUBLIC const char *Source_get_name(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return String_get(self->names[index]);
}



/// @brief Extract unit of parameter by index
///
/// Public method for returning a pointer to the unit string of the
/// specified parameter.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter the unit of which is to
///               be returned.
///
/// @return Pointer to the unit string of the specified parameter.

PUBLIC const char *Source_get_unit(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return String_get(self->units[index]);
}



/// @brief Extract type of parameter by index
///
/// Public method for returning the data type of the specified
/// parameter, where 0 means integer and 1 means floating-point.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter the type of which is to
///               be returned.
///
/// @return Data type of the specified parameter.

PUBLIC unsigned char Source_get_type(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return self->types[index];
}



/// @brief Extract UCD of parameter by index
///
/// Public method for returning a pointer to the Unified Content
/// Descriptor (UCD) string of the specified parameter.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter the UCD of which is to
///               be returned.
///
/// @return Pointer to the UCD string of the specified parameter.

PUBLIC const char *Source_get_ucd(const Source *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Source parameter index out of range.");
	return String_get(self->ucds[index]);
}



/// @brief Get identifier of specified source
///
/// Public method for returning the identifier string of the specified
/// source. If no identifier has been set, an empty string will be
/// returned instead.
///
/// @param self  Object self-reference.
///
/// @return Identifier of the specified source.

PUBLIC const char *Source_get_identifier(const Source *self)
{
	check_null(self);
	return String_get(self->identifier);
}



/// @brief Get number of parameters for specified source
///
/// Public method for returning the number of parameters currently
/// defined for the specified source.
///
/// @param self  Object self-reference.
///
/// @return Number of parameters currently defined.

PUBLIC size_t Source_get_num_par(const Source *self)
{
	check_null(self);
	return self->n_par;
}



/// @brief Reallocate memory for one additional parameter
///
/// Private method for allocating additional memory for one more
/// parameter in the specified source. Note that this will not
/// create a new parameter yet, but just allocate the memory needed
/// to append a parameter at the end of the parameter list. The
/// function should be called from public member functions that will
/// add parameters to a source prior to assigning the new parameter
/// values.
///
/// @param self  Object self-reference.

PRIVATE void Source_append_memory(Source *self)
{
	self->n_par += 1;
	self->values = (SourceValue *)   memory_realloc(self->values, self->n_par, sizeof(SourceValue));
	self->types  = (unsigned char *) memory_realloc(self->types,  self->n_par, sizeof(unsigned char));
	self->names  = (String **)       memory_realloc(self->names,  self->n_par, sizeof(String *));
	self->units  = (String **)       memory_realloc(self->units,  self->n_par, sizeof(String *));
	self->ucds   = (String **)       memory_realloc(self->ucds,   self->n_par, sizeof(String *));
	
	// Call constructor on all strings
	self->names[self->n_par - 1] = String_new("");
	self->units[self->n_par - 1] = String_new("");
	self->ucds [self->n_par - 1] = String_new("");
	
	return;
}
