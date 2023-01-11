// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (Parameter.c) - Source Finding Application               //
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

/// @file   Parameter.c
/// @author Tobias Westmeier
/// @date   23/11/2021
/// @brief  Class for storing and handling SoFiA parameter settings.


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "Parameter.h"
#include "String.h"



/// @brief Class for storing and handling SoFiA parameter settings
///
/// The purpose of this class is to provide a structure for handling
/// SoFiA parameter settings. Settings can be loaded from a file and
/// then read or updated as needed. All parameter settings are treated
/// as strings, and several methods are available for extracting the
/// parameter value as a specific data type.

CLASS Parameter
{
	size_t   n_par;      ///< Number of parameter settings.
	String **keys;       ///< Pointer to array of String objects holding the parameter names.
	String **values;     ///< Pointer to array of String objects holding the parameter values.
	int      verbosity;  ///< Verbosity level.
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new and empty Parameter object
/// and return a pointer to the newly created object. No memory will
/// be allocated other than for the object itself. Note that the
/// destructor will need to be called explicitly once the object is no
/// longer required to release any memory allocated during the
/// lifetime of the object.
///
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created Parameter object.

PUBLIC Parameter *Parameter_new(const bool verbosity)
{
	Parameter *self = (Parameter *)memory(MALLOC, 1, sizeof(Parameter));
	
	self->n_par  = 0;
	self->keys   = NULL;
	self->values = NULL;
	
	self->verbosity = verbosity;
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the memory
/// occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void Parameter_delete(Parameter *self)
{
	if(self != NULL)
	{
		if(self->keys   != NULL) for(size_t i = self->n_par; i--;) String_delete(self->keys[i]);
		if(self->values != NULL) for(size_t i = self->n_par; i--;) String_delete(self->values[i]);
		
		free(self->keys);
		free(self->values);
		free(self);
	}
	
	return;
}



/// @brief Return the number of parameters in the list
///
/// Public method for returning the number of parameters currently
/// stored in the parameter list object.
///
/// @param self  Object self-reference.
///
/// @return Number of currently stored parameters.

PUBLIC size_t Parameter_get_size(const Parameter *self)
{
	check_null(self);
	return self->n_par;
}



/// @brief Set parameter to given value
///
/// Public method for setting a parameter to a specified value. If
/// a parameter of that name already exists, its previous value will
/// be overwritten. Otherwise, a new parameter will be appended to
/// the current parameter list.
///
/// @param self   Object self-reference.
/// @param key    Name of the parameter to be set.
/// @param value  String containing the value of the parameter.

PUBLIC void Parameter_set(Parameter *self, const char *key, const char *value)
{
	// Sanity checks
	check_null(self);
	check_null(key);
	check_null(value);
	
	size_t index;
	
	// Check if parameter already exists
	if(Parameter_exists(self, key, &index))
	{
		warning_verb(self->verbosity, "Parameter \'%s\' already exists.\n         Replacing existing definition.", key);
	}
	else
	{
		Parameter_append_memory(self);
		index = self->n_par - 1;
		String_set(self->keys[index], key);
	}
	
	String_set(self->values[index], value);
	
	return;
}



/// @brief Check if parameter exists
///
/// Public method for checking if the specified parameter name
/// already exists in the current parameter list. If the parameter
/// is found, the function returns`true`, otherwise `false`. The
/// variable `index` will be set to the index of the parameter if
/// found and otherwise be left untouched. If the index is not
/// required, a `NULL` pointer can instead be provided.
///
/// @param self   Object self-reference.
/// @param key    Name of the parameter to be checked.
/// @param index  Pointer to an index variable that will be set to
///               the index of the parameter, if found. This can
///               be `NULL` if the index is not required.
///
/// @return Returns `true` if the specified parameter exists and
///         `false` otherwise.

PUBLIC bool Parameter_exists(const Parameter *self, const char *key, size_t *index)
{
	// Sanity checks
	check_null(self);
	check_null(key);
	ensure(strlen(key), ERR_USER_INPUT, "Empty parameter keyword provided.");
	
	for(size_t i = self->n_par; i--;)
	{
		if(String_compare(self->keys[i], key))
		{
			if(index != NULL) *index = i;
			return true;
		}
	}
	
	return false;
}



/// @brief Extract parameter value as floating-point number
///
/// Public method for returning the value of the specified parameter
/// as a double-precision floating-point number. If the parameter
/// does not exist, a value of `NaN` will instead be returned.
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns the value of the specified parameter.

PUBLIC double Parameter_get_flt(const Parameter *self, const char *key)
{
	const char *value_raw = Parameter_get_raw(self, key);
	if(value_raw == NULL) return NAN;
	return strtod(value_raw, NULL);
}



/// @brief Extract parameter value as integer number
///
/// Public method for returning the value of the specified parameter
/// as a long integer number. If the parameter does not exist, a
/// value of 0 will instead be returned.
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns the value of the specified parameter.

PUBLIC long int Parameter_get_int(const Parameter *self, const char *key)
{
	const char *value_raw = Parameter_get_raw(self, key);
	if(value_raw == NULL) return 0L;
	return strtol(value_raw, NULL, 10);
}



/// @brief Extract parameter value as unsigned integer number
///
/// Public method for returning the value of the specified parameter
/// as an unsigned integer number. If the parameter does not exist,
/// a value of 0 will instead be returned.
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns the value of the specified parameter.

PUBLIC unsigned long int Parameter_get_uint(const Parameter *self, const char *key)
{
	const char *value_raw = Parameter_get_raw(self, key);
	if(value_raw == NULL) return 0UL;
	return (unsigned long int)strtol(value_raw, NULL, 10);
}



/// @brief Extract parameter value as Boolean value
///
/// Public method for returning the value of the specified parameter
/// as a Boolean value (true or false). If the parameter does not
/// exist, a value of `false` will be returned.
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns the value of the specified parameter.

PUBLIC bool Parameter_get_bool(const Parameter *self, const char *key)
{
	const char *value_raw = Parameter_get_raw(self, key);
	if(value_raw == NULL) return false;
	return strcmp(value_raw, "true") == 0;
}



/// @brief Extract parameter value as String
///
/// Public method for returning the value of the specified parameter
/// as a pointer to a string value. If the parameter does not exist,
/// a `NULL` pointer will be returned. Note that this is just a
/// convenience function that will simply return the output of
/// Parameter_get_raw().
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns a pointer to the string representing the value.

PUBLIC const char *Parameter_get_str(const Parameter *self, const char *key)
{
	return Parameter_get_raw(self, key);
}



/// @brief Extract parameter value as String by index
///
/// Public method for returning the value of the parameter specified
/// by `index` as a pointer to a string value. If the index is out
/// of range, the process will be terminated.
///
/// @param self   Object self-reference.
/// @param index  Index of the parameter to be extracted.
///
/// @return Returns a pointer to the string representing the value.

PUBLIC const char *Parameter_get_str_index(const Parameter *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Parameter list index out of range.");
	return String_get(self->values[index]);
}



/// @brief Extract key at specified index position
///
/// Public method for returning the name of the parameter at the
/// specified index position as a pointer to a string value. The
/// process will be terminated if the index is out of range.
///
/// @param self   Object self-reference.
/// @param index  Index of the key to be returned.
///
/// @return Returns a pointer to the string representing the key.

PUBLIC const char *Parameter_get_key(const Parameter *self, const size_t index)
{
	check_null(self);
	ensure(index < self->n_par, ERR_INDEX_RANGE, "Parameter list index out of range.");
	return String_get(self->keys[index]);
}



/// @brief Extract parameter value as raw string
///
/// Private method for returning the value of the specified parameter
/// as a pointer to the raw value string. If the parameter does not
/// exist, a `NULL` pointer will instead be returned.
///
/// @param self  Object self-reference.
/// @param key   Name of the parameter to be extracted.
///
/// @return Returns the value of the specified parameter as a raw string.

PRIVATE const char *Parameter_get_raw(const Parameter *self, const char *key)
{
	size_t index;
	if(Parameter_exists(self, key, &index)) return String_get(self->values[index]);
	return NULL;
}



/// @brief Load parameter settings from file
///
/// Public method for loading parameter settings from an input file
/// of the name `filename`. The parameter settings in the file must
/// be of the form
///
/// \code
///   key = [value] [# comment]
/// \endcode
///
/// where both the value and comment are optional (indicated by the
/// brackets). Empty lines and lines starting with # (comments) will be
/// ignored. If mode = `PARAMETER_APPEND`, the function will update any
/// existing parameter or append a new parameter setting if the parameter
/// name does not yet exist. If mode = `PARAMETER_UPDATE`, the function
/// will only update existing parameters and discard any non-existing
/// parameter settings. If the parameter `pipeline.pedantic` is found
/// and set to true, an error message will be produced if a non-existing
/// parameter name is found, and the process will be terminated in this
/// case (only for mode = `PARAMETER_UPDATE`).
///
/// @param self      Object self-reference.
/// @param filename  Name of the input file.
/// @param mode      Mode of operation; can be `PARAMETER_APPEND` or
///                  `PARAMETER_UPDATE`.

PUBLIC void Parameter_load(Parameter *self, const char *filename, const int mode)
{
	// Sanity checks
	check_null(self);
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "Empty file name provided.");
	ensure(mode == PARAMETER_APPEND
		|| mode == PARAMETER_UPDATE, ERR_USER_INPUT, "Mode must be \'PARAMETER_APPEND\' or \'PARAMETER_UPDATE\'.");
	
	// Try to open file
	FILE *fp = fopen(filename, "r");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open input file: %s.", filename);
	
	// Allocate memory for a single line
	char *line = (char *)memory(MALLOC, PARAMETER_MAX_LINE_SIZE, sizeof(char));
	
	// Record if unknown parameters are encountered
	bool unknown_parameter = false;
	
	// Read lines from file
	while(fgets(line, PARAMETER_MAX_LINE_SIZE, fp))
	{
		// Trim line and check for comments and empty lines
		char *trimmed = trim_string(line);
		if(strlen(trimmed) == 0 || !isalnum(trimmed[0])) continue;
		
		// Extract keyword
		char *key = trim_string(strtok(trimmed, "="));
		if(key == NULL || strlen(key) == 0)
		{
			warning_verb(self->verbosity, "Failed to parse the following setting:\n         %s", trimmed);
			continue;
		}
		
		// Check if keyword already exists
		if(mode == PARAMETER_UPDATE && !Parameter_exists(self, key, NULL))
		{
			message("  Unknown parameter: \'%s\'", key);
			unknown_parameter = true;
		}
		else
		{
			// Extract value and insert into parameter list
			char *value = trim_string(strtok(NULL, "#"));
			if(value == NULL || strlen(value) == 0)
			{
				warning_verb(self->verbosity, "Parameter \'%s\' has no value.", key);
				Parameter_set(self, key, "");
			}
			else
			{
				Parameter_set(self, key, value);
			}
			
			// Check for verbosity keyword
			if(strcmp(key, "pipeline.verbose") == 0) self->verbosity = Parameter_get_bool(self, "pipeline.verbose");
		}
	}
	
	// De-allocate memory and close file
	free(line);
	fclose(fp);
	
	// Check pedantic keyword
	ensure(!unknown_parameter || !Parameter_get_bool(self, "pipeline.pedantic"), ERR_USER_INPUT, "Unknown parameter settings encountered. Please check\n       your input or change \'pipeline.pedantic\' to \'false\'.");
	
	return;
}



/// @brief Set SoFiA default parameters
///
/// Public method for setting the SoFiA default parameters. Parameters
/// that do not yet exist will be created and any existing ones reset
/// to their default value, so this method can either be used on a new
/// Parameter object to establish default settings or on an already
/// populated Parameter object to reset all values to their default.
///
/// @param self  Object self-reference.

PUBLIC void Parameter_default(Parameter *self)
{
	// Sanity checks
	check_null(self);
	
	// Global settings
	Parameter_set(self, "pipeline.verbose"         , "false");
	Parameter_set(self, "pipeline.pedantic"        , "true");
	Parameter_set(self, "pipeline.threads"         , "0");
	
	// Input
	Parameter_set(self, "input.data"               , "");
	Parameter_set(self, "input.region"             , "");
	Parameter_set(self, "input.gain"               , "");
	Parameter_set(self, "input.noise"              , "");
	Parameter_set(self, "input.weights"            , "");
	Parameter_set(self, "input.mask"               , "");
	Parameter_set(self, "input.invert"             , "false");
	
	// Flagging
	Parameter_set(self, "flag.region"              , "");
	Parameter_set(self, "flag.catalog"             , "");
	Parameter_set(self, "flag.radius"              , "5");
	Parameter_set(self, "flag.auto"                , "false");
	Parameter_set(self, "flag.threshold"           , "5.0");
	Parameter_set(self, "flag.log"                 , "false");
	
	// Continuum subtraction
	Parameter_set(self, "contsub.enable"           , "false");
	Parameter_set(self, "contsub.order"            , "0");
	Parameter_set(self, "contsub.threshold"        , "2.0");
	Parameter_set(self, "contsub.shift"            , "4");
	Parameter_set(self, "contsub.padding"          , "3");
	
	// Noise scaling
	Parameter_set(self, "scaleNoise.enable"        , "false");
	Parameter_set(self, "scaleNoise.mode"          , "spectral");
	Parameter_set(self, "scaleNoise.statistic"     , "mad");
	Parameter_set(self, "scaleNoise.fluxRange"     , "negative");
	Parameter_set(self, "scaleNoise.windowXY"      , "25");
	Parameter_set(self, "scaleNoise.windowZ"       , "15");
	Parameter_set(self, "scaleNoise.gridXY"        , "0");
	Parameter_set(self, "scaleNoise.gridZ"         , "0");
	Parameter_set(self, "scaleNoise.interpolate"   , "false");
	Parameter_set(self, "scaleNoise.scfind"        , "false");
	
	// Ripple filter
	Parameter_set(self, "rippleFilter.enable"      , "false");
	Parameter_set(self, "rippleFilter.statistic"   , "median");
	Parameter_set(self, "rippleFilter.windowXY"    , "31");
	Parameter_set(self, "rippleFilter.windowZ"     , "15");
	Parameter_set(self, "rippleFilter.gridXY"      , "0");
	Parameter_set(self, "rippleFilter.gridZ"       , "0");
	Parameter_set(self, "rippleFilter.interpolate" , "false");
	
	// S+C finder
	Parameter_set(self, "scfind.enable"            , "true");
	Parameter_set(self, "scfind.kernelsXY"         , "0, 3, 6");
	Parameter_set(self, "scfind.kernelsZ"          , "0, 3, 7, 15");
	Parameter_set(self, "scfind.threshold"         , "5.0");
	Parameter_set(self, "scfind.replacement"       , "2.0");
	Parameter_set(self, "scfind.statistic"         , "mad");
	Parameter_set(self, "scfind.fluxRange"         , "negative");
	
	// Threshold finder
	Parameter_set(self, "threshold.enable"         , "false");
	Parameter_set(self, "threshold.threshold"      , "5.0");
	Parameter_set(self, "threshold.mode"           , "relative");
	Parameter_set(self, "threshold.statistic"      , "mad");
	Parameter_set(self, "threshold.fluxRange"      , "negative");
	
	// Linker
	Parameter_set(self, "linker.enable"            , "true");
	Parameter_set(self, "linker.radiusXY"          , "1");
	Parameter_set(self, "linker.radiusZ"           , "1");
	Parameter_set(self, "linker.minSizeXY"         , "5");
	Parameter_set(self, "linker.minSizeZ"          , "5");
	Parameter_set(self, "linker.maxSizeXY"         , "0");
	Parameter_set(self, "linker.maxSizeZ"          , "0");
	Parameter_set(self, "linker.minPixels"         , "0");
	Parameter_set(self, "linker.maxPixels"         , "0");
	Parameter_set(self, "linker.minFill"           , "0.0");
	Parameter_set(self, "linker.maxFill"           , "0.0");
	Parameter_set(self, "linker.positivity"        , "false");
	Parameter_set(self, "linker.keepNegative"      , "false");
	
	// Reliability
	Parameter_set(self, "reliability.enable"       , "false");
	Parameter_set(self, "reliability.parameters"   , "peak, sum, mean");
	Parameter_set(self, "reliability.threshold"    , "0.9");
	Parameter_set(self, "reliability.scaleKernel"  , "0.4");
	Parameter_set(self, "reliability.minSNR"       , "3.0");
	Parameter_set(self, "reliability.minPixels"    , "0");
	Parameter_set(self, "reliability.autoKernel"   , "false");
	Parameter_set(self, "reliability.iterations"   , "30");
	Parameter_set(self, "reliability.tolerance"    , "0.05");
	Parameter_set(self, "reliability.catalog"      , "");
	Parameter_set(self, "reliability.plot"         , "true");
	Parameter_set(self, "reliability.debug"        , "false");
	
	// Mask dilation
	Parameter_set(self, "dilation.enable"          , "false");
	Parameter_set(self, "dilation.iterationsXY"    , "10");
	Parameter_set(self, "dilation.iterationsZ"     , "5");
	Parameter_set(self, "dilation.threshold"       , "0.001");
	
	// Parameterisation
	Parameter_set(self, "parameter.enable"         , "true");
	Parameter_set(self, "parameter.wcs"            , "true");
	Parameter_set(self, "parameter.physical"       , "false");
	Parameter_set(self, "parameter.prefix"         , "SoFiA");
	Parameter_set(self, "parameter.offset"         , "false");
	
	// Output
	Parameter_set(self, "output.directory"         , "");
	Parameter_set(self, "output.filename"          , "");
	Parameter_set(self, "output.writeCatASCII"     , "true");
	Parameter_set(self, "output.writeCatXML"       , "true");
	Parameter_set(self, "output.writeCatSQL"       , "false");
	Parameter_set(self, "output.writeNoise"        , "false");
	Parameter_set(self, "output.writeFiltered"     , "false");
	Parameter_set(self, "output.writeMask"         , "false");
	Parameter_set(self, "output.writeMask2d"       , "false");
	Parameter_set(self, "output.writeRawMask"      , "false");
	Parameter_set(self, "output.writeMoments"      , "false");
	Parameter_set(self, "output.writeCubelets"     , "false");
	Parameter_set(self, "output.marginCubelets"    , "10");
	Parameter_set(self, "output.thresholdMom12"    , "0.0");
	Parameter_set(self, "output.overwrite"         , "true");
	
	return;
}



/// @brief Reallocate memory for one additional parameter setting
///
/// Private method for allocating additional memory for one more
/// parameter in the specified parameter list. Note that this will
/// not create a new parameter setting yet, but will just allocate
/// the memory needed to append a parameter at the end of the
/// current list. The function should be called from public member
/// functions that will add parameters to a parameter list prior to
/// inserting the new parameter name and value at the end.
///
/// @param self  Object self-reference.

PRIVATE void Parameter_append_memory(Parameter *self)
{
	// Extend memory for parameter settings
	self->n_par += 1;
	self->keys   = (String **)memory_realloc(self->keys,   self->n_par, sizeof(String *));
	self->values = (String **)memory_realloc(self->values, self->n_par, sizeof(String *));
	
	// Call constructor for strings
	self->keys  [self->n_par - 1] = String_new("");
	self->values[self->n_par - 1] = String_new("");
	
	return;
}
