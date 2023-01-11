// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (DataCube.c) - Source Finding Application                //
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

/// @file   DataCube.c
/// @author Tobias Westmeier
/// @date   08/12/2021
/// @brief  Class for storage, source finding and parameterisation of FITS data cubes.


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "DataCube.h"
#include "Table.h"
#include "Source.h"
#include "statistics_flt.h"
#include "statistics_dbl.h"



// ----------------------------------------------------------------- //
// Compile-time checks to ensure that                                //
//                                                                   //
//   1. the number of bits per byte is 8 (CHAR_BIT),                 //
//                                                                   //
//   2. sizeof(int) is greater than sizeof(char),                    //
//                                                                   //
//   3. the size of int8_t, int16_t, int32_t and int64_t is exactly  //
//      1, 2, 4 and 8, respectively, and                             //
//                                                                   //
//   4. the size of float and double is 4 and 8, respectively.       //
//                                                                   //
// Without these conditions the code would not function properly. If //
// any of these conditions is not met, a compiler error should be    //
// raised, e.g.: "size of array is negative". This is to ensure that //
// the code can only be run on compliant architectures, as it would  //
// not work correctly on machines that don't fulfil these require-   //
// ments. For that reason these checks must not be disabled to en-   //
// force compilation on non-compliant systems.                       //
// ----------------------------------------------------------------- //

/// \cond doxy_ignore
COMPILE_TIME_CHECK ( CHAR_BIT == 8,              FATAL_Number_of_bits_per_byte_is_not_equal_to_8 );
COMPILE_TIME_CHECK ( sizeof(int) > sizeof(char), FATAL_Size_of_int_is_not_greater_than_size_of_char );
COMPILE_TIME_CHECK ( sizeof(int8_t)  == 1,       FATAL_Size_of_uint8_is_not_equal_to_1 );
COMPILE_TIME_CHECK ( sizeof(int16_t) == 2,       FATAL_Size_of_int16_is_not_equal_to_2 );
COMPILE_TIME_CHECK ( sizeof(int32_t) == 4,       FATAL_Size_of_int32_is_not_equal_to_4 );
COMPILE_TIME_CHECK ( sizeof(int64_t) == 8,       FATAL_Size_of_int64_is_not_equal_to_8 );
COMPILE_TIME_CHECK ( sizeof(float)   == 4,       FATAL_Size_of_float_is_not_equal_to_4 );
COMPILE_TIME_CHECK ( sizeof(double)  == 8,       FATAL_Size_of_double_is_not_equal_to_8 );
/// \endcond



/// @brief Class for storage, source finding and parameterisation of FITS data cubes.
///
/// The purpose of this class is to handle up to three-dimensional
/// astronomical data cubes. The class is intended for reading and
/// manipulating FITS data cubes by providing methods for loading and
/// saving FITS files and manipulating the header and data units of a
/// FITS file. In addition, methods for filtering, source finding and
/// parameterisation are available. Currently, only single-HDU files
/// are supported.

CLASS DataCube
{
	char   *data;          ///< Pointer to `char` array containing data values.
	size_t  data_size;     ///< Total number of data values in cube.
	Header *header;        ///< Pointer to Header object containing header information.
	int     data_type;     ///< FITS data type (-64, -32, 8, 16, 32 or 64).
	int     word_size;     ///< Size of a single datum in multiples of `sizeof(char)`.
	size_t  dimension;     ///< Dimension of the cube (1-4).
	size_t  axis_size[4];  ///< Size of the up-to-four axes of the cube.
	bool    verbosity;     ///< Verbosity level (0 or 1).
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new and empty DataCube
/// object and return a pointer to the newly created object. No
/// memory will be allocated other than for the object itself. Note
/// that the destructor will need to be called explicitly once the
/// object is no longer required to release any memory allocated
/// during the lifetime of the object.
///
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created DataCube object.

PUBLIC DataCube *DataCube_new(const bool verbosity)
{
	DataCube *self = (DataCube*)memory(MALLOC, 1, sizeof(DataCube));
	
	// Initialise properties
	self->data         = NULL;
	self->data_size    = 0;
	self->header       = NULL;
	self->data_type    = 0;
	self->word_size    = 0;
	self->dimension    = 0;
	self->axis_size[0] = 0;
	self->axis_size[1] = 0;
	self->axis_size[2] = 0;
	self->axis_size[3] = 0;
	
	self->verbosity = verbosity;
	
	return self;
}



/// @brief Copy constructor
///
/// Copy constructor. Will create a new DataCube object that is a
/// physical copy of the object pointed to by source. A pointer
/// to the newly created object will be returned. Note that the
/// destructor will need to be called explicitly once the object
/// is no longer required to release any memory allocated to the
/// object.
///
/// @param source  Pointer to DataCube object to be copied.
///
/// @return Pointer to newly created DataCube object.

PUBLIC DataCube *DataCube_copy(const DataCube *source)
{
	// Sanity checks
	check_null(source);
	
	DataCube *self = DataCube_new(source->verbosity);
	
	// Copy header
	self->header = Header_copy(source->header);
	
	// Copy data
	self->data = (char *)memory(MALLOC, source->data_size, source->word_size * sizeof(char));
	memcpy(self->data, source->data, source->word_size * source->data_size);
	
	// Copy remaining properties
	self->data_size    = source->data_size;
	self->data_type    = source->data_type;
	self->word_size    = source->word_size;
	self->dimension    = source->dimension;
	self->axis_size[0] = source->axis_size[0];
	self->axis_size[1] = source->axis_size[1];
	self->axis_size[2] = source->axis_size[2];
	self->axis_size[3] = source->axis_size[3];
	
	return self;
}



/// @brief Variant of standard constructor
///
/// Alternative standard constructor. Will create a new DataCube
/// object with the dimensions and data type specified. Memory for
/// the data array and basic header will be allocated. The array
/// will be initialised with a value of 0. A pointer to the newly
/// created object will be returned. Note that the destructor will
/// need to be called explicitly once the object is no longer
/// required to release any memory allocated to the object.
///
/// @param nx         Size of first axis of data array.
/// @param ny         Size of second axis of data array.
/// @param nz         Size of third axis of data array.
/// @param type       Standard FITS data type (-64, -32, 8, 16,
///                   32 or 64).
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created DataCube object.

PUBLIC DataCube *DataCube_blank(const size_t nx, const size_t ny, const size_t nz, const int type, const bool verbosity)
{
	// Sanity checks
	ensure(nx > 0 && ny > 0 && nz > 0, ERR_USER_INPUT, "Illegal data cube size of (%zu, %zu, %zu) requested.", nx, ny, nz);
	ensure(abs(type) == 64 || abs(type) == 32 || type == 8 || type == 16, ERR_USER_INPUT, "Invalid FITS data type of %d requested.", type);
	
	DataCube *self = DataCube_new(verbosity);
	
	// Set up properties
	self->data_size    = nx * ny * nz;
	self->data_type    = type;
	self->word_size    = abs(type / 8);
	self->dimension    = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
	self->axis_size[0] = nx;
	self->axis_size[1] = ny;
	self->axis_size[2] = nz;
	self->axis_size[3] = 0;
	
	// Get current date and time
	char current_time_string[32];
	time_t current_time = time(NULL);
	strftime(current_time_string, 32, "%Y-%m-%dT%H:%M:%S", gmtime(&current_time));
	
	// Create data array filled with 0
	self->data = (char *)memory(CALLOC, self->data_size, self->word_size * sizeof(char));
	
	// Create empty header (single block with just the END keyword)
	self->header = Header_blank(verbosity);
	
	// Insert required header elements
	Header_set_bool(self->header, "SIMPLE", true);
	Header_set_int(self->header, "BITPIX", self->data_type);
	Header_set_int(self->header, "NAXIS",  self->dimension);
	Header_set_int(self->header, "NAXIS1", self->axis_size[0]);
	if(self->dimension > 1) Header_set_int (self->header, "NAXIS2", self->axis_size[1]);
	if(self->dimension > 2) Header_set_int (self->header, "NAXIS3", self->axis_size[2]);
	
	Header_set_str(self->header, "CTYPE1", " ");
	Header_set_flt(self->header, "CRPIX1", 1.0);
	Header_set_flt(self->header, "CDELT1", 1.0);
	Header_set_flt(self->header, "CRVAL1", 1.0);
	
	if(self->dimension > 1)
	{
		Header_set_str(self->header, "CTYPE2", " ");
		Header_set_flt(self->header, "CRPIX2", 1.0);
		Header_set_flt(self->header, "CDELT2", 1.0);
		Header_set_flt(self->header, "CRVAL2", 1.0);
	}
	
	if(self->dimension > 2)
	{
		Header_set_str(self->header, "CTYPE3", " ");
		Header_set_flt(self->header, "CRPIX3", 1.0);
		Header_set_flt(self->header, "CDELT3", 1.0);
		Header_set_flt(self->header, "CRVAL3", 1.0);
	}
	
	Header_set_str(self->header, "ORIGIN", SOFIA_VERSION_FULL " (" SOFIA_CREATION_DATE ")");
	Header_set_str(self->header, "DATE", current_time_string);
	
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void DataCube_delete(DataCube *self)
{
	if(self != NULL)
	{
		Header_delete(self->header);
		free(self->data);
		free(self);
	}
	
	return;
}



/// @brief Get data array size
///
/// Public method for retrieving the size of the data array of the
/// specified data cube, i.e. the total number of data samples. If
/// a `NULL` pointer is provided, 0 will be returned.
///
/// @param self  Object self-reference.
///
/// @return Size of the data array (number of elements).

PUBLIC size_t DataCube_get_size(const DataCube *self)
{
	return self == NULL ? 0 : self->data_size;
}



/// @brief Get data axis size
///
/// Public method for retrieving the size of the specified axis of
/// the data array. Note that axis must be in the range of 0 to 3.
///
/// @param self  Object self-reference.
/// @param axis  Index of the axis the size of which is needed.
///
/// @return Size of the requested axis in pixels.

PUBLIC size_t DataCube_get_axis_size(const DataCube *self, const size_t axis)
{
	ensure(axis < 4, ERR_USER_INPUT, "Axis must be in the range of 0 to 3.");
	return self == NULL ? 0 : self->axis_size[axis];
}



/// @brief Read data cube from FITS file
///
/// Public method for reading a data cube from a FITS file. The
/// data cube must have between 1 and 3 dimensions. 4-dimensional
/// FITS cubes are also supported as long as the 4th axis is of
/// size 1 (e.g. Stokes I). A region can be specified to read only
/// a portion of the image. The region must be of the form `x_min`,
/// `x_max`, `y_min`, `y_max`, `z_min`, `z_max`. If `NULL`, the
/// full cube will be read in.
///
/// @param self      Object self-reference.
/// @param filename  Name of the input FITS file.
/// @param region    Array of 6 values denoting a region of the cube
///                  to be read in (format: `x_min`, `x_max`, `y_min`,
///                  `y_max`, `z_min`, `z_max`). Set to `NULL` to read
///                  entire data cube.

PUBLIC void DataCube_load(DataCube *self, const char *filename, const Array_siz *region)
{
	// Sanity checks
	check_null(self);
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "Empty file name provided.");
	
	// Check region specification
	if(region != NULL)
	{
		ensure(Array_siz_get_size(region) == 6, ERR_USER_INPUT, "Invalid region supplied; must contain 6 values.");
		for(size_t i = 0; i < Array_siz_get_size(region); i += 2) ensure(Array_siz_get(region, i) <= Array_siz_get(region, i + 1), ERR_USER_INPUT, "Invalid region supplied; minimum greater than maximum.");
	}
	
	// Open FITS file
	message("Opening FITS file \'%s\'.", filename);
	FILE *fp = fopen(filename, "rb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open FITS file \'%s\'.", filename);
		
	// Read entire header into temporary array
	char *header = NULL;
	size_t header_size = 0;
	bool end_reached = false;
	
	while(!end_reached)
	{
		// (Re-)allocate memory as needed
		header = (char *)memory_realloc(header, header_size + FITS_HEADER_BLOCK_SIZE, sizeof(char));
		
		// Read header block
		ensure(fread(header + header_size, 1, FITS_HEADER_BLOCK_SIZE, fp) == FITS_HEADER_BLOCK_SIZE, ERR_FILE_ACCESS, "FITS file ended unexpectedly while reading header.");
		
		// Check if we have reached the end of the header
		char *ptr = header + header_size;
		
		while(!end_reached && ptr < header + header_size + FITS_HEADER_BLOCK_SIZE)
		{
			if(strncmp(ptr, "END", 3) == 0) end_reached = true;
			else ptr += FITS_HEADER_LINE_SIZE;
		}
		
		// Set header size parameter
		header_size += FITS_HEADER_BLOCK_SIZE;
	}
	
	// Check if valid FITS file
	ensure(strncmp(header, "SIMPLE", 6) == 0, ERR_USER_INPUT, "Missing \'SIMPLE\' keyword; file does not appear to be a FITS file.");
	
	// Create Header object and de-allocate memory again
	self->header = Header_new(header, header_size, self->verbosity);
	free(header);
	
	// Extract crucial header elements
	self->data_type    = Header_get_int(self->header, "BITPIX");
	self->dimension    = Header_get_int(self->header, "NAXIS");
	self->axis_size[0] = Header_get_int(self->header, "NAXIS1");
	self->axis_size[1] = Header_get_int(self->header, "NAXIS2");
	self->axis_size[2] = Header_get_int(self->header, "NAXIS3");
	self->axis_size[3] = Header_get_int(self->header, "NAXIS4");
	self->word_size    = abs(self->data_type) / 8;             // WARNING: Assumes 8 bits per char; see CHAR_BIT in limits.h.
	self->data_size    = self->axis_size[0];
	for(size_t i = 1; i < self->dimension; ++i) self->data_size *= self->axis_size[i];
	
	// Sanity checks
	ensure(self->data_type == -64
		|| self->data_type == -32
		|| self->data_type == 8
		|| self->data_type == 16
		|| self->data_type == 32
		|| self->data_type == 64,
		ERR_USER_INPUT, "Invalid BITPIX keyword encountered.");
	
	ensure(self->dimension > 0
		&& self->dimension < 5,
		ERR_USER_INPUT, "Only FITS files with 1-4 dimensions are supported.");
	
	ensure(self->dimension < 4
		|| self->axis_size[3] == 1
		|| self->axis_size[2] == 1,
		ERR_USER_INPUT, "The size of the 3rd or 4th axis must be 1.");
	
	ensure(self->data_size > 0,
		ERR_USER_INPUT, "Invalid NAXISn keyword encountered.");
	
	if(self->dimension < 3) self->axis_size[2] = 1;
	if(self->dimension < 2) self->axis_size[1] = 1;
	
	// Swap third and fourth axis header keywords if necessary
	if(self->dimension == 4 && self->axis_size[2] == 1 && self->axis_size[3] > 1)
	{
		warning("Swapping order of 3rd and 4th axis of 4D cube.");
		
		double tmp;
		size_t tmp2;
		String *str3, *str4;
		
		tmp2 = self->axis_size[2];
		self->axis_size[2] = self->axis_size[3];
		self->axis_size[3] = tmp2;
		
		Header_set_int(self->header, "NAXIS3", Header_get_int(self->header, "NAXIS4"));
		Header_set_int(self->header, "NAXIS4", 1);
		
		tmp = Header_get_flt(self->header, "CRPIX3");
		Header_set_flt(self->header, "CRPIX3", Header_get_flt(self->header, "CRPIX4"));
		Header_set_flt(self->header, "CRPIX4", tmp);
		
		tmp = Header_get_flt(self->header, "CRVAL3");
		Header_set_flt(self->header, "CRVAL3", Header_get_flt(self->header, "CRVAL4"));
		Header_set_flt(self->header, "CRVAL4", tmp);
		
		tmp = Header_get_flt(self->header, "CDELT3");
		Header_set_flt(self->header, "CDELT3", Header_get_flt(self->header, "CDELT4"));
		Header_set_flt(self->header, "CDELT4", tmp);
		
		str3 = Header_get_string(self->header, "CTYPE3");
		str4 = Header_get_string(self->header, "CTYPE4");
		Header_set_str(self->header, "CTYPE3", String_get(str4));
		Header_set_str(self->header, "CTYPE4", String_get(str3));
		String_delete(str3);
		String_delete(str4);
		
		str3 = Header_get_string(self->header, "CUNIT3");
		str4 = Header_get_string(self->header, "CUNIT4");
		Header_set_str(self->header, "CUNIT3", String_get(str4));
		Header_set_str(self->header, "CUNIT4", String_get(str3));
		String_delete(str3);
		String_delete(str4);
	}
	
	// Work out region
	const size_t x_min = (region != NULL && Array_siz_get(region, 0) > 0) ? Array_siz_get(region, 0) : 0;
	const size_t x_max = (region != NULL && Array_siz_get(region, 1) < self->axis_size[0] - 1) ? Array_siz_get(region, 1) : self->axis_size[0] - 1;
	const size_t y_min = (region != NULL && Array_siz_get(region, 2) > 0) ? Array_siz_get(region, 2) : 0;
	const size_t y_max = (region != NULL && Array_siz_get(region, 3) < self->axis_size[1] - 1) ? Array_siz_get(region, 3) : self->axis_size[1] - 1;
	const size_t z_min = (region != NULL && Array_siz_get(region, 4) > 0) ? Array_siz_get(region, 4) : 0;
	const size_t z_max = (region != NULL && Array_siz_get(region, 5) < self->axis_size[2] - 1) ? Array_siz_get(region, 5) : self->axis_size[2] - 1;
	
	ensure(x_min <= x_max && y_min <= y_max && z_min <= z_max, ERR_USER_INPUT, "Invalid data cube region requested.");
	
	const size_t region_nx = x_max - x_min + 1;
	const size_t region_ny = y_max - y_min + 1;
	const size_t region_nz = z_max - z_min + 1;
	const size_t region_size = region_nx * region_ny * region_nz;
	
	// Print status information
	message("Reading FITS data with the following specifications:");
	message("  Data type:    %d", self->data_type);
	message("  No. of axes:  %zu", self->dimension);
	message("  Axis sizes:   %zu, %zu, %zu", self->axis_size[0], self->axis_size[1], self->axis_size[2]);
	message("  Region:       %zu-%zu, %zu-%zu, %zu-%zu", x_min, x_max, y_min, y_max, z_min, z_max);
	message("  Memory used:  %.1f MB", (double)(region_size * self->word_size) / MEGABYTE);
	
	// Allocate memory for data array
	self->data = (char *)memory_realloc(self->data, region_size, self->word_size * sizeof(char));
	
	// Read data
	if(region == NULL)
	{
		// No region supplied -> read full cube
		ensure(fread(self->data, self->word_size, self->data_size, fp) == self->data_size, ERR_FILE_ACCESS, "FITS file ended unexpectedly while reading data.");
	}
	else
	{
		// Region supplied -> read sub-cube
		char *ptr_data = self->data;
		const size_t fp_start = (size_t)ftell(fp); // Start position of data array in file
		
		// Create buffer for a single set of image rows
		const size_t buffer_size = self->axis_size[0] * region_ny;
		char *buffer = (char *)memory(MALLOC, buffer_size, self->word_size * sizeof(char));
		
		const size_t bytes_per_data_row = self->axis_size[0] * self->word_size;
		const size_t bytes_per_region_row = region_nx * self->word_size;
		
		// Read relevant data segments
		for(size_t z = z_min; z <= z_max; ++z)
		{
			progress_bar("Progress: ", z - z_min, z_max - z_min);
			
			// Point file pointer to start of image segment to be read
			ensure(!fseek(fp, fp_start + DataCube_get_index(self, 0, y_min, z) * self->word_size, SEEK_SET), ERR_FILE_ACCESS, "Error while reading FITS file.");
			
			// Read image plane into buffer
			ensure(fread(buffer, self->word_size, buffer_size, fp) == buffer_size, ERR_FILE_ACCESS, "FITS file ended unexpectedly while reading data.");
			
			// Point to start of first row to be extracted
			char *ptr_buffer = buffer + x_min * self->word_size;
			
			// Copy relevant segments into data array
			for(size_t y = y_min; y <= y_max; ++y)
			{
				memcpy(ptr_data, ptr_buffer, bytes_per_region_row);
				ptr_buffer += bytes_per_data_row;
				ptr_data   += bytes_per_region_row;
			}
		}
		
		// Delete buffer again
		free(buffer);
		
		// Update object properties
		// NOTE: This must happen after reading the sub-cube, as the full
		///     cube dimensions must be known during data extraction.
		self->data_size = region_size;
		self->axis_size[0] = region_nx;
		self->axis_size[1] = region_ny;
		self->axis_size[2] = region_nz;
		
		// Adjust WCS information in header
		Header_adjust_wcs_to_subregion(self->header, x_min, x_max, y_min, y_max, z_min, z_max);
	}
	
	// Close FITS file
	fclose(fp);
	
	// Swap byte order if required
	DataCube_swap_byte_order(self);
	
	// Handle BSCALE and BZERO if necessary
	const double bscale = Header_get_flt(self->header, "BSCALE");
	const double bzero  = Header_get_flt(self->header, "BZERO");
	
	if((IS_NOT_NAN(bscale) && bscale != 1.0) || (IS_NOT_NAN(bzero) && bzero != 0.0))
	{
		// Scaling required
		if(self->data_type < 0.0)
		{
			// Floating-point data; simply print warning...
			warning("Applying non-trivial BSCALE and BZERO to floating-point data.");
			
			// ...and scale data
			if(bscale != 1.0) DataCube_multiply_const(self, bscale);
			if(bzero  != 0.0) DataCube_add_const(self, bzero);
			
			// Update header
			Header_remove(self->header, "BSCALE");
			Header_remove(self->header, "BZERO");
		}
		else
		{
			// Integer data; conversion to 32-bit floating-point data required
			warning("Applying non-trivial BSCALE and BZERO to integer data\n         and converting to 32-bit floating-point type.");
			
			// Create 32-bit array
			float *data_copy = (float *)memory(MALLOC, self->axis_size[0] * self->axis_size[1] * self->axis_size[2], sizeof(float));
			float *ptr = data_copy;
			
			// Check for blanking value
			const bool blanking_required = (Header_check(self->header, "BLANK") > 0);
			const long int blanking_value = blanking_required ? Header_get_int(self->header, "BLANK") : 0;
			long int value;
			
			// Copy scaled data over
			for(size_t z = 0; z < self->axis_size[2]; ++z)
			{
				for(size_t y = 0; y < self->axis_size[1]; ++y)
				{
					for(size_t x = 0; x < self->axis_size[0]; ++x)
					{
						value = DataCube_get_data_int(self, x, y, z);
						if(blanking_required && blanking_value == value) *ptr = NAN;
						else *ptr = bzero + bscale * value;
						++ptr;
					}
				}
			}
			
			// Delete original array and point to new copy instead
			free(self->data);
			self->data = (char *)data_copy;
			
			// Update header
			Header_set_int(self->header, "BITPIX", -32);
			Header_remove(self->header, "BSCALE");
			Header_remove(self->header, "BZERO");
			Header_remove(self->header, "BLANK");
			
			// Update object properties
			self->data_type = -32;
			self->word_size = 4;
		}
	}
	
	return;
}



/// @brief Write data cube into FITS file
///
/// Public method for writing the current data cube object
/// (referenced by `self`) into a FITS file. The function will
/// terminate the current programme execution if an error is
/// encountered during the write process. If the output file
/// already exists, it will be overwritten only if overwrite
/// is set to true.
///
/// @param self        Object self-reference.
/// @param filename    Name of output FITS file.
/// @param overwrite   If `true`, overwrite existing file. Otherwise
///                    terminate if the file already exists.
/// @param preserve    If `true`, ensure that the data array is in the
///                    correct byte order after returning. If `false`,
///                    the byte order may be corrupted, which may be
///                    faster and acceptable if the data are no
///                    longer needed afterwards.

PUBLIC void DataCube_save(const DataCube *self, const char *filename, const bool overwrite, const bool preserve)
{
	// Sanity checks
	check_null(self);
	check_null(filename);
	ensure(strlen(filename), ERR_USER_INPUT, "Empty file name provided.");
	
	// Open FITS file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to create new FITS file: %s\n       Does the destination exist and is writeable?", filename);
	
	message("Creating FITS file: %s", strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/') + 1);
	
	// Write entire header
	ensure(fwrite(Header_get(self->header), 1, Header_get_size(self->header), fp) == Header_get_size(self->header), ERR_FILE_ACCESS, "Failed to write header to FITS file.");
	
	// Swap byte order of array in memory if necessary
	DataCube_swap_byte_order(self);
	
	// Write entire data array
	ensure(fwrite(self->data, self->word_size, self->data_size, fp) == self->data_size, ERR_FILE_ACCESS, "Failed to write data to FITS file.");
	
	// Fill file with 0x00 if necessary
	const size_t size_footer = ((self->data_size * self->word_size) % FITS_HEADER_BLOCK_SIZE);
	if(size_footer)
	{
		const char footer = '\0';
		for(size_t counter = FITS_HEADER_BLOCK_SIZE - size_footer; counter--;) fwrite(&footer, 1, 1, fp);
	}
	
	// Close file
	fclose(fp);
	
	// Revert to original byte order if necessary and requested
	if(preserve) DataCube_swap_byte_order(self);
	
	return;
}



/// @brief Wrapper around commonly needed Header method
///
/// See Header_get_int() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC long int DataCube_gethd_int(const DataCube *self, const char *key) {
	return Header_get_int(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_get_flt() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC double DataCube_gethd_flt(const DataCube *self, const char *key) {
	return Header_get_flt(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_get_bool() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC bool DataCube_gethd_bool(const DataCube *self, const char *key) {
	return Header_get_bool(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_get_string() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC String *DataCube_gethd_string(const DataCube *self, const char *key) {
	return Header_get_string(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_get_str() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_gethd_str(const DataCube *self, const char *key, char *value) {
	return Header_get_str(self->header, key, value);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_set_int() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_puthd_int(DataCube *self, const char *key, const long int value) {
	return Header_set_int(self->header, key, value);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_set_flt() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_puthd_flt(DataCube *self, const char *key, const double value) {
	return Header_set_flt(self->header, key, value);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_set_bool() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_puthd_bool(DataCube *self, const char *key, const bool value) {
	return Header_set_bool(self->header, key, value);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_set_str() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_puthd_str(DataCube *self, const char *key, const char *value) {
	return Header_set_str(self->header, key, value);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_check() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC size_t DataCube_chkhd(const DataCube *self, const char *key) {
	return Header_check(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_compare() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC bool DataCube_cmphd(const DataCube *self, const char *key, const char *value, const size_t n) {
	return Header_compare(self->header, key, value, n);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_remove() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC int DataCube_delhd(DataCube *self, const char *key) {
	return Header_remove(self->header, key);
}

/// @brief Wrapper around commonly needed Header method
///
/// See Header_copy_wcs() in class Header for detailed information on
/// the respective method and its arguments and return value.

PUBLIC void DataCube_copy_wcs(const DataCube *source, DataCube *target) {
	Header_copy_wcs(source->header, target->header);
	return;
}



/// @brief Add SoFiA settings to FITS file history
///
/// Public method for writing the SoFiA parameter settings stored
/// in 'par' as FITS header HISTORY items. If `par` is `NULL`, then
/// no history will be written.
///
/// @param self  Object self-reference.
/// @param par   SoFiA parameter settings object.

PUBLIC void DataCube_add_history(DataCube *self, const Parameter *par)
{
	// Sanity checks
	check_null(self);
	if(par == NULL) return;
	
	String *value = String_new("");
	Header_comment(self->header, "SoFiA control parameters", true);
	
	for(size_t i = 0; i < Parameter_get_size(par); ++i)
	{
		// Construct value string
		String_set(value, Parameter_get_key(par, i));
		String_append(value, " = ");
		String_append(value, Parameter_get_str_index(par, i));
		
		// Sanitise value string
		for(size_t j = 0; j < String_size(value); ++j)
		{
			const char c = String_at(value, j);
			if(c < '\x20' || c > '\x7E') String_set_char(value, j, '\x23');  // Replace with '#'
		}
		
		// Write history entry
		Header_comment(self->header, String_get(value), true);
	}
	
	String_delete(value);
	return;
}



/// @brief Read data value as double-precision floating-point number
///
/// Public method to extract the data value at the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. The function will return
/// the result as a double-precision floating-point value irrespective
/// of the native data type of the FITS file.
///
/// @param self  Object self-reference.
/// @param x     First coordinate.
/// @param y     Second coordinate.
/// @param z     Third coordinate.
///
/// @return Returns the value of the data array at the given position
///         as a double-precision floating-point value.

PUBLIC double DataCube_get_data_flt(const DataCube *self, const size_t x, const size_t y, const size_t z)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position (%zu, %zu, %zu) outside of image boundaries.", x, y, z);
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type)
	{
		case -64:
			return *((double *)(self->data + i * self->word_size));
		case -32:
			return (double)(*((float *)(self->data + i * self->word_size)));
		case 8:
			return (double)(*((uint8_t *)(self->data + i * self->word_size)));
		case 16:
			return (double)(*((int16_t *)(self->data + i * self->word_size)));
		case 32:
			return (double)(*((int32_t *)(self->data + i * self->word_size)));
		case 64:
			return (double)(*((int64_t *)(self->data + i * self->word_size)));
	}
	
	return NAN;
}



/// @brief Read data value as long integer number
///
/// Public method to extract the data value at the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. The function will return
/// the result as a long integer value irrespective of the native data
/// type of the FITS file.
///
/// @param self  Object self-reference.
/// @param x     First coordinate.
/// @param y     Second coordinate.
/// @param z     Third coordinate.
///
/// @return Returns the value of the data array at the given position as a
///         long integer value.

PUBLIC long int DataCube_get_data_int(const DataCube *self, const size_t x, const size_t y, const size_t z)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position (%zu, %zu, %zu) outside of image boundaries.", x, y, z);
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type)
	{
		case -64:
			return (long int)(*((double *)(self->data + i * self->word_size)));
		case -32:
			return (long int)(*((float *)(self->data + i * self->word_size)));
		case 8:
			return (long int)(*((uint8_t *)(self->data + i * self->word_size)));
		case 16:
			return (long int)(*((int16_t *)(self->data + i * self->word_size)));
		case 32:
			return (long int)(*((int32_t *)(self->data + i * self->word_size)));
		case 64:
			return (long int)(*((int64_t *)(self->data + i * self->word_size)));
	}
	
	return 0;
}



/// @brief Set data value as double-precision floating-point number
///
/// Public method to write the data value to the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. Note that the data value
/// will be cast to the native data type of the array before being written.
///
/// @param self   Object self-reference.
/// @param x      First coordinate.
/// @param y      Second coordinate.
/// @param z      Third coordinate.
/// @param value  Data value to be written to array.

PUBLIC void DataCube_set_data_flt(DataCube *self, const size_t x, const size_t y, const size_t z, const double value)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position outside of image boundaries.");
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type)
	{
		case -64:
			*((double *)(self->data + i * self->word_size))  = value;
			break;
		case -32:
			*((float *)(self->data + i * self->word_size))   = (float)value;
			break;
		case 8:
			*((uint8_t *)(self->data + i * self->word_size)) = (uint8_t)value;
			break;
		case 16:
			*((int16_t *)(self->data + i * self->word_size)) = (int16_t)value;
			break;
		case 32:
			*((int32_t *)(self->data + i * self->word_size)) = (int32_t)value;
			break;
		case 64:
			*((int64_t *)(self->data + i * self->word_size)) = (int64_t)value;
			break;
	}
	
	return;
}



/// @brief Add data value as double-precision floating-point number
///
/// Public method to add the data value to the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. Note that the data value
/// will be cast to the native data type of the array before being written.
///
/// @param self   Object self-reference.
/// @param x      First coordinate.
/// @param y      Second coordinate.
/// @param z      Third coordinate.
/// @param value  Data value to be added to array.

PUBLIC void DataCube_add_data_flt(DataCube *self, const size_t x, const size_t y, const size_t z, const double value)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position outside of image boundaries.");
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type)
	{
		case -64:
			*((double *)(self->data + i * self->word_size))  += value;
			break;
		case -32:
			*((float *)(self->data + i * self->word_size))   += (float)value;
			break;
		case 8:
			*((uint8_t *)(self->data + i * self->word_size)) += (uint8_t)value;
			break;
		case 16:
			*((int16_t *)(self->data + i * self->word_size)) += (int16_t)value;
			break;
		case 32:
			*((int32_t *)(self->data + i * self->word_size)) += (int32_t)value;
			break;
		case 64:
			*((int64_t *)(self->data + i * self->word_size)) += (int64_t)value;
			break;
	}
	
	return;
}



/// @brief Set data value as long integer number
///
/// Public method to write the data value to the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. Note that the data value
/// will be cast to the native data type of the array before being written.
///
/// @param self   Object self-reference.
/// @param x      First coordinate.
/// @param y      Second coordinate.
/// @param z      Third coordinate.
/// @param value  Data value to be written to array.

PUBLIC void DataCube_set_data_int(DataCube *self, const size_t x, const size_t y, const size_t z, const long int value)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position outside of image boundaries.");
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type) {
		case -64:
			*((double *)(self->data + i * self->word_size))  = (double)value;
			break;
		case -32:
			*((float *)(self->data + i * self->word_size))   = (float)value;
			break;
		case 8:
			*((uint8_t *)(self->data + i * self->word_size)) = (uint8_t)value;
			break;
		case 16:
			*((int16_t *)(self->data + i * self->word_size)) = (int16_t)value;
			break;
		case 32:
			*((int32_t *)(self->data + i * self->word_size)) = (int32_t)value;
			break;
		case 64:
			*((int64_t *)(self->data + i * self->word_size)) = (int64_t)value;
			break;
	}
	
	return;
}



/// @brief Add data value as long integer number
///
/// Public method to add the data value to the specified position
/// (`x`, `y`, `z`), where `x` indexes the first axis, `y` the second
/// axis and `z` the third axis of the cube. Note that the data value
/// will be cast to the native data type of the array before being written.
///
/// @param self   Object self-reference.
/// @param x      First coordinate.
/// @param y      Second coordinate.
/// @param z      Third coordinate.
/// @param value  Data value to be added to array.

PUBLIC void DataCube_add_data_int(DataCube *self, const size_t x, const size_t y, const size_t z, const long int value)
{
	//check_null(self);
	//check_null(self->data);
	//ensure(x < self->axis_size[0] && y < self->axis_size[1] && z < self->axis_size[2], ERR_INDEX_RANGE, "Position outside of image boundaries.");
	const size_t i = DataCube_get_index(self, x, y, z);
	
	switch(self->data_type) {
		case -64:
			*((double *)(self->data + i * self->word_size))  += (double)value;
			break;
		case -32:
			*((float *)(self->data + i * self->word_size))   += (float)value;
			break;
		case 8:
			*((uint8_t *)(self->data + i * self->word_size)) += (uint8_t)value;
			break;
		case 16:
			*((int16_t *)(self->data + i * self->word_size)) += (int16_t)value;
			break;
		case 32:
			*((int32_t *)(self->data + i * self->word_size)) += (int32_t)value;
			break;
		case 64:
			*((int64_t *)(self->data + i * self->word_size)) += (int64_t)value;
			break;
	}
	
	return;
}



/// @brief Fill data cube with floating-point value
///
/// Public method for filling the data cube with the specified
/// floating-point value. Note that this will only be possible if
/// the data cube is of 32 or 64-bit floating-point type.
///
/// @param self   Object self-reference.
/// @param value  Data value to be written to array.

PUBLIC void DataCube_fill_flt(DataCube *self, const double value)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot fill integer array with floating-point value.");
	
	if(self->data_type == -32)
	{
		float *ptr = (float *)(self->data) + self->data_size;
		while(ptr --> (float *)(self->data)) *ptr = value;
	}
	else
	{
		double *ptr = (double *)(self->data) + self->data_size;
		while(ptr --> (double *)(self->data)) *ptr = value;
	}
	
	return;
}



/// @brief Divide a data cube by another cube
///
/// Public method for dividing a data cube by another one. Both
/// cubes must be of floating-point type and need to have the same
/// size. The dividend will be set to `NaN` in places where the divisor
/// is zero.
///
/// @param self     Object self-reference.
/// @param divisor  Data cube to divide by.

PUBLIC void DataCube_divide(DataCube *self, const DataCube *divisor)
{
	// Sanity checks
	check_null(self);
	check_null(divisor);
	check_null(self->data);
	check_null(divisor->data);
	ensure((self->data_type == -32 || self->data_type == -64) && (divisor->data_type == -32 || divisor->data_type == -64), ERR_USER_INPUT, "Dividend and divisor cubes must be of floating-point type.");
	ensure(self->axis_size[0] == divisor->axis_size[0] && self->axis_size[1] == divisor->axis_size[1] && self->axis_size[2] == divisor->axis_size[2], ERR_USER_INPUT, "Dividend and divisor cubes have different sizes.");
	
	if(self->data_type == -32)
	{
		if(divisor->data_type == -32)
		{
			float *ptr_data    = (float *)(self->data);
			float *ptr_divisor = (float *)(divisor->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				if(*(ptr_divisor + i) != 0.0) *(ptr_data + i) /= *(ptr_divisor + i);
				else *(ptr_data + i) = NAN;
			}
		}
		else
		{
			float  *ptr_data    = (float *)(self->data);
			double *ptr_divisor = (double *)(divisor->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				if(*(ptr_divisor + i) != 0.0) *(ptr_data + i) /= *(ptr_divisor + i);
				else *(ptr_data + i) = NAN;
			}
		}
	}
	else
	{
		if(divisor->data_type == -32)
		{
			double *ptr_data    = (double *)(self->data);
			float  *ptr_divisor = (float *)(divisor->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				if(*(ptr_divisor + i) != 0.0) *(ptr_data + i) /= *(ptr_divisor + i);
				else *(ptr_data + i) = NAN;
			}
		}
		else
		{
			double *ptr_data    = (double *)(self->data);
			double *ptr_divisor = (double *)(divisor->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				if(*(ptr_divisor + i) != 0.0) *(ptr_data + i) /= *(ptr_divisor + i);
				else *(ptr_data + i) = NAN;
			}
		}
	}
	
	return;
}



/// @brief Multiply by square root of weights cube
///
/// Public method for multiplying a data cube by the square root of
/// a weights cube. Both cubes must be of floating-point type and
/// need to have the same size. The cube will be set to `NaN` in places
/// where the weights cube is `NaN`.
///
/// @param self     Object self-reference.
/// @param weights  Weights cube to be applied.

PUBLIC void DataCube_apply_weights(DataCube *self, const DataCube *weights)
{
	// Sanity checks
	check_null(self);
	check_null(weights);
	check_null(self->data);
	check_null(weights->data);
	ensure((self->data_type == -32 || self->data_type == -64) && (weights->data_type == -32 || weights->data_type == -64), ERR_USER_INPUT, "Data and weights cubes must be of floating-point type.");
	ensure(self->axis_size[0] == weights->axis_size[0] && self->axis_size[1] == weights->axis_size[1] && self->axis_size[2] == weights->axis_size[2], ERR_USER_INPUT, "Data and weights cubes have different sizes.");
	
	if(self->data_type == -32)
	{
		if(weights->data_type == -32)
		{
			float *ptr_data    = (float *)(self->data);
			float *ptr_weights = (float *)(weights->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i) *(ptr_data + i) *= sqrt(*(ptr_weights + i));
		}
		else
		{
			float  *ptr_data    = (float *)(self->data);
			double *ptr_weights = (double *)(weights->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i) *(ptr_data + i) *= sqrt(*(ptr_weights + i));
		}
	}
	else
	{
		if(weights->data_type == -32)
		{
			double *ptr_data    = (double *)(self->data);
			float  *ptr_weights = (float *)(weights->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i) *(ptr_data + i) *= sqrt(*(ptr_weights + i));
		}
		else
		{
			double *ptr_data    = (double *)(self->data);
			double *ptr_weights = (double *)(weights->data);
			
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i) *(ptr_data + i) *= sqrt(*(ptr_weights + i));
		}
	}
	
	return;
}



/// @brief Multiply data cube by constant factor
///
/// Public method for multiplying a data cube by a constant factor.
/// The data cube must be of floating-point type.
///
/// @param self    Object self-reference.
/// @param factor  Factor to multiply by.

PUBLIC void DataCube_multiply_const(DataCube *self, const double factor)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cube must be of floating-point type for multiplication.");
	
	if(self->data_type == -32) {
		for(float *ptr = (float *)(self->data) + self->data_size; ptr --> (float *)(self->data);) *ptr *= factor;
	}
	else {
		for(double *ptr = (double *)(self->data) + self->data_size; ptr --> (double *)(self->data);) *ptr *= factor;
	}
	
	return;
}



/// @brief Add constant to data cube
///
/// Public method for adding a constant to the entire data cube.
/// The data cube must be of floating-point type.
///
/// @param self     Object self-reference.
/// @param summand  Constant to be added to data cube.

PUBLIC void DataCube_add_const(DataCube *self, const double summand)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cube must be of floating-point type for addition.");
	
	if(self->data_type == -32) {
		for(float *ptr = (float *)(self->data) + self->data_size; ptr --> (float *)(self->data);) *ptr += summand;
	}
	else {
		for(double *ptr = (double *)(self->data) + self->data_size; ptr --> (double *)(self->data);) *ptr += summand;
	}
	
	return;
}



/// @brief Calculate the standard deviation about a value
///
/// Public method for calculating the standard deviation of the
/// data array about a specified value. The cadence specifies which
/// fraction of the elements in the array will be used in the
/// calculation; it can be set to > 1 for large arrays in order to
/// reduce the processing time of the algorithm. The range defines
/// the flux range to be used; if set to 0, all pixels will be used,
/// whereas a positive or negative value indicates that only positive
/// or negative pixels should be used, respectively. This is useful
/// for increasing the robustness of the standard deviation in the
/// presence of negative or positive flux or artefacts in the data.
///
/// @param self     Object self-reference.
/// @param value    Value about which to calculate the standard
///                 deviation.
/// @param cadence  Cadence used in the calculation, i.e. a cadence
///                 of N will calculate the standard deviation using
///                 every N-th element from the array.
/// @param range    Flux range to be used in the calculation. Options
///                 are 0 (entire flux range), -1 (negative fluxes
///                 only) and +1 (positive fluxes only).
///
/// @return Standard deviation about the specified value.

PUBLIC double DataCube_stat_std(const DataCube *self, const double value, const size_t cadence, const int range)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot evaluate standard deviation for integer array.");
	
	if(self->data_type == -32) return std_dev_val_flt((float *)self->data, self->data_size, value, cadence ? cadence : 1, range);
	else return std_dev_val_dbl((double *)self->data, self->data_size, value, cadence ? cadence : 1, range);
}



/// @brief Calculate the median absolute deviation of the array
///
/// Public method for calculating the median absolute deviation re
/// relative to the specified value.
///
/// @param self     Object self-reference.
/// @param value    Value relative to which to calculate the MAD.
/// @param cadence  Cadence used in the calculation, i.e. a cadence
///                 of N will calculate the standard deviation using
///                 every N-th element from the array.
/// @param range    Flux range to be used in the calculation. Options
///                 are 0 (entire flux range), -1 (negative fluxes
///                 only) and +1 (positive fluxes only).
///
/// @return Median absolute deviation of the data array.
///         were found.
///
/// @note A copy of (parts of) the data array will need to be made
///       in order to calculate the median of the data as part of
///       this process.

PUBLIC double DataCube_stat_mad(const DataCube *self, const double value, const size_t cadence, const int range)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot evaluate MAD for integer array.");
	
	// Derive MAD of data copy
	if(self->data_type == -32) return mad_val_flt((float *)self->data, self->data_size, value, cadence ? cadence : 1, range);
	return mad_val_dbl((double *)self->data, self->data_size, value, cadence ? cadence : 1, range);
}



/// @brief Calculate the noise via Gaussian fitting to flux histogram
///
/// Public method for determining the noise level in the data array
/// by fitting a Gaussian function to the flux histogram and determining
/// the standard deviation of that Gaussian. The cadence specifies which
/// fraction of the elements in the array will be used in the
/// calculation; it can be set to > 1 for large arrays in order to
/// reduce the processing time of the algorithm. The range defines the
/// flux range to be used; if set to 0, all pixels will be used, whereas
/// a positive or negative value indicates that only positive or
/// negative pixels should be used, respectively. This is useful for
/// increasing the robustness of the standard deviation in the presence
/// of negative or positive flux or artefacts in the data.
///
/// @param self     Object self-reference.
/// @param cadence  Cadence used in the calculation, i.e. a cadence
///                 of N will calculate the flux histogram using
///                 every N-th element from the array.
/// @param range    Flux range to be used in the calculation. Options
///                 are 0 (entire flux range), -1 (negative fluxes
///                 only) and +1 (positive fluxes only).
///
/// @return Standard deviation of the Gaussian fitted to the histogram.

PUBLIC double DataCube_stat_gauss(const DataCube *self, const size_t cadence, const int range)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot evaluate standard deviation for integer array.");
	
	if(self->data_type == -32) return gaufit_flt((float *)self->data, self->data_size, cadence ? cadence : 1, range);
	else return gaufit_dbl((double *)self->data, self->data_size, cadence ? cadence : 1, range);
}



/// @brief Global noise scaling along spectral axis
///
/// Public method for dividing the data cube by the global noise
/// level as a function of frequency as measured in each spatial
/// plane of the cube. The statistic and flux range used in the
/// noise measurement can be selected to ensure a robust noise
/// measurement. This method should be applied prior to source
/// finding on data cubes where the noise level varies with
/// frequency, but is constant along the two spatial axes in each
/// channel.
///
/// @param self       Object self-reference.
/// @param statistic  Statistic to use in noise measurement. Can be
///                   `NOISE_STAT_STD` for standard deviation,
///                   `NOISE_STAT_MAD` for median absolute deviation or
///                   `NOISE_STAT_GAUSS` for Gaussian fitting to the
///                   flux histogram.
/// @param range      Flux range to be used in noise measurement. Can
///                   be -1, 0 or +1 for negative range, full range
///                   or positive range, respectively.
///
/// @return Pointer to Array object containing the measured noise
///         values for each channel.
///
/// @note The caller will be responsible for calling the destructor
///       on the returned Array object once it is no longer needed.

PUBLIC Array_dbl *DataCube_scale_noise_spec(const DataCube *self, const noise_stat statistic, const int range)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot run noise scaling on integer array.");
	
	// A few settings
	const size_t size_xy = self->axis_size[0] * self->axis_size[1];
	const size_t size_z  = self->axis_size[2];
	double rms;
	Array_dbl *noise_spectrum = Array_dbl_new(size_z);
	
	message("Dividing by noise in each image plane.");
	size_t progress = 0;
	const size_t progress_max = size_z - 1;
	
	#pragma omp parallel for schedule(static) private(rms)
	for(size_t i = 0; i < size_z; ++i)
	{
		#pragma omp critical
		progress_bar("Progress: ", progress++, progress_max);
		
		if(self->data_type == -32)
		{
			float *ptr_start = (float *)(self->data) + i * size_xy;
			
			if(statistic == NOISE_STAT_STD) rms = std_dev_val_flt(ptr_start, size_xy, 0.0, 1, range);
			else if(statistic == NOISE_STAT_MAD) rms = MAD_TO_STD * mad_val_flt(ptr_start, size_xy, 0.0, 1, range);
			else rms = gaufit_flt(ptr_start, size_xy, 1, range);
			
			for(float *ptr = ptr_start + size_xy; ptr --> ptr_start;) *ptr /= rms;
			Array_dbl_set(noise_spectrum, i, rms);
		}
		else
		{
			double *ptr_start = (double *)(self->data) + i * size_xy;
			
			if(statistic == NOISE_STAT_STD) rms = std_dev_val_dbl(ptr_start, size_xy, 0.0, 1, range);
			else if(statistic == NOISE_STAT_MAD) rms = MAD_TO_STD * mad_val_dbl(ptr_start, size_xy, 0.0, 1, range);
			else rms = gaufit_dbl(ptr_start, size_xy, 1, range);
			
			for(double *ptr = ptr_start + size_xy; ptr --> ptr_start;) *ptr /= rms;
			Array_dbl_set(noise_spectrum, i, rms);
		}
	}
	
	return noise_spectrum;
}



/// @brief Local noise scaling within running window
///
/// Public method for dividing the data cube by the local noise
/// level in a running window throughout the entire data cube. The
/// size of the window and the size of the grid across which the
/// window is moved within the cube can be specified by the user.
/// If set to 0, default values will instead apply, with the grid
/// size being set to half the window size. Nearest-neighbour
/// interpolation will be used by default to fill the grid cells with
/// the noise measurement, unless `interpolation` is set to `true`,
/// in which case bilinear interpolation will instead be used for
/// positions in between the grid points. Once completed, the
/// method will return a DataCube object that contains the measured
/// noise values by which the cube was divided.
///
/// Alternatively, this method can also be used to measure and
/// subtract either the mean or the median across the specified window
/// by setting the argument `statistic` to either `NOISE_STAT_MEAN`
/// or `NOISE_STAT_MEDIAN`. This can be useful to remove DC offsets
/// or bandpass ripples from the data.
///
/// @param self         Object self-reference.
/// @param statistic    Statistic to use in noise measurement. Can
///                     be `NOISE_STAT_STD` for standard deviation,
///                     `NOISE_STAT_MAD` for median absolute deviation
///                     or `NOISE_STAT_GAUSS` for Gaussian fitting to
///                     the flux histogram. Alternatively, `NOISE_STAT_MEAN`
///                     or `NOISE_STAT_MEDIAN` can be used to measure and
///                     subtract the mean or median value.
/// @param range        Flux range to be used in noise measurement.
///                     Can be -1, 0 or +1 for negative range, full
///                     range or positive range, respectively. Note
///                     that this has no effect if the mean or median
///                     is to be calculated rather than the noise.
/// @param window_spat  Spatial window size in pixels; must be odd.
/// @param window_spec  Spectral window size in chan.; must be odd.
/// @param grid_spat    Spatial grid size in pixels; must be odd.
/// @param grid_spec    Spectral grid size in chan.; must be odd.
/// @param interpolate  If `true`, the noise or mean/median values will
///                     be interpolated in between the grid points.
///
/// @return Returns a data cube containing the measured noise values
///         (or the mean/median value if requested).

PUBLIC DataCube *DataCube_scale_noise_local(DataCube *self, const noise_stat statistic, const int range, size_t window_spat, size_t window_spec, size_t grid_spat, size_t grid_spec, const bool interpolate)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot run filter on integer array.");
	
	// Make window sizes integers >= 1
	window_spat = window_spat ? window_spat : 25;
	window_spec = window_spec ? window_spec : 15;
	
	// Ensure that window sizes are odd
	window_spat += 1 - window_spat % 2;
	window_spec += 1 - window_spec % 2;
	
	// Set grid to half the window size if not set
	grid_spat = grid_spat ? grid_spat : window_spat / 2;
	grid_spec = grid_spec ? grid_spec : window_spec / 2;
	
	// Make grid sizes integers >= 1
	grid_spat = grid_spat ? grid_spat : 1;
	grid_spec = grid_spec ? grid_spec : 1;
	
	// Ensure that grid sizes are odd
	grid_spat += 1 - grid_spat % 2;
	grid_spec += 1 - grid_spec % 2;
	
	// Print adopted grid and window sizes
	message("- Grid size:    %zu x %zu", grid_spat, grid_spec);
	message("- Window size:  %zu x %zu\n", window_spat, window_spec);
	
	// Divide grid and window sizes by 2 to get radii
	const size_t radius_grid_spat = grid_spat / 2;
	const size_t radius_grid_spec = grid_spec / 2;
	const size_t radius_window_spat = window_spat / 2;
	const size_t radius_window_spec = window_spec / 2;
	
	// Define starting point of grid
	const size_t grid_start_x = (self->axis_size[0] - grid_spat * (size_t)(ceil((double)(self->axis_size[0]) / (double)(grid_spat)) - 1.0)) / 2;
	const size_t grid_start_y = (self->axis_size[1] - grid_spat * (size_t)(ceil((double)(self->axis_size[1]) / (double)(grid_spat)) - 1.0)) / 2;
	const size_t grid_start_z = (self->axis_size[2] - grid_spec * (size_t)(ceil((double)(self->axis_size[2]) / (double)(grid_spec)) - 1.0)) / 2;
	
	// Define end point of grid
	const size_t grid_end_x = self->axis_size[0] - ((self->axis_size[0] - grid_start_x - 1) % grid_spat) - 1;
	const size_t grid_end_y = self->axis_size[1] - ((self->axis_size[1] - grid_start_y - 1) % grid_spat) - 1;
	const size_t grid_end_z = self->axis_size[2] - ((self->axis_size[2] - grid_start_z - 1) % grid_spec) - 1;
	
	// Create empty cube (filled with NaN) to hold noise/mean/median values
	DataCube *noiseCube = DataCube_blank(self->axis_size[0], self->axis_size[1], self->axis_size[2], self->data_type, self->verbosity);
	Header_copy_wcs(self->header, noiseCube->header);
	Header_copy_misc(self->header, noiseCube->header, true, true);
	DataCube_fill_flt(noiseCube, NAN);
	
	if(statistic == NOISE_STAT_STD || statistic == NOISE_STAT_MAD || statistic == NOISE_STAT_GAUSS) message("Measuring noise in running window.");
	else message("Measuring %s in running window.", statistic == NOISE_STAT_MEAN ? "mean" : "median");
	size_t progress = 0;
	const size_t progress_max = (grid_end_z - grid_start_z) / grid_spec;
	
	// Determine RMS/mean/median across window centred on grid cell
	#pragma omp parallel for schedule(static)
	for(size_t z = grid_start_z; z <= grid_end_z; z += grid_spec)
	{
		#pragma omp critical
		progress_bar("Progress: ", progress++, progress_max);
		
		for(size_t y = grid_start_y; y < self->axis_size[1]; y += grid_spat)
		{
			for(size_t x = grid_start_x; x < self->axis_size[0]; x += grid_spat)
			{
				// Determine extent of grid cell (inclusive of end point)
				const size_t grid[6] = {
					x < radius_grid_spat ? 0 : x - radius_grid_spat,
					x + radius_grid_spat >= self->axis_size[0] ? self->axis_size[0] - 1 : x + radius_grid_spat,
					y < radius_grid_spat ? 0 : y - radius_grid_spat,
					y + radius_grid_spat >= self->axis_size[1] ? self->axis_size[1] - 1 : y + radius_grid_spat,
					z < radius_grid_spec ? 0 : z - radius_grid_spec,
					z + radius_grid_spec >= self->axis_size[2] ? self->axis_size[2] - 1 : z + radius_grid_spec
				};
				
				// Determine extent of window (inclusive of end point)
				const size_t window[6] = {
					x < radius_window_spat ? 0 : x - radius_window_spat,
					x + radius_window_spat >= self->axis_size[0] ? self->axis_size[0] - 1 : x + radius_window_spat,
					y < radius_window_spat ? 0 : y - radius_window_spat,
					y + radius_window_spat >= self->axis_size[1] ? self->axis_size[1] - 1 : y + radius_window_spat,
					z < radius_window_spec ? 0 : z - radius_window_spec,
					z + radius_window_spec >= self->axis_size[2] ? self->axis_size[2] - 1 : z + radius_window_spec
				};
				
				// Create temporary array
				// NOTE: The use of float is faster and more memory-efficient than double.
				float *array = (float *)memory(MALLOC, (window[5] - window[4] + 1) * (window[3] - window[2] + 1) * (window[1] - window[0] + 1), sizeof(float));
				
				// Copy values from window into temporary array
				size_t counter = 0;
				for(size_t zz = window[4]; zz <= window[5]; ++zz)
				{
					for(size_t yy = window[2]; yy <= window[3]; ++yy)
					{
						for(size_t xx = window[0]; xx <= window[1]; ++xx)
						{
							const double value = DataCube_get_data_flt(self, xx, yy, zz);
							if(IS_NOT_NAN(value)) array[counter++] = value;
						}
					}
				}
				
				// Move on if not enough finite values found
				// NOTE: The threshold of 10 is somewhat arbitrary.
				if(counter < 10)
				{
					free(array);
					continue;
				}
				
				// Determine noise/mean/median in temporary array
				double rms;
				if(statistic == NOISE_STAT_STD) rms = std_dev_val_flt(array, counter, 0.0, 1, range);
				else if(statistic == NOISE_STAT_MAD) rms = MAD_TO_STD * mad_val_flt(array, counter, 0.0, 1, range);
				else if(statistic == NOISE_STAT_GAUSS) rms = gaufit_flt(array, counter, 1, range);
				else if(statistic == NOISE_STAT_MEAN) rms = mean_flt(array, counter);
				else rms = median_flt(array, counter, false);
				
				// Delete temporary array again
				free(array);
				
				// Fill entire grid cell with noise/mean/median value
				for(size_t zz = grid[4]; zz <= grid[5]; ++zz)
				{
					for(size_t yy = grid[2]; yy <= grid[3]; ++yy)
					{
						for(size_t xx = grid[0]; xx <= grid[1]; ++xx)
						{
							DataCube_set_data_flt(noiseCube, xx, yy, zz, rms);
						}
					}
				}
			}
		}
	}
	
	// Apply bilinear interpolation if requested
	if(interpolate && (grid_spat > 1 || grid_spec > 1))
	{
		if(statistic == NOISE_STAT_STD || statistic == NOISE_STAT_MAD || statistic == NOISE_STAT_GAUSS) message("Interpolating noise values.");
		else message("Interpolating %s values.", statistic == NOISE_STAT_MEAN ? "mean" : "median");
		
		// First interpolate along z-axis if necessary
		// (No need for multi-threading, as this is super-fast anyway)
		if(grid_spec > 1)
		{
			for(size_t y = grid_start_y; y <= grid_end_y; y += grid_spat)
			{
				progress_bar("Spectral: ", y - grid_start_y, grid_end_y - grid_start_y);
				
				for(size_t x = grid_start_x; x <= grid_end_x; x += grid_spat)
				{
					for(size_t z = grid_start_z; z < grid_end_z; z += grid_spec)
					{
						const size_t z0 = z;
						const size_t z2 = z + grid_spec;
						const double s0 = DataCube_get_data_flt(noiseCube, x, y, z0);
						const double s2 = DataCube_get_data_flt(noiseCube, x, y, z2);
						
						if(IS_NAN(s0) || IS_NAN(s2)) continue;
						
						for(size_t i = 1; i < grid_spec; ++i)
						{
							const size_t z1 = z0 + i;
							DataCube_set_data_flt(noiseCube, x, y, z1, s0 + (s2 - s0) * (double)(z1 - z0) / (double)(z2 - z0));
						}
					}
				}
			}
		}
		
		// Then interpolate across each spatial plane if necessary
		if(grid_spat > 1)
		{
			progress = 0;
			
			#pragma omp parallel for schedule(static)
			for(size_t z = 0; z < self->axis_size[2]; ++z)
			{
				#pragma omp critical
				progress_bar("Spatial:  ", progress++, self->axis_size[2] - 1);
				
				// Interpolate along y-axis
				for(size_t x = grid_start_x; x <= grid_end_x; x += grid_spat)
				{
					for(size_t y = grid_start_y; y < grid_end_y; y += grid_spat)
					{
						const size_t y0 = y;
						const size_t y2 = y + grid_spat;
						const double s0 = DataCube_get_data_flt(noiseCube, x, y0, z);
						const double s2 = DataCube_get_data_flt(noiseCube, x, y2, z);
						
						if(IS_NAN(s0) || IS_NAN(s2)) continue;
						
						for(size_t i = 1; i < grid_spat; ++i)
						{
							const size_t y1 = y0 + i;
							DataCube_set_data_flt(noiseCube, x, y1, z, s0 + (s2 - s0) * (double)(y1 - y0) / (double)(y2 - y0));
						}
					}
				}
				
				// Interpolate along x-axis
				for(size_t y = grid_start_y; y <= grid_end_y; ++y)
				{
					for(size_t x = grid_start_x; x < grid_end_x; x += grid_spat)
					{
						const size_t x0 = x;
						const size_t x2 = x + grid_spat;
						const double s0 = DataCube_get_data_flt(noiseCube, x0, y, z);
						const double s2 = DataCube_get_data_flt(noiseCube, x2, y, z);
						
						if(IS_NAN(s0) || IS_NAN(s2)) continue;
						
						for(size_t i = 1; i < grid_spat; ++i)
						{
							const size_t x1 = x0 + i;
							DataCube_set_data_flt(noiseCube, x1, y, z, s0 + (s2 - s0) * (double)(x1 - x0) / (double)(x2 - x0));
						}
					}
				}
			}
		}
	}
	
	if(statistic == NOISE_STAT_STD || statistic == NOISE_STAT_MAD || statistic == NOISE_STAT_GAUSS)
	{
		// Divide data cube by noise cube
		if(self->data_type == -32)
		{
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				float *ptr_data  = (float *)(self->data) + i;
				float *ptr_noise = (float *)(noiseCube->data) + i;
				
				if(*ptr_noise > 0.0) *ptr_data /= *ptr_noise;
				else *ptr_data = NAN;
			}
		}
		else
		{
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				double *ptr_data  = (double *)(self->data) + i;
				double *ptr_noise = (double *)(noiseCube->data) + i;
				
				if(*ptr_noise > 0.0) *ptr_data /= *ptr_noise;
				else *ptr_data = NAN;
			}
		}
	}
	else
	{
		// Subtract mean/median value
		if(self->data_type == -32)
		{
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				float *ptr_data  = (float *)(self->data) + i;
				float *ptr_noise = (float *)(noiseCube->data) + i;
				*ptr_data -= *ptr_noise;
			}
		}
		else
		{
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < self->data_size; ++i)
			{
				double *ptr_data  = (double *)(self->data) + i;
				double *ptr_noise = (double *)(noiseCube->data) + i;
				*ptr_data -= *ptr_noise;
			}
		}
	}
	
	return noiseCube;
}



/// @brief Apply boxcar filter to spectral axis
///
/// Public method for convolving each spectrum of the data cube
/// with a boxcar filter of size `2 * radius + 1`.
///
/// @note The algorithm is `NaN`-safe by setting all `NaN` values
///       to 0 prior to filtering. Any pixel outside of the cube's
///       spectral range is also assumed to be 0.
///
/// @param self    Object self-reference.
/// @param radius  Filter radius in channels.

PUBLIC void DataCube_boxcar_filter(DataCube *self, size_t radius)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot run boxcar filter on integer array.");
	if(radius < 1) return;
	
	if(self->data_type == -32)
	{
		// Single-precision floating-point type
		#pragma omp parallel
		{
			// Allocate memory for a single spectrum
			float  *spectrum = (float *)memory(MALLOC, self->axis_size[2], sizeof(float));
			
			// Request memory for boxcar filter to operate on
			float  *data_box = (float *) memory(MALLOC, self->axis_size[2] + 2 * radius, sizeof(float));
			
			#pragma omp for schedule(static)
			for(size_t y = 0; y < self->axis_size[1]; ++y)
			{
				for(size_t x = 0; x < self->axis_size[0]; ++x)
				{
					// Extract spectrum
					for(size_t z = self->axis_size[2]; z--;) *(spectrum + z) = DataCube_get_data_flt(self, x, y, z);
					
					// Apply filter
					filter_boxcar_1d_flt(spectrum, data_box, self->axis_size[2], radius);
					
					// Copy filtered spectrum back into array
					for(size_t z = self->axis_size[2]; z--;) DataCube_set_data_flt(self, x, y, z, *(spectrum + z));
				}
			}
			
			// Release memory
			free(spectrum);
			free(data_box);
		}
	}
	else
	{
		// Double-precision floating-point type
		#pragma omp parallel
		{
			// Allocate memory for a single spectrum
			double *spectrum = (double *)memory(MALLOC, self->axis_size[2], sizeof(double));
			
			// Request memory for boxcar filter to operate on
			double *data_box = (double *)memory(MALLOC, self->axis_size[2] + 2 * radius, sizeof(double));
			
			#pragma omp for schedule(static)
			for(size_t y = 0; y < self->axis_size[1]; ++y)
			{
				for(size_t x = 0; x < self->axis_size[0]; ++x)
				{
					// Extract spectrum
					for(size_t z = self->axis_size[2]; z--;) *(spectrum + z) = DataCube_get_data_flt(self, x, y, z);
					
					// Apply filter
					filter_boxcar_1d_dbl(spectrum, data_box, self->axis_size[2], radius);
					
					// Copy filtered spectrum back into array
					for(size_t z = self->axis_size[2]; z--;) DataCube_set_data_flt(self, x, y, z, *(spectrum + z));
				}
			}
			
			// Release memory
			free(spectrum);
			free(data_box);
		}
	}
	
	return;
}



/// @brief Apply 2D Gaussian filter to spatial planes
///
/// Public method for convolving each spatial image plane (x-y) of
/// the data cube with a Gaussian function of standard deviation
/// sigma. The Gaussian convolution is approximated through a set
/// of 1D boxcar filters, which makes the algorithm extremely fast.
/// Limitations from this approach are that the resulting convolution
/// kernel is only an approximation of a Gaussian (although a fairly
/// accurate one) and the value of `sigma` can only be approximated
/// (typically within +/- 0.2 sigma) and must be at least 1.5 pixels.
///
/// @note The algorithm is `NaN`-safe by setting all `NaN` values to 0.
///       Any pixel outside of the image boundaries is also assumed
///       to be 0.
///
/// @param self   Object self-reference.
/// @param sigma  Standard deviation of the Gaussian in pixels.

PUBLIC void DataCube_gaussian_filter(DataCube *self, const double sigma)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot run boxcar filter on integer array.");
	
	// Set up parameters required for boxcar filter
	// NOTE: We don't need to extract a copy of each image plane, as
	///     x-y planes are contiguous in memory.
	size_t n_iter;
	size_t filter_radius;
	optimal_filter_size_dbl(sigma, &filter_radius, &n_iter);
	const size_t size = self->axis_size[0] * self->axis_size[1] * self->word_size;
	
	if(self->data_type == -32)
	{
		#pragma omp parallel
		{
			// Memory for one column
			float  *column   = (float *)memory(MALLOC, self->axis_size[1], sizeof(float));
			
			// Memory for boxcar filter to operate on
			float  *data_row = (float *)memory(MALLOC, self->axis_size[0] + 2 * filter_radius, sizeof(float));
			float  *data_col = (float *)memory(MALLOC, self->axis_size[1] + 2 * filter_radius, sizeof(float));
			
			// Apply filter
			#pragma omp for schedule(static)
			for(char *ptr = self->data; ptr < self->data + self->data_size * self->word_size; ptr += size)
			{
				filter_gauss_2d_flt((float *)ptr, column, data_row, data_col, self->axis_size[0], self->axis_size[1], n_iter, filter_radius);
			}
			
			// Release memory
			free(data_row);
			free(data_col);
			free(column);
		}
	}
	else
	{
		#pragma omp parallel
		{
			// Memory for boxcar filter to operate on
			double *data_row = (double *)memory(MALLOC, self->axis_size[0] + 2 * filter_radius, sizeof(double));
			double *data_col = (double *)memory(MALLOC, self->axis_size[1] + 2 * filter_radius, sizeof(double));
			
			// Memory for one column
			double *column   = (double *)memory(MALLOC, self->axis_size[1], sizeof(double));
			
			// Apply filter
			#pragma omp for schedule(static)
			for(char *ptr = self->data; ptr < self->data + self->data_size * self->word_size; ptr += size)
			{
				filter_gauss_2d_dbl((double *)ptr, column, data_row, data_col, self->axis_size[0], self->axis_size[1], n_iter, filter_radius);
			}
			
			// Release memory
			free(data_row);
			free(data_col);
			free(column);
		}
	}
	
	return;
}



/// @brief Subtract residual continuum emission
///
/// Public method for fitting and subtracting a polynomial from the
/// spectrum at each spatial position of the data cube. Currently,
/// `order` = 0 (constant offset) and 1 (offset + linear slope) are
/// implemented and supported.
/// The algorithm works by subtracting the spectrum shifted by
/// `-shift` from the same spectrum shifted by `+shift`. It then uses
/// a robust algorithm for measuring the noise in the shifted and
/// subtracted spectrum and flags all channels in the original
/// spectrum where the flux density exceeds threshold times the
/// noise level. A certain amount of padding can be applied by setting
/// the padding parameter to > 0. In a last step, a polynomial of
/// order 0 or 1 is fitted and subtracted from all channels in the
/// original data cube.
/// For this algorithm to work correctly, the data must be a 3D
/// cube with a sufficiently large number of channels. In addition,
/// the continuum residual must be a simple offset + slope, and
/// higher-order variation cannot be handled at the moment. It is
/// also crucial that any emission lines in the data cube do not
/// cover more than about 20% of the spectral band, as otherwise
/// their presence may start to influence the fit.
///
/// @param self       Object self-reference.
/// @param order      Order of polynomial fit (0 or 1).
/// @param shift      Amount by which to shift and subtract spectrum.
/// @param padding    Padding around flagged channels with emission.
/// @param threshold  Threshold for flagging of emission.

PUBLIC void DataCube_contsub(DataCube *self, unsigned int order, size_t shift, const size_t padding, double threshold)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cannot subtract continuum from integer data.");
	ensure(self->axis_size[2] > 5 * shift, ERR_USER_INPUT, "Continuum subtraction requires 3D data cube with > %zu channels.", 5 * shift);
	
	if(order > 1)
	{
		order = 1;
		warning("Adjusting value of polynomial order to 1.");
	}
	
	if(shift < 1)
	{
		shift = 1;
		warning("Adjusting value of shift to 1.");
	}
	
	if(threshold < 0.0)
	{
		threshold = -threshold;
		warning("Adjusting threshold to be positive.");
	}
	
	// Basic settings
	const size_t nx = self->axis_size[0];
	const size_t ny = self->axis_size[1];
	const size_t nz = self->axis_size[2];
	const size_t nxy = nx * ny;
	size_t progress = 0;
	
	#pragma omp parallel
	{
		double *spectrum     = (double *)memory(MALLOC, nz, sizeof(double));
		double *spectrum_tmp = (double *)memory(MALLOC, nz, sizeof(double));
		
		// Loop over all spatial pixels
		#pragma omp for schedule(static)
		for(size_t y = 0; y < ny; ++y)
		{
			#pragma omp critical
			progress_bar("Progress: ", progress++, ny - 1);
			
			for(size_t x = 0; x < nx; ++x)
			{
				// Extract spectrum
				if(self->data_type == -32)
				{
					// 32-bit float
					float *ptr = (float *)(self->data) + x + nx * y;
					for(size_t i = 0; i < nz; ++i)
					{
						spectrum[i] = *ptr;
						ptr += nxy;
					}
				}
				else
				{
					// 64-bit float
					double *ptr = (double *)(self->data) + x + nx * y;
					for(size_t i = 0; i < nz; ++i)
					{
						spectrum[i] = *ptr;
						ptr += nxy;
					}
				}
				
				// Shift and subtract spectrum from itself
				for(size_t i = 0; i < nz; ++i)
				{
					spectrum_tmp[i] = (i < shift || i >= nz - shift) ? NAN : spectrum[i - shift] - spectrum[i + shift];
				}
				
				// Robust noise measurement
				const double rms = threshold * robust_noise_2_dbl(spectrum_tmp, nz);
				
				// Mask everything > rms
				for(size_t i = 0; i < nz; ++i)
				{
					if(fabs(spectrum_tmp[i]) > rms)
					{
						const size_t j_min = (i > padding) ? i - padding : 0;
						const size_t j_max = (i + padding < nz) ? i + padding : nz - 1;
						for(size_t j = j_min; j <= j_max; ++j) spectrum[j] = NAN;
					}
				}
				
				// Measure means
				double x_mean = 0.0;
				double y_mean = 0.0;
				size_t counter = 0;
				
				for(size_t i = 0; i < nz; ++i)
				{
					if(IS_NOT_NAN(spectrum[i]))
					{
						x_mean += i;
						y_mean += spectrum[i];
						++counter;
					}
				}
				
				if(counter == 0) continue;  // Cannot fit, as nothing left after filtering
				x_mean /= counter;
				y_mean /= counter;
				
				if(order)
				{
					// Fit and subtract 1st-order polynomial
					double alpha = 0.0;
					double beta = 0.0;
					
					for(size_t i = 0; i < nz; ++i)
					{
						if(IS_NOT_NAN(spectrum[i]))
						{
							alpha += (x_mean - i) * (x_mean - i);
							beta  += (x_mean - i) * (y_mean - spectrum[i]);
						}
					}
					
					if(alpha == 0)
					{
						// Cannot fit for some reason
						warning("Polynomial fit failed at position (%zu, %zu).", x, y);
						continue;
					}
					
					beta /= alpha;
					alpha = y_mean - beta * x_mean;
					
					if(self->data_type == -32)
					{
						// 32-bit float
						float *ptr = (float *)(self->data) + x + nx * y;
						for(size_t i = 0; i < nz; ++i)
						{
							*ptr -= (alpha + beta * i);
							ptr += nxy;
						}
					}
					else
					{
						// 64-bit float
						double *ptr = (double *)(self->data) + x + nx * y;
						for(size_t i = 0; i < nz; ++i)
						{
							*ptr -= (alpha + beta * i);
							ptr += nxy;
						}
					}
				}
				else
				{
					// Subtract mean
					if(self->data_type == -32)
					{
						// 32-bit float
						float *ptr = (float *)(self->data) + x + nx * y;
						for(size_t i = 0; i < nz; ++i)
						{
							*ptr -= y_mean;
							ptr += nxy;
						}
					}
					else
					{
						// 64-bit float
						double *ptr = (double *)(self->data) + x + nx * y;
						for(size_t i = 0; i < nz; ++i)
						{
							*ptr -= y_mean;
							ptr += nxy;
						}
					}
				}
			}
		}
		
		// Clean-up
		free(spectrum);
		free(spectrum_tmp);
	}
	
	return;
}



/// @brief Mask pixels of abs(value) > threshold
///
/// Public method for setting pixels in the mask cube to 1 when
/// their absolute value in the data cube is greater than the
/// specified threshold. Similar to DataCube_mask_8(), but for
/// mask cubes of any integer type.
///
/// @param self       Object self-reference.
/// @param maskCube   Pointer to mask cube.
/// @param threshold  Flux threshold for masking operation.

PUBLIC void DataCube_mask(const DataCube *self, DataCube *maskCube, const double threshold)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 8 || maskCube->data_type == 16 || maskCube->data_type == 32 || maskCube->data_type == 64, ERR_USER_INPUT, "Mask cube must be of integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(threshold > 0.0, ERR_USER_INPUT, "Negative threshold provided.");
	
	// Declaration of variables
	char *ptr_data = self->data + self->data_size * self->word_size;
	char *ptr_mask = maskCube->data + maskCube->data_size * maskCube->word_size;
	const float  thresh_flt_pos =  threshold;
	const float  thresh_flt_neg = -threshold;
	const double thresh_dbl_pos =  threshold;
	const double thresh_dbl_neg = -threshold;
	float  *ptr_flt;
	double *ptr_dbl;
	
	while(ptr_data > self->data)
	{
		ptr_data -= self->word_size;
		ptr_mask -= maskCube->word_size;
		
		if(self->data_type == -32)
		{
			ptr_flt = (float *)ptr_data;
			if(*ptr_flt > thresh_flt_pos || *ptr_flt < thresh_flt_neg)
			{
				if     (maskCube->data_type ==  8) *((uint8_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 16) *((int16_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 32) *((int32_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 64) *((int64_t *)ptr_mask) = 1;
			}
		}
		else
		{
			ptr_dbl = (double *)ptr_data;
			if(*ptr_dbl > thresh_dbl_pos || *ptr_dbl < thresh_dbl_neg)
			{
				if     (maskCube->data_type ==  8) *((uint8_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 16) *((int16_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 32) *((int32_t *)ptr_mask) = 1;
				else if(maskCube->data_type == 64) *((int64_t *)ptr_mask) = 1;
			}
		}
	}
	
	return;
}



/// @brief Mask pixels of abs(value) > threshold
///
/// Public method for setting pixels in the mask cube to `value`
/// when their absolute value in the data cube is greater than the
/// specified threshold. Similar to DataCube_mask(), but for 8-bit
/// masks.
///
/// @param self       Object self-reference.
/// @param maskCube   Pointer to mask cube.
/// @param threshold  Flux threshold for masking operation.
/// @param value      Value to set the mask to.

PUBLIC void DataCube_mask_8(const DataCube *self, DataCube *maskCube, const double threshold, const uint8_t value)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 8, ERR_USER_INPUT, "Mask cube must be of 8-bit integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(threshold > 0.0, ERR_USER_INPUT, "Threshold must be positive.");
	
	uint8_t *ptr_mask = (uint8_t *)(maskCube->data);
	
	if(self->data_type == -32)
	{
		const float *ptr_data = (float *)(self->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(fabs(*(ptr_data + i)) > threshold) *(ptr_mask + i) = value;
		}
	}
	else
	{
		const double *ptr_data = (double *)(self->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(fabs(*(ptr_data + i)) > threshold) *(ptr_mask + i) = value;
		}
	}
	
	return;
}



/// @brief Set masked pixels to constant value
///
/// Public method for replacing the values of all pixels in the
/// data cube that are non-zero in the mask cube to their signum
/// multiplied by the specified value. Same as DataCube_set_masked_8(),
/// but for masks of any integer type.
///
/// @param self      Object self-reference.
/// @param maskCube  Pointer to mask cube.
/// @param value     Flux value to replace pixels with.

PUBLIC void DataCube_set_masked(DataCube *self, const DataCube *maskCube, const double value)
{
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 8 || maskCube->data_type == 16 || maskCube->data_type == 32 || maskCube->data_type == 64, ERR_USER_INPUT, "Mask cube must be of integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	
	char *ptr_data = self->data + self->data_size * self->word_size;
	char *ptr_mask = maskCube->data + maskCube->data_size * maskCube->word_size;
	const float  value_flt =  value;
	float *ptr_flt;
	double *ptr_dbl;
	const int mask_type = maskCube->data_type;
	
	while(ptr_data > self->data)
	{
		ptr_data -= self->word_size;
		ptr_mask -= maskCube->word_size;
		
		if((mask_type == 8 && *((uint8_t *)ptr_mask)) || (mask_type == 16 && *((int16_t *)ptr_mask)) || (mask_type == 32 && *((int32_t *)ptr_mask)) || (mask_type == 64 && *((int64_t *)ptr_mask)))
		{
			if(self->data_type == -32)
			{
				ptr_flt  = (float *)ptr_data;
				*ptr_flt = copysign(value_flt, *ptr_flt);
			}
			else
			{
				ptr_dbl  = (double *)ptr_data;
				*ptr_dbl = copysign(value, *ptr_dbl);
			}
		}
	}
	
	return;
}



/// @brief Set masked pixels to constant value
///
/// Public method for replacing the values of all pixels in the
/// data cube that are non-zero in the mask cube to their signum
/// multiplied by the specified value. Same as DataCube_set_masked(),
/// but for 8-bit masks.
///
/// @param self      Object self-reference.
/// @param maskCube  Pointer to mask cube.
/// @param value     Flux value to replace pixels with.

PUBLIC void DataCube_set_masked_8(DataCube *self, const DataCube *maskCube, const double value)
{
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 8, ERR_USER_INPUT, "Mask cube must be of 8-bit integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	
	const uint8_t *ptr_mask = (uint8_t *)(maskCube->data);
	
	if(self->data_type == -32)
	{
		float *ptr_data = (float *)(self->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptr_mask + i)) *(ptr_data + i) = copysign(value, *(ptr_data + i));
		}
	}
	else
	{
		double *ptr_data = (double *)(self->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptr_mask + i)) *(ptr_data + i) = copysign(value, *(ptr_data + i));
		}
	}
	
	return;
}



/// @brief Replace masked pixels with the specified value
///
/// Public method for replacing the values of all pixels in the
/// mask cube that are non-zero with the specified value. The mask
/// cube must be of 32-bit integer type.
///
/// @param self   Object self-reference.
/// @param value  Mask value to replace pixels with.

PUBLIC void DataCube_reset_mask_32(DataCube *self, const int32_t value)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == 32, ERR_USER_INPUT, "Mask cube must be of 32-bit integer type.");
	
	#pragma omp parallel for schedule(static)
	for(int32_t *ptr = (int32_t *)(self->data); ptr < (int32_t *)(self->data) + self->data_size; ++ptr)
	{
		if(*ptr) *ptr = value;
	}
	
	return;
}



/// @brief Remove unreliable sources from mask and relabel remaining ones
///
/// Public method for removing unreliable sources from the mask.
/// This is done by comparing the pixel values with a list of old
/// and new source labels and, if present in that list, replace the
/// pixel value with its new label. Pixel values not present in the
/// list will be discarded by setting them to 0. If an empty filter
/// is supplied, a warning message appears and no filtering will be
/// done.
///
/// @param self    Object self-reference.
/// @param filter  Map object with old and new label pairs of all
///                reliable sources.

PUBLIC void DataCube_filter_mask_32(DataCube *self, const Map *filter)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == 32, ERR_USER_INPUT, "Mask cube must be of 32-bit integer type.");
	
	check_null(filter);
	if(Map_get_size(filter) == 0)
	{
		warning("Empty filter provided. Cannot filter mask.");
		return;
	}
	
	#pragma omp parallel for schedule(static)
	for(int32_t *ptr = (int32_t *)(self->data); ptr < (int32_t *)(self->data) + self->data_size; ++ptr)
	{
		if(*ptr > 0)
		{
			if(Map_key_exists(filter, *ptr)) *ptr = Map_get_value(filter, *ptr);
			else *ptr = 0;
		}
	}
	
	return;
}



/// @brief Copy masked pixels from any integer mask to 32-bit mask
///
/// Public method for setting all of the pixels that are not equal
/// to zero in the integer source mask to the specified value in
/// the signed 32-bit target mask. This can be used to copy masked
/// pixels from any integer mask to a 32-bit mask.
///
/// @param self    32-bit target mask.
/// @param source  integer source mask.
/// @param value   Mask value to set in target mask.
///
/// @return Number of masked pixels.

PUBLIC size_t DataCube_copy_mask_32(DataCube *self, const DataCube *source, const int32_t value)
{
	// Sanity checks
	check_null(self);
	check_null(source);
	ensure(self->data_type == 32, ERR_USER_INPUT, "Target mask cube must be of 32-bit integer type.");
	ensure(source->data_type == 8 || source->data_type == 16 || source->data_type == 32 || source->data_type == 64, ERR_USER_INPUT, "Source mask cube must be of integer type.");
	
	int32_t *ptrTarget = (int32_t *)(self->data);
	size_t counter = 0;
	
	if(source->data_type == 8)
	{
		uint8_t *ptrSource = (uint8_t *)(source->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptrSource + i) != 0)
			{
				*(ptrTarget + i) = value;
				#pragma omp atomic update
				++counter;
			}
		}
	}
	else if(source->data_type == 16)
	{
		int16_t *ptrSource = (int16_t *)(source->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptrSource + i) != 0)
			{
				*(ptrTarget + i) = value;
				#pragma omp atomic update
				++counter;
			}
		}
	}
	else if(source->data_type == 32)
	{
		int32_t *ptrSource = (int32_t *)(source->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptrSource + i) != 0)
			{
				*(ptrTarget + i) = value;
				#pragma omp atomic update
				++counter;
			}
		}
	}
	else
	{
		int64_t *ptrSource = (int64_t *)(source->data);
		
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size; ++i)
		{
			if(*(ptrSource + i) != 0)
			{
				*(ptrTarget + i) = value;
				#pragma omp atomic update
				++counter;
			}
		}
	}
	
	return counter;
}



/// @brief Grow mask and update basic source parameters
///
/// Private method for growing the mask of the source identified by
/// `src_id` radially outwards using a circular kernel of the specified
/// radius (in pixels). The algorithm will update the relevant source
/// parameters (which must be provided as pointers holding the old
/// values prior to mask growth), including fluxes, flags and source
/// bounding boxes. In addition, the user can specify a mask value to
/// which all newly added pixels will be set in the source mask.
///
/// @param self        Data cube.
/// @param mask        Mask cube.
/// @param src_id      ID of source to be grown.
/// @param radius      Radius by which to grow (in pixels).
/// @param mask_value  Value to set grown pixels to in the mask.
/// @param f_sum       Summed flux density prior to mask growth.
/// @param f_min       Minimum flux density prior to mask growth.
/// @param f_max       Maximum flux density prior to mask growth.
/// @param n_pix       Number of pixels in mask prior to mask growth.
/// @param flag        Source flag prior to mask growth.
/// @param x_min       Source bounding box prior to mask growth.
/// @param x_max       Source bounding box prior to mask growth.
/// @param y_min       Source bounding box prior to mask growth.
/// @param y_max       Source bounding box prior to mask growth.
/// @param z_min       Source bounding box prior to mask growth.
/// @param z_max       Source bounding box prior to mask growth.
///
/// @note As the algorithm only grows the mask in x and y, the
///       bounding box in z will not change and will need to be
///       provided by value rather than by reference.

PRIVATE void DataCube_grow_mask_xy(const DataCube *self, DataCube *mask, const long int src_id, const size_t radius, const long int mask_value, double *f_sum, double *f_min, double *f_max, size_t *n_pix, long int *flag, size_t *x_min, size_t *x_max, size_t *y_min, size_t *y_max, const size_t z_min, const size_t z_max)
{
	// Set loop boundaries from bounding box
	const size_t x1 = *x_min;
	const size_t x2 = *x_max;
	const size_t y1 = *y_min;
	const size_t y2 = *y_max;
	
	// Additional variables
	const double radius2 = (double)(radius * radius);
	size_t xx_min = 0;
	size_t xx_max = 0;
	size_t yy_min = 0;
	size_t yy_max = 0;
	
	// Loop over source bounding box
	for(size_t z = z_min; z <= z_max; ++z)
	{
		for(size_t y = y1; y <= y2; ++y)
		{
			for(size_t x = x1; x <= x2; ++x)
			{
				const long int id = DataCube_get_data_int(mask, x, y, z);
				
				if(id == src_id)
				{
					// Define boundaries of box to be processed
					if(x < radius) xx_min = 0, *flag |= 1L;
					else xx_min = x - radius;
					
					if(x + radius >= self->axis_size[0]) xx_max = self->axis_size[0] - 1, *flag |= 1L;
					else xx_max = x + radius;
					
					if(y < radius) yy_min = 0, *flag |= 1L;
					else yy_min = y - radius;
					
					if(y + radius >= self->axis_size[1]) yy_max = self->axis_size[1] - 1, *flag |= 1L;
					else yy_max = y + radius;
					
					// Loop over that box
					for(size_t yy = yy_min; yy <= yy_max; ++yy)
					{
						for(size_t xx = xx_min; xx <= xx_max; ++xx)
						{
							const double tmp1 = (double)xx - (double)x;
							const double tmp2 = (double)yy - (double)y;
							if(tmp1 * tmp1 + tmp2 * tmp2 > radius2) continue;
							
							const long int id_new = DataCube_get_data_int(mask, xx, yy, z);
							if(id_new == 0)
							{
								const double value = DataCube_get_data_flt(self, xx, yy, z);
								
								if(IS_NOT_NAN(value))
								{
									DataCube_set_data_int(mask, xx, yy, z, -1);
									*f_sum += value;
									if(value < *f_min) *f_min = value;
									if(value > *f_max) *f_max = value;
									if(xx < *x_min) *x_min = xx;
									if(xx > *x_max) *x_max = xx;
									if(yy < *y_min) *y_min = yy;
									if(yy > *y_max) *y_max = yy;
									*n_pix += 1;
								}
								else *flag |= 4L;
							}
							else if(id_new > 0 && id_new != src_id) *flag |= 8L;
						} // END loop over xx
					} // END loop over yy
				}
			} // END loop over x
		} // END loop over y
	} // END loop over z
	
	// Replace mask values as per user request
	for(size_t z = z_min; z <= z_max; ++z)
	{
		for(size_t y = *y_min; y <= *y_max; ++y)
		{
			for(size_t x = *x_min; x <= *x_max; ++x)
			{
				if(DataCube_get_data_int(mask, x, y, z) == -1) DataCube_set_data_int(mask, x, y, z, mask_value);
			}
		}
	}
	
	return;
}



/// @brief Mask dilation in the spatial plane
///
/// Public method for dilating the masks of all sources found in
/// catalogue in the spatial plane. Dilation will occur iteratively
/// until either the maximum number of iterations is reached or
/// the relative increase in source flux drops below the specified
/// threshold. If the threshold is negative, then the mask will be
/// dilated by the maximum number of iterations regardless of the
/// resulting flux change.
///
/// Mask dilation will work correctly for sources with positive or
/// negative flux; in the latter case the integrated flux is expected
/// to decrease in each iteration. The source parameters in the
/// catalogue will be updated with the new, dilated values.
/// Dilation will be carried out using a circular mask and iteratively
/// progress outwards by increasing the radius of the mask by
/// 1 pixel in each iteration. The source mask should therefore
/// approach a circle for a large number of iterations.
///
/// @param self       Data cube.
/// @param mask       Mask cube.
/// @param cat        Source catalogue.
/// @param iter_max   Maximum number of iterations.
/// @param threshold  Threshold for relative flux increase.

PUBLIC void DataCube_dilate_mask_xy(const DataCube *self, DataCube *mask, Catalog *cat, const size_t iter_max, const double threshold)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	check_null(cat);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(mask->data_type > 8, ERR_USER_INPUT, "Mask must be of signed integer type.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(iter_max < self->axis_size[0] || iter_max < self->axis_size[1], ERR_USER_INPUT, "Maximum number of iterations exceeds spatial axis size.");
	
	// Determine catalogue size
	const size_t cat_size = Catalog_get_size(cat);
	if(cat_size == 0)
	{
		warning("No sources in catalogue; skipping mask dilation.");
		return;
	}
	
	// Loop over all sources in catalogue
	for(size_t i = 0; i < cat_size; ++i)
	{
		// Extract source
		Source *src = Catalog_get_source(cat, i);
		
		// Get source ID
		const long int src_id = Source_get_par_by_name_int(src, "id");
		ensure(src_id, ERR_USER_INPUT, "Source ID missing from catalogue; mask dilation failed.");
		
		// Get source bounding box
		const size_t x_min = Source_get_par_by_name_int(src, "x_min");
		const size_t x_max = Source_get_par_by_name_int(src, "x_max");
		const size_t y_min = Source_get_par_by_name_int(src, "y_min");
		const size_t y_max = Source_get_par_by_name_int(src, "y_max");
		const size_t z_min = Source_get_par_by_name_int(src, "z_min");
		const size_t z_max = Source_get_par_by_name_int(src, "z_max");
		ensure(x_min <= x_max && y_min <= y_max && z_min <= z_max, ERR_INDEX_RANGE, "Illegal source bounding box: min > max!");
		ensure(x_max < self->axis_size[0] && y_max < self->axis_size[1] && z_max < self->axis_size[2], ERR_INDEX_RANGE, "Source bounding box outside data cube boundaries.");
		
		// Get fluxes and check if negative
		const double f_sum = Source_get_par_by_name_flt(src, "f_sum");
		const double f_min = Source_get_par_by_name_flt(src, "f_min");
		const double f_max = Source_get_par_by_name_flt(src, "f_max");
		const bool is_negative = (f_sum < 0.0);
		
		// Get other relevant source parameters
		const size_t n_pix = Source_get_par_by_name_int(src, "n_pix");
		const long int flag = Source_get_par_by_name_int(src, "flag");
		
		// If threshold negative, simply dilate and move on
		if(threshold < 0.0)
		{
			double f_sum_new = f_sum;
			double f_min_new = f_min;
			double f_max_new = f_max;
			long int flag_new = flag;
			size_t n_pix_new = n_pix;
			size_t x_min_new = x_min;
			size_t x_max_new = x_max;
			size_t y_min_new = y_min;
			size_t y_max_new = y_max;
			
			// Dilate by iter_max
			DataCube_grow_mask_xy(self, mask, src_id, iter_max, src_id, &f_sum_new, &f_min_new, &f_max_new, &n_pix_new, &flag_new, &x_min_new, &x_max_new, &y_min_new, &y_max_new, z_min, z_max);
			
			// Update source parameters with new values
			Source_set_par_flt(src, "f_min", f_min_new, NULL, NULL);
			Source_set_par_flt(src, "f_max", f_max_new, NULL, NULL);
			Source_set_par_flt(src, "f_sum", f_sum_new, NULL, NULL);
			Source_set_par_int(src, "x_min", x_min_new, NULL, NULL);
			Source_set_par_int(src, "x_max", x_max_new, NULL, NULL);
			Source_set_par_int(src, "y_min", y_min_new, NULL, NULL);
			Source_set_par_int(src, "y_max", y_max_new, NULL, NULL);
			Source_set_par_int(src, "n_pix", n_pix_new, NULL, NULL);
			Source_set_par_int(src, "flag",  flag_new,  NULL, NULL);
			
			// Update progress bar
			progress_bar("Progress: ", i + 1, cat_size);
			
			// Next source
			continue;
		}
		else
		{
			// Iterative dilation required
			message_verb(self->verbosity, "Source %zu", i);
			
			// Iterate
			size_t iter = 1;
			double f_sum_old = f_sum;
			while(iter <= iter_max)
			{
				double f_sum_new = f_sum;
				double f_min_new = f_min;
				double f_max_new = f_max;
				long int flag_new = flag;
				size_t n_pix_new = n_pix;
				size_t x_min_new = x_min;
				size_t x_max_new = x_max;
				size_t y_min_new = y_min;
				size_t y_max_new = y_max;
				
				// Dilate mask by radius 'iter'
				DataCube_grow_mask_xy(self, mask, src_id, iter, 0, &f_sum_new, &f_min_new, &f_max_new, &n_pix_new, &flag_new, &x_min_new, &x_max_new, &y_min_new, &y_max_new, z_min, z_max);
				
				// Print information
				message_verb(self->verbosity, " - Iteration %zu: df = %.3f (%.3f%%)", iter, f_sum_new - f_sum_old, 100.0 * (f_sum_new - f_sum_old) / f_sum_old);
				
				// Check if flux increased within boundaries
				if((is_negative && (f_sum_new - f_sum_old) < threshold * f_sum_old) || (!is_negative && (f_sum_new - f_sum_old) > threshold * f_sum_old))
				{
					// Mask should be grown further
					f_sum_old = f_sum_new;
					++iter;
					continue;
				}
				else
				{
					// No significant improvement; stop iterating
					break;
				}
			} // END iteration loop
			
			// Apply final dilation if required
			if(iter > 1)
			{
				double f_sum_new = f_sum;
				double f_min_new = f_min;
				double f_max_new = f_max;
				long int flag_new = flag;
				size_t n_pix_new = n_pix;
				size_t x_min_new = x_min;
				size_t x_max_new = x_max;
				size_t y_min_new = y_min;
				size_t y_max_new = y_max;
				
				DataCube_grow_mask_xy(self, mask, src_id, iter - 1, src_id, &f_sum_new, &f_min_new, &f_max_new, &n_pix_new, &flag_new, &x_min_new, &x_max_new, &y_min_new, &y_max_new, z_min, z_max);
				
				// Update source parameters with new values
				Source_set_par_flt(src, "f_min", f_min_new, NULL, NULL);
				Source_set_par_flt(src, "f_max", f_max_new, NULL, NULL);
				Source_set_par_flt(src, "f_sum", f_sum_new, NULL, NULL);
				Source_set_par_int(src, "x_min", x_min_new, NULL, NULL);
				Source_set_par_int(src, "x_max", x_max_new, NULL, NULL);
				Source_set_par_int(src, "y_min", y_min_new, NULL, NULL);
				Source_set_par_int(src, "y_max", y_max_new, NULL, NULL);
				Source_set_par_int(src, "n_pix", n_pix_new, NULL, NULL);
				Source_set_par_int(src, "flag",  flag_new,  NULL, NULL);
			}
			
			// Update progress bar
			if(!self->verbosity) progress_bar("Progress: ", i + 1, cat_size);
		}
	}  // END source loop
	
	return;
}



/// @brief Mask dilation along spectral axis
///
/// Public method for dilating the masks of all sources found in
/// catalogue along the spectral axis. Dilation will occur iteratively
/// until either the maximum number of iterations is reached or the
/// relative increase in source flux drops below the specified
/// threshold.
///
/// Mask dilation will work correctly for sources with positive or
/// negative flux; in the latter case the integrated flux is expected
/// to decrease in each iteration. The source parameters in the
/// catalogue will be updated with the new, dilated values.
/// Dilation will progress by 1 channel per iteration in the
/// directions directly adjacent to a pixel along the spectral axis.
///
/// @param self       Data cube.
/// @param mask       Mask cube.
/// @param cat        Source catalogue.
/// @param iter_max   Maximum number of iterations.
/// @param threshold  Threshold for relative flux increase.

PUBLIC void DataCube_dilate_mask_z(const DataCube *self, DataCube *mask, Catalog *cat, const size_t iter_max, const double threshold)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	check_null(cat);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	ensure(mask->data_type > 8, ERR_USER_INPUT, "Mask must be of signed integer type.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(iter_max < self->axis_size[2], ERR_USER_INPUT, "Maximum number of iterations exceeds spectral axis size.");
		
	// Determine catalogue size
	const size_t cat_size = Catalog_get_size(cat);
	if(cat_size == 0)
	{
		warning("No sources in catalogue; skipping mask dilation.");
		return;
	}
	
	// Loop over all sources in catalogue
	for(size_t i = 0; i < cat_size; ++i)
	{
		// Extract source
		Source *src = Catalog_get_source(cat, i);
		message_verb(self->verbosity, "Source %zu", i + 1);
		
		// Get source ID & flag
		const long int src_id = Source_get_par_by_name_int(src, "id");
		ensure(src_id, ERR_USER_INPUT, "Source ID missing from catalogue; mask dilation failed.");
		long int flag = Source_get_par_by_name_int(src, "flag");
		
		// Get source bounding box
		const size_t x_min = Source_get_par_by_name_int(src, "x_min");
		const size_t x_max = Source_get_par_by_name_int(src, "x_max");
		const size_t y_min = Source_get_par_by_name_int(src, "y_min");
		const size_t y_max = Source_get_par_by_name_int(src, "y_max");
		size_t z_min = Source_get_par_by_name_int(src, "z_min");
		size_t z_max = Source_get_par_by_name_int(src, "z_max");
		ensure(x_min <= x_max && y_min <= y_max && z_min <= z_max, ERR_INDEX_RANGE, "Illegal source bounding box: min > max!");
		ensure(x_max < self->axis_size[0] && y_max < self->axis_size[1] && z_max < self->axis_size[2], ERR_INDEX_RANGE, "Source bounding box outside data cube boundaries.");
		
		// Get flux and check if source has negative flux
		double f_sum = Source_get_par_by_name_flt(src, "f_sum");
		double f_min = Source_get_par_by_name_flt(src, "f_min");
		double f_max = Source_get_par_by_name_flt(src, "f_max");
		size_t n_pix = Source_get_par_by_name_int(src, "n_pix");
		const bool is_negative = (f_sum < 0.0);
		
		// Additional variables
		double df_sum = 0.0;
		size_t z_min_new = z_min;
		size_t z_max_new = z_max;
		
		// Iterate
		for(size_t iter = 0; iter < iter_max; ++iter)
		{
			df_sum = 0.0;
			
			// Loop over source bounding box
			for(size_t z = z_min; z <= z_max; ++z)
			{
				for(size_t y = y_min; y <= y_max; ++y)
				{
					for(size_t x = x_min; x <= x_max; ++x)
					{
						const long int id = DataCube_get_data_int(mask, x, y, z);
						
						if(id == src_id)
						{
							// Check lower z
							if(z > 0)
							{
								const long int id_new = DataCube_get_data_int(mask, x, y, z - 1);
								if(id_new == 0)
								{
									if(IS_NOT_NAN(DataCube_get_data_flt(self, x, y, z - 1)))
									{
										DataCube_set_data_int(mask, x, y, z - 1, -1);
										df_sum += DataCube_get_data_flt(self, x, y, z - 1);
										if(z - 1 < z_min_new) z_min_new = z - 1;
									}
									else flag |= 4L;
								}
								else if(id_new > 0 && id_new != src_id) flag |= 8L;
							}
							else flag |= 2L;
							
							// Check higher z
							if(z < self->axis_size[2] - 1)
							{
								const long int id_new = DataCube_get_data_int(mask, x, y, z + 1);
								if(id_new == 0)
								{
									if(IS_NOT_NAN(DataCube_get_data_flt(self, x, y, z + 1)))
									{
										DataCube_set_data_int(mask, x, y, z + 1, -1);
										df_sum += DataCube_get_data_flt(self, x, y, z + 1);
										if(z + 1 > z_max_new) z_max_new = z + 1;
									}
									else flag |= 4L;
								}
								else if(id_new > 0 && id_new != src_id) flag |= 8L;
							}
							else flag |= 2L;
						}
					} // END loop over x
				} // END loop over y
			} // END loop over z
			
			// Check if flux increased within boundaries
			if(threshold < 0.0 || (is_negative && df_sum < threshold * f_sum) || (!is_negative && df_sum > threshold * f_sum))
			{
				// Mask should be grown
				f_sum += df_sum;
				z_min = z_min_new;
				z_max = z_max_new;
				
				// Loop over new bounding box
				for(size_t z = z_min; z <= z_max; ++z)
				{
					for(size_t y = y_min; y <= y_max; ++y)
					{
						for(size_t x = x_min; x <= x_max; ++x)
						{
							// If mask value is -1...
							if(DataCube_get_data_int(mask, x, y, z) == -1)
							{
								// ...switch to source ID...
								DataCube_set_data_int(mask, x, y, z, src_id);
								
								// ...and update n_pix, f_min and f_max if necessary
								const double value = DataCube_get_data_flt(self, x, y, z);
								if(value < f_min) f_min = value;
								if(value > f_max) f_max = value;
								++n_pix;
							}
						}
					}
				}
				
				message_verb(self->verbosity, " - Iteration %zu: df = %.3f (%.3f%%)", iter + 1, df_sum, 100.0 * df_sum / (f_sum - df_sum));
			}
			else
			{
				// No significant improvement; clean up again and exit
				for(size_t z = z_min_new; z <= z_max_new; ++z)
				{
					for(size_t y = y_min; y <= y_max; ++y)
					{
						for(size_t x = x_min; x <= x_max; ++x)
						{
							// Reset mask value again
							if(DataCube_get_data_int(mask, x, y, z) == -1) DataCube_set_data_int(mask, x, y, z, 0);
						}
					}
				}
				
				// Stop iterating
				break;
			}
		} // END iteration loop
		
		// Update source parameters
		Source_set_par_flt(src, "f_min", f_min, NULL, NULL);
		Source_set_par_flt(src, "f_max", f_max, NULL, NULL);
		Source_set_par_flt(src, "f_sum", f_sum, NULL, NULL);
		Source_set_par_int(src, "z_min", z_min, NULL, NULL);
		Source_set_par_int(src, "z_max", z_max, NULL, NULL);
		Source_set_par_int(src, "n_pix", n_pix, NULL, NULL);
		Source_set_par_int(src, "flag",  flag,  NULL, NULL);
		
		// Update progress bar
		progress_bar("Progress: ", i + 1, cat_size);
	}  // END source loop
	
	return;
}



/// @brief 2-D projection of 3-D mask
///
/// Public method for projecting a 3-D mask cube onto a 2-D mask
/// image along the spectral axis. It is implicitly assumed that
/// the mask cube is of integer type. The projected 2-D image will
/// be returned. Note that it is the caller's responsibility to run
/// the destructor on the returned image when it is no longer
/// required to release its memory again.
///
/// @param self  Object self-reference.
///
/// @return Projected 2-D mask image.

PUBLIC DataCube *DataCube_2d_mask(const DataCube *self)
{
	check_null(self);
	
	DataCube *maskImage = DataCube_blank(self->axis_size[0], self->axis_size[1], 1, self->data_type, self->verbosity);
	Header_copy_wcs(self->header, maskImage->header);
	Header_copy_misc(self->header, maskImage->header, true, true);
	
	// Loop over all pixels to project mask cube onto image
	#pragma omp parallel for schedule(static)
	for(size_t y = 0; y < self->axis_size[1]; ++y)
	{
		for(size_t x = self->axis_size[0]; x--;)
		{
			for(size_t z = self->axis_size[2]; z--;)
			{
				size_t value = DataCube_get_data_int(self, x, y, z);
				if(value)
				{
					DataCube_set_data_int(maskImage, x, y, 0, value);
					break;
				}
			}
		}
	}
	
	return maskImage;
}



/// @brief Flag regions in data cube
///
/// @param self    Object self-reference.
/// @param region  Array containing the regions to be flagged.
///                Must be of the form `x_min`, `x_max`, `y_min`,
///                `y_max`, `z_min`, `z_max`, ... where the
///                boundaries are inclusive.
///
/// Public method for flagging the specified regions in the data
/// cube. If the data cube is of floating-point type, then all
/// pixels to be flagged will be set to `NaN`. For integer cubes a
/// value of 0 will instead be used. The region must contain a
/// multiple of 6 entries of the form `x_min`, `x_max`, `y_min`, `y_max`,
/// `z_min`, `z_max`. Boundaries extending beyond the boundaries of the
/// cube will be automatically adjusted.

PUBLIC void DataCube_flag_regions(DataCube *self, const Array_siz *region)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(region);
	
	const size_t size = Array_siz_get_size(region);
	ensure(size % 6 == 0, ERR_USER_INPUT, "Flagging regions must contain a multiple of 6 entries.");
	
	message("Applying flags.");
	
	// Loop over regions
	for(size_t i = 0; i < size; i += 6)
	{
		// Establish boundaries
		size_t x_min = Array_siz_get(region, i + 0);
		size_t x_max = Array_siz_get(region, i + 1);
		size_t y_min = Array_siz_get(region, i + 2);
		size_t y_max = Array_siz_get(region, i + 3);
		size_t z_min = Array_siz_get(region, i + 4);
		size_t z_max = Array_siz_get(region, i + 5);
		
		// Adjust boundaries if necessary
		if(x_max >= self->axis_size[0]) x_max = self->axis_size[0] - 1;
		if(y_max >= self->axis_size[1]) y_max = self->axis_size[1] - 1;
		if(z_max >= self->axis_size[2]) z_max = self->axis_size[2] - 1;
		
		if(x_min > x_max) x_min = x_max;
		if(y_min > y_max) y_min = y_max;
		if(z_min > z_max) z_min = z_max;
		
		message_verb(self->verbosity, "  Region: [%zu, %zu, %zu, %zu, %zu, %zu]", x_min, x_max, y_min, y_max, z_min, z_max);
		
		for(size_t z = z_min; z <= z_max; ++z)
		{
			for(size_t y = y_min; y <= y_max; ++y)
			{
				for(size_t x = x_min; x <= x_max; ++x)
				{
					if(self->data_type < 0.0) DataCube_set_data_flt(self, x, y, z, NAN);
					else DataCube_set_data_int(self, x, y, z, 0);
				}
			}
		}
	}
	
	return;
}



/// @brief Flagging based on catalogue of positions
///
/// Public method for flagging positions specified in an external
/// source catalogue file. The file must contain just two columns
/// specifying the longitude and latitude of the positions to be
/// flagged. No other content (e.g. comments) is allowed. The
/// coordinates can either be specified in pixels (`coord_sys` = 1) or
/// in the native world coordinate system of the data cube (e.g. RA
/// and declination in units of degrees). Lastly, a radius of > 0
/// can be specified, in which case a circular region of that radius
/// (in pixels) will be flagged around the central coordinate.
///
/// @param self       Object self-reference.
/// @param filename   Name of the source catalogue file.
/// @param coord_sys  Pixel (0) or world (1) coordinates.
/// @param radius     Radius of flagging area in pixels.

PUBLIC void DataCube_continuum_flagging(DataCube *self, const char *filename, const int coord_sys, const long int radius)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type.");
	check_null(filename);
	
	// Set up a few variables
	const long int radius2 = radius * radius;
	const long int axis_size_x = (long int)(self->axis_size[0]);
	const long int axis_size_y = (long int)(self->axis_size[1]);
	const long int axis_size_z = (long int)(self->axis_size[2]);
	size_t counter = 0;
	
	// Read catalogue file
	Table *cont_cat = Table_from_file(filename, " \t,|");
	if(Table_rows(cont_cat) == 0 || Table_cols(cont_cat) != 2)
	{
		warning("Continuum catalogue non-compliant; must contain 2 data columns.\n         Flagging catalogue file will be ignored.");
		Table_delete(cont_cat);
		return;
	}
	
	// Set up WCS if requested
	WCS *wcs = NULL;
	if(coord_sys == 1)
	{
		wcs = DataCube_extract_wcs(self);
		if(wcs == NULL)
		{
			warning("WCS conversion failed; cannot apply flagging catalogue.");
			Table_delete(cont_cat);
			return;
		}
	}
	
	// Process catalogue line-by-line
	for(size_t i = 0; i < Table_rows(cont_cat); ++i)
	{
		// Extract longitude and latitude
		double lon = -1e+30;
		double lat = -1e+30;
		
		// Convert position from WCS to pixels if needed
		if(coord_sys == 1) WCS_convertToPixel(wcs, Table_get(cont_cat, i, 0), Table_get(cont_cat, i, 1), 0.0, &lon, &lat, NULL);
		
		// Ensure that source is within cube boundaries
		const long int pos_x = (long int)(lon + 0.5);
		const long int pos_y = (long int)(lat + 0.5);
		if(pos_x < 0 || pos_y < 0 || pos_x >= axis_size_x || pos_y >= axis_size_y) continue;
		++counter;
		
		// Establish bounding box
		const long int x_min = pos_x > radius ? pos_x - radius : 0;
		const long int y_min = pos_y > radius ? pos_y - radius : 0;
		const long int x_max = pos_x + radius < axis_size_x ? pos_x + radius : axis_size_x - 1;
		const long int y_max = pos_y + radius < axis_size_y ? pos_y + radius : axis_size_y - 1;
		
		for(long int z = 0; z < axis_size_z; ++z)
		{
			for(long int y = y_min; y <= y_max; ++y)
			{
				for(long int x = x_min; x <= x_max; ++x)
				{
					if((x - pos_x) * (x - pos_x) + (y - pos_y) * (y - pos_y) <= radius2) DataCube_set_data_flt(self, x, y, z, NAN);
				}
			}
		}
	}
	
	// Print some statistics
	message("Flagged %zu out of %zu positions from catalogue.", counter, Table_rows(cont_cat));
	
	// Clean up
	Table_delete(cont_cat);
	WCS_delete(wcs);
	
	return;
}



/// @brief Automatic flagging module
///
/// Public method for automatically determining spectral channels
/// and/or spatial pixels to be flagged. The algorithm first determines
/// the RMS in each channel or pixel. It then calculates the median of
/// the RMS values across all channels or pixels to determine the
/// typical RMS. Next, the median absolute deviation will be determined
/// as a measure of the scatter of the individual RMS values about the
/// median. Lastly, any pixels or channels with
///
/// \code{.c}
///   |rms - median| > threshold * MAD
/// \endcode
///
/// will be added to the array of regions to be flagged. The order of the
/// region specification is `x_min`, `x_max`, `y_min`, `y_max`, `z_min`,
/// `z_max`.
///
/// @param self       Object self-reference.
/// @param threshold  Threshold for flagging (see description below).
/// @param mode       Flagging mode; 0 = no flagging, 1 = channels,
///                   2 = pixels, 3 = channels + pixels.
/// @param region     Array containing the regions to be flagged.
///                   New regions to be flagged will be appended to
///                   any existing region specifications.

PUBLIC void DataCube_autoflag(const DataCube *self, const double threshold, const unsigned int mode, Array_siz *region)
{
	// Sanity checks
	check_null(self);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Automatic flagging will only work on floating-point data.");
	ensure(mode < 4, ERR_USER_INPUT, "Flagging mode must be 0 (false), 1 (channels), 2 (pixels) or 3 (true).");
	check_null(region);
	
	const char *mode_labels[] = {"disabled", "channels", "pixels", "channels + pixels"};
	const unsigned int id_spectral = 1;
	const unsigned int id_spatial  = 2;
	
	message("Running auto-flagger with the following settings:");
	message("  Mode:       %s", mode_labels[mode]);
	message("  Threshold:  %.1f * rms\n", threshold);
	
	// A few settings
	const size_t size_x  = self->axis_size[0];
	const size_t size_y  = self->axis_size[1];
	const size_t size_z  = self->axis_size[2];
	const size_t size_xy = size_x * size_y;
	
	// Auto-flag spectral channels if requested
	if(mode & id_spectral)
	{
		message("Auto-flagging of spectral channels:");
		size_t counter = 0;
		
		if(self->data_type == -32)
		{
			// 32-bit single-precision
			float *noise_array = (float *)memory(MALLOC, size_z, sizeof(float));
			const float *ptr = (float *)self->data;
			
			// Measure noise in each channel
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < size_z; ++i) noise_array[i] = robust_noise_2_flt(ptr + i * size_xy, size_xy);
			
			// Determine median of noise measurements
			const float median = median_safe_flt(noise_array, size_z, false);
			
			// Determine RMS via MAD
			const double rms = MAD_TO_STD * mad_val_flt(noise_array, size_z, median, 1, 0);
			
			// Check which channels exceed threshold
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < size_z; ++i)
			{
				if(fabs(noise_array[i] - median) > threshold * rms)
				{
					// Add channel to flagging regions
					#pragma omp critical
					{
						Array_siz_push(region, 0);
						Array_siz_push(region, size_x - 1);
						Array_siz_push(region, 0);
						Array_siz_push(region, size_y - 1);
						Array_siz_push(region, i);
						Array_siz_push(region, i);
						++counter;
					}
				}
			}
			
			// Delete noise array
			free(noise_array);
		}
		else
		{
			// 64-bit double-precision
			double *noise_array = (double *)memory(MALLOC, size_z, sizeof(double));
			const double *ptr = (double *)self->data;
			
			// Measure noise in each channel
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < size_z; ++i) noise_array[i] = robust_noise_2_dbl(ptr + i * size_xy, size_xy);
			
			// Determine median of noise measurements
			const double median = median_safe_dbl(noise_array, size_z, false);
			
			// Determine RMS via MAD
			const double rms = MAD_TO_STD * mad_val_dbl(noise_array, size_z, median, 1, 0);
			
			// Check which channels exceed threshold
			#pragma omp parallel for schedule(static)
			for(size_t i = 0; i < size_z; ++i)
			{
				if(fabs(noise_array[i] - median) > threshold * rms)
				{
					// Add channel to flagging regions
					#pragma omp critical
					{
						Array_siz_push(region, 0);
						Array_siz_push(region, size_x - 1);
						Array_siz_push(region, 0);
						Array_siz_push(region, size_y - 1);
						Array_siz_push(region, i);
						Array_siz_push(region, i);
						++counter;
					}
				}
			}
			
			// Delete noise array
			free(noise_array);
		}
		
		message("  %zu spectral channel%s marked for flagging.\n", counter, counter == 1 ? "" : "s");
	}
	
	// Auto-flag spatial pixels if requested
	if(mode & id_spatial)
	{
		message("Auto-flagging of spatial pixels:");
		size_t counter = 0;
		
		// 64-bit double-precision
		DataCube *noise_array = DataCube_blank(size_x, size_y, 1, -64, self->verbosity);
		
		// Loop over all pixels
		#pragma omp parallel
		{
			double *spectrum = (double *)memory(MALLOC, size_z, sizeof(double));
			
			#pragma omp for collapse(2) schedule(static)
			for(size_t y = 0; y < size_y; ++y)
			{
				for(size_t x = 0; x < size_x; ++x)
				{
					// Loop over all channels
					for(size_t z = size_z; z--;)
					{
						// Copy values into spectrum
						spectrum[z] = DataCube_get_data_flt(self, x, y, z);
					}
					
					// Measure noise across spectrum
					DataCube_set_data_flt(noise_array, x, y, 0, robust_noise_2_dbl(spectrum, size_z));
				}
			}
			
			// Delete spectrum again
			free(spectrum);
		}
		
		// Determine median of noise measurements
		const double median = median_safe_dbl((double *)(noise_array->data), size_xy, false);
		
		// Determine RMS via MAD
		const double rms = MAD_TO_STD * mad_val_dbl((double *)(noise_array->data), size_xy, median, 1, 0);
		
		// Check which pixels exceed threshold
		#pragma omp parallel for collapse(2) schedule(static)
		for(size_t y = 0; y < size_y; ++y)
		{
			for(size_t x = 0; x < size_x; ++x)
			{
				if(fabs(DataCube_get_data_flt(noise_array, x, y, 0) - median) > threshold * rms)
				{
					// Add pixel to flagging regions
					#pragma omp critical
					{
						Array_siz_push(region, x);
						Array_siz_push(region, x);
						Array_siz_push(region, y);
						Array_siz_push(region, y);
						Array_siz_push(region, 0);
						Array_siz_push(region, size_z - 1);
					}
					++counter;
				}
			}
		}
		
		// Delete noise array
		DataCube_delete(noise_array);
		
		message("  %zu spatial pixel%s marked for flagging.\n", counter, counter == 1 ? "" : "s");
	}
	
	return;
}



/// @brief Identify pixels with value of infinity for flagging
///
/// Public method for searching for values of infinity in
/// floating-point data cubes. Pixels with infinite values will
/// be added to the specified flagging region array for later
/// flagging.
///
/// @param self    Object self-reference.
/// @param region  Array for storage of flagging regions.
///
/// @return Number of pixels containing a value of infinity.

PUBLIC size_t DataCube_flag_infinity(const DataCube *self, Array_siz *region)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Flagging of infinity only possible for floating-point data.");
	check_null(region);
	
	message("Searching for values of infinity.");
	size_t counter = 0;
	
	// Loop over entire array
	#pragma omp parallel for schedule(static)
	for(size_t z = 0; z < self->axis_size[2]; ++z)
	{
		for(size_t y = self->axis_size[1]; y--;)
		{
			for(size_t x = self->axis_size[0]; x--;)
			{
				// Check for Inf
				if(isinf(DataCube_get_data_flt(self, x, y, z)))
				{
					// Add Inf pixel to flagging regions
					#pragma omp critical
					{
						Array_siz_push(region, x);
						Array_siz_push(region, x);
						Array_siz_push(region, y);
						Array_siz_push(region, y);
						Array_siz_push(region, z);
						Array_siz_push(region, z);
						++counter;
					}
				}
			}
		}
	}
	
	if(counter) message("  Found and flagged %zu infinite %s.", counter, counter == 1 ? "pixel" : "pixels");
	else message("  No infinite pixel values found.");
	
	return counter;
}



/// @brief Copy blanked pixels from one cube to another
///
/// Public method for copying blanked pixels from one cube to the
/// other. Both cubes need to be of the same size in x, y and z and
/// must be of floating-point type. Blanked pixels are assumed to
/// be represented by `NaN` (not a number).
///
/// @param self    Object self-reference.
/// @param source  Data cube from which to copy blanked pixels.

PUBLIC void DataCube_copy_blanked(DataCube *self, const DataCube *source)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(source);
	check_null(source->data);
	ensure((self->data_type == -32 || self->data_type == -64) && (source->data_type == -32 || source->data_type == -64), ERR_USER_INPUT, "Cannot copy blanked pixels; both data cubes must be floating-point.");
	ensure(self->axis_size[0] == source->axis_size[0] && self->axis_size[1] == source->axis_size[1] && self->axis_size[2] == source->axis_size[2], ERR_USER_INPUT, "Cannot copy blanked pixels; data cubes differ in size.");
	
	// Loop over entire array and copy blanks
	if(self->data_type == -32)
	{
		float *ptr_dst = (float *)(self->data) + self->data_size;
		
		if(source->data_type == -32)
		{
			for(float *ptr_src = (float *)(source->data) + source->data_size; ptr_src --> (float *)(source->data);)
			{
				--ptr_dst;
				if(IS_NAN(*ptr_src)) *ptr_dst = NAN;
			}
		}
		else
		{
			for(double *ptr_src = (double *)(source->data) + source->data_size; ptr_src --> (double *)(source->data);)
			{
				--ptr_dst;
				if(IS_NAN(*ptr_src)) *ptr_dst = NAN;
			}
		}
	}
	else
	{
		double *ptr_dst = (double *)(self->data) + self->data_size;
		
		if(source->data_type == -32)
		{
			for(float *ptr_src = (float *)(source->data) + source->data_size; ptr_src --> (float *)(source->data);)
			{
				--ptr_dst;
				if(IS_NAN(*ptr_src)) *ptr_dst = NAN;
			}
		}
		else
		{
			for(double *ptr_src = (double *)(source->data) + source->data_size; ptr_src --> (double *)(source->data);)
			{
				--ptr_dst;
				if(IS_NAN(*ptr_src)) *ptr_dst = NAN;
			}
		}
	}
	
	return;
}



/// @brief Return array index from x, y and z
///
/// Private method to turn a 3-D pixel coordinate into a 1-D array
/// index position under the fundamental assumption that the cube
/// is three-dimensional. Note that this function will still work
/// for 2-D arrays by simply setting `z` = 0 and the size of the
/// third axis to 1 (likewise for 1-D arrays).
///
/// @param self  Object self-reference.
/// @param x     First coordinate.
/// @param y     Second coordinate.
/// @param z     Third coordinate.
///
/// @return Returns the 1-D array index corresponding to the 3 coordinate
/// values specified under the assumption that the cube is 3-D.

PRIVATE inline size_t DataCube_get_index(const DataCube *self, const size_t x, const size_t y, const size_t z)
{
	return x + self->axis_size[0] * (y + self->axis_size[1] * z);
}



/// @brief Calculate x, y and z from array index
///
/// Private method for determining the array coordinates x, y and z
/// from the specified array index. The results will be written to
/// the specified pointers to `x`, `y` and `z`. Note that this will
/// also work for 2-D arrays for which the size of the third axis is 1
/// (and likewise for 1-D arrays); the resulting `z` (and/or y) will
/// be 0 in that case.
///
/// @param self   Object self-reference.
/// @param index  Index for which `x`, `y` and `z` are to be determined.
/// @param x      First coordinate.
/// @param y      Second coordinate.
/// @param z      Third coordinate.

PRIVATE void DataCube_get_xyz(const DataCube *self, const size_t index, size_t *x, size_t *y, size_t *z)
{
	*z = index / (self->axis_size[0] * self->axis_size[1]);
	const size_t ixy = index - self->axis_size[0] * self->axis_size[1] * *z;
	*y = ixy / self->axis_size[0];
	*x = ixy - self->axis_size[0] * *y;
	
	return;
}



/// @brief Run Smooth + Clip (S+C) finder on data cube
///
/// Public method for running the **Smooth + Clip** (S+C) finder on
/// the specified data cube. The S+C finder will smooth the data on
/// the specified spatial and spectral scales, applying a Gaussian
/// filter in the spatial domain and a boxcar filter in the spectral
/// domain. It will then measure the noise level in each iteration
/// and mark all pixels with absolute values greater than or equal
/// to the specified threshold (relative to the noise level) as 1
/// in the specified mask cube, which must be of 8-bit integer
/// type, while non-detected pixels will be set to a value of 0.
/// Pixels already detected in a previous iteration will be set to
/// `maskScaleXY` times the original rms noise level of the data
/// before smoothing. If the value of `maskScaleXY` is negative, no
/// replacement will be carried out.
///
/// The input data cube must be a 32 or 64-bit floating point data
/// array. The spatial kernel sizes must be positive floating point
/// values that represent the FWHM of the Gaussian kernels to be
/// applied in the spatial domain. The spectral kernel sizes must
/// be positive, odd integer numbers representing the widths of the
/// boxcar filters to be applied in the spectral domain. The threshold
/// is relative to the noise level and should be a floating
/// point number greater than about 3.0. Lastly, the value of
/// `maskScaleXY` times the original rms of the data will be used to
/// replace pixels in the data cube that were already detected in a
/// previous iteration. This is to ensure that any sources in the
/// data will not be smeared out beyond the extent of the source
/// when convolving with large kernel sizes. It will, however,
/// create a positive bias in the flux measurement of the source and
/// can therefore be disabled by setting `maskScaleXY` to a negative
/// value.
///
/// Several methods are available for **measuring the noise** in the
/// data cube, including the standard deviation, median absolute
/// deviation and a Gaussian fit to the flux histogram. These
/// differ in their speed and robustness. In addition, the flux range
/// used in the noise measurement can be restricted to negative or
/// positive pixels only to reduce the impact or actual emission or
/// absorption featured on the noise measurement.
///
/// It is also possible to **renormalise the noise level** after each
/// smoothing operation by setting `scaleNoise` to a value of 1 or 2.
/// It should be noted that local noise normalisation on smoothed data
/// is risky, as the number of statistically independent data samples
/// may have become too low for a reliable measurement of the noise.
///
/// @param self          Data cube to run the S+C finder on.
/// @param maskCube      Mask cube for recording detected pixels.
/// @param kernels_spat  List of spatial smoothing lengths corresponding
///                      to the FWHM of the Gaussian kernels to be
///                      applied; 0 = no smoothing.
/// @param kernels_spec  List of spectral smoothing lengths corresponding
///                      to the widths of the boxcar filters to be
///                      applied. Must be odd or 0.
/// @param threshold     Relative flux threshold to be applied.
/// @param maskScaleXY   Already detected pixels will be set to this
///                      value times the original rms of the data
///                      before smoothing the data again. If negative,
///                      no replacement will be carried out.
/// @param method        Method to use for measuring the noise in
///                      the smoothed copies of the cube; can be
///                      `NOISE_STAT_STD`, `NOISE_STAT_MAD` or
///                      `NOISE_STAT_GAUSS` for standard deviation,
///                      median absolute deviation and Gaussian fit
///                      to flux histogram, respectively.
/// @param range         Flux range to used in noise measurement, Can
///                      be -1, 0 or 1 for negative only, all or
///                      positive only.
/// @param scaleNoise    0 = no noise scaling; 1 = global noise scaling;
///                      2 = local noise scaling. Applied after each
///                      smoothing operation.
/// @param snStatistic   Statistic to use in the noise scaling. Can
///                      be `NOISE_STAT_STD` for standard deviation,
///                      `NOISE_STAT_MAD` for median absolute deviation
///                      or `NOISE_STAT_GAUSS` for Gaussian fitting to
///                      the flux histogram.
/// @param snRange       Flux range to be used in the noise scaling.
///                      Can be -1, 0 or +1 for negative range, full
///                      range or positive range, respectively.
/// @param snWindowXY    Spatial window size for local noise scaling.
///                      See DataCube_scale_noise_local() for details.
/// @param snWindowZ     Spectral window size for local noise scaling
///                      See DataCube_scale_noise_local() for details.
/// @param snGridXY      Spatial grid size for local noise scaling.
///                      See DataCube_scale_noise_local() for details.
/// @param snGridZ       Spectral grid size for local noise scaling.
///                      See DataCube_scale_noise_local() for details.
/// @param snInterpol    Enable interpolation for local noise scaling
///                      if true. See DataCube_scale_noise_local()
///                      for details.
/// @param start_time    Arbitrary time stamp; progress time of the
///                      algorithm will be calculated and printed
///                      relative to `start_time`.
/// @param start_clock   Arbitrary clock count; progress time of the
///                      algorithm in term of CPU time will be
///                      calculated and printed relative to `clock_time`.

PUBLIC void DataCube_run_scfind(const DataCube *self, DataCube *maskCube, const Array_dbl *kernels_spat, const Array_siz *kernels_spec, const double threshold, const double maskScaleXY, const noise_stat method, const int range, const int scaleNoise, const noise_stat snStatistic, const int snRange, const size_t snWindowXY, const size_t snWindowZ, const size_t snGridXY, const size_t snGridZ, const bool snInterpol, const time_t start_time, const clock_t start_clock)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	ensure(self->data_type < 0, ERR_USER_INPUT, "The S+C finder can only be applied to floating-point data.");
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(maskCube->data_type == 8, ERR_USER_INPUT, "Mask cube must be of 8-bit integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	check_null(kernels_spat);
	check_null(kernels_spec);
	ensure(Array_dbl_get_size(kernels_spat) && Array_siz_get_size(kernels_spec), ERR_USER_INPUT, "Invalid spatial or spectral kernel list encountered.");
	ensure(threshold >= 0.0, ERR_USER_INPUT, "Negative flux threshold encountered.");
	ensure(method == NOISE_STAT_STD || method == NOISE_STAT_MAD || method == NOISE_STAT_GAUSS, ERR_USER_INPUT, "Invalid noise measurement method: %d.", method);
	
	// A few additional settings
	const double FWHM_CONST = 2.0 * sqrt(2.0 * log(2.0));  // Conversion between sigma and FWHM of Gaussian function
	size_t cadence = self->data_size / NOISE_SAMPLE_SIZE;  // Stride for noise calculation
	if(cadence < 2) cadence = 1;
	else if(cadence % self->axis_size[0] == 0) cadence -= 1;    // Ensure stride is not equal to multiple of x-axis size
	message("Using a stride of %zu in noise measurement.\n", cadence);
	
	// Measure noise in original cube with sampling "cadence"
	double rms;
	double rms_smooth;
	
	if(method == NOISE_STAT_STD)      rms = DataCube_stat_std(self, 0.0, cadence, range);
	else if(method == NOISE_STAT_MAD) rms = MAD_TO_STD * DataCube_stat_mad(self, 0.0, cadence, range);
	else                              rms = DataCube_stat_gauss(self, cadence, range);
	
	// Run S+C finder for all smoothing kernels
	for(size_t i = 0; i < Array_dbl_get_size(kernels_spat); ++i)
	{
		for(size_t j = 0; j < Array_siz_get_size(kernels_spec); ++j)
		{
			message("Smoothing kernel:  [%.1f] x [%zu]", Array_dbl_get(kernels_spat, i), Array_siz_get(kernels_spec, j));
			
			// Check if any smoothing requested
			if(Array_dbl_get(kernels_spat, i) || Array_siz_get(kernels_spec, j))
			{
				// Smoothing required; create a copy of the original cube
				DataCube *smoothedCube = DataCube_copy(self);
				
				// Set flux of already detected pixels to maskScaleXY * rms
				if(maskScaleXY >= 0.0) DataCube_set_masked_8(smoothedCube, maskCube, maskScaleXY * rms);
				
				// Spatial and spectral smoothing
				if(Array_dbl_get(kernels_spat, i) > 0.0) DataCube_gaussian_filter(smoothedCube, Array_dbl_get(kernels_spat, i) / FWHM_CONST);
				if(Array_siz_get(kernels_spec, j) > 0)   DataCube_boxcar_filter(smoothedCube, Array_siz_get(kernels_spec, j) / 2);
				
				// Copy original blanks into smoothed cube again
				// (these were set to 0 during smoothing)
				DataCube_copy_blanked(smoothedCube, self);
				
				// Scale noise if requested
				if(scaleNoise == 1)
				{
					message("Correcting for noise variations along spectral axis.\n");
					DataCube_scale_noise_spec(smoothedCube, snStatistic, snRange);
				}
				else if(scaleNoise == 2)
				{
					message("Correcting for local noise variations.");
					DataCube *noiseCube = DataCube_scale_noise_local(
						smoothedCube,
						snStatistic,
						snRange,
						snWindowXY,
						snWindowZ,
						snGridXY,
						snGridZ,
						snInterpol
					);
					DataCube_delete(noiseCube);
				}
				
				// Calculate the RMS of the smoothed cube
				if(method == NOISE_STAT_STD)      rms_smooth = DataCube_stat_std(smoothedCube, 0.0, cadence, range);
				else if(method == NOISE_STAT_MAD) rms_smooth = MAD_TO_STD * DataCube_stat_mad(smoothedCube, 0.0, cadence, range);
				else                              rms_smooth = DataCube_stat_gauss(smoothedCube, cadence, range);
				
				message("Noise level:       %.3e", rms_smooth);
				
				// Add pixels above threshold to mask
				DataCube_mask_8(smoothedCube, maskCube, threshold * rms_smooth, 1);
				
				// Delete smoothed cube again
				DataCube_delete(smoothedCube);
			}
			else
			{
				// No smoothing required; apply threshold to original cube
				message("Noise level:       %.3e", rms);
				DataCube_mask_8(self, maskCube, threshold * rms, 1);
			}
			
			// Print time
			timestamp(start_time, start_clock);
		}
	}
	
	return;
}



/// @brief Run simple threshold finder on data cube
///
/// Public method for running a simple threshold finder on the data
/// cube specified by the user. Detected pixels will be added to
/// the mask cube provided, which must be of 8-bit integer type.
/// The specified flux threshold can either be absolute or relative
/// depending on the value of the `absolute` parameter. In the latter
/// case, the threshold will be multiplied by the noise level across
/// the data cube as measured using the method and range specified
/// by the user. In both cases, pixels with an absolute flux value
/// greater than the threshold will be added to the mask cube.
///
/// @param self       Data cube to run the threshold finder on.
/// @param maskCube   Mask cube for recording detected pixels.
/// @param absolute   If true, apply absolute threshold; otherwise
///                   multiply threshold by noise level.
/// @param threshold  Absolute or relative flux threshold.
/// @param method     Method to use for measuring the noise in
///                   the cube; can be `NOISE_STAT_STD`,
///                   `NOISE_STAT_MAD` or `NOISE_STAT_GAUSS` for
///                   standard deviation, median absolute deviation
///                   and Gaussian fit to flux histogram, respectively.
/// @param range      Flux range to used in noise measurement, Can
///                   be -1, 0 or 1 for negative only, all or
///                   positive only, respectively.

PUBLIC void DataCube_run_threshold(const DataCube *self, DataCube *maskCube, const bool absolute, double threshold, const noise_stat method, const int range)
{
	// Sanity checks
	check_null(self);
	ensure(self->data_type < 0, ERR_USER_INPUT, "The S+C finder can only be applied to floating-point data.");
	check_null(maskCube);
	ensure(maskCube->data_type == 8, ERR_USER_INPUT, "Mask cube must be of 8-bit integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(threshold >= 0.0, ERR_USER_INPUT, "Negative flux threshold encountered.");
	ensure(method == NOISE_STAT_STD || method == NOISE_STAT_MAD || method == NOISE_STAT_GAUSS, ERR_USER_INPUT, "Invalid noise measurement method: %d.", method);
	
	// Set threshold relative to noise level if requested
	if(!absolute)
	{
		// Set stride for noise calculation
		size_t cadence = self->data_size / NOISE_SAMPLE_SIZE;
		if(cadence < 2) cadence = 1;
		else if(cadence % self->axis_size[0] == 0) cadence -= 1;  // Ensure stride is not equal to multiple of x-axis size
		
		// Multiply threshold by rms
		double rms = 0.0;
		if(method == NOISE_STAT_STD)      rms = DataCube_stat_std(self, 0.0, cadence, range);
		else if(method == NOISE_STAT_MAD) rms = DataCube_stat_mad(self, 0.0, cadence, range) * MAD_TO_STD;
		else                              rms = DataCube_stat_gauss(self, cadence, range);
		message("- Noise level:      %.3e  (using stride of %zu)", rms, cadence);
		threshold *= rms;
	}
	
	// Apply threshold
	DataCube_mask_8(self, maskCube, threshold, 1);
	
	return;
}




/// @brief Link objects in an integer mask
///
/// Public method for linking objects recorded in an integer mask
/// within the specified merging radii. The mask must be a 32-bit
/// integer array with a background value of 0, while objects can
/// have any value `!= 0`. If values `!= 0` are present, they will be
/// set to `-1` at the start. The linker will first give objects that
/// are connected within the specified radii a unique label.
/// Objects that fall outside of the minimum or maximum size
/// requirements will be removed on the fly. If `pos_src` is set to
/// `true`, sources with negative total flux will also be removed.
/// Likewise, if `pos_pix` is `true`, only positive pixels will be
/// linked and negative ones discarded.
///
/// Several thresholds can be provided by the user. Objects that
/// fall outside these thresholds will be discarded by the linker.
/// Thresholds are specified with the `min_` or `max_` prefix and
/// include size in x, y and z, total number of pixels (`npix`) and
/// filling factor (`fill`), with the latter defined as the total
/// number of pixels divided by the number of pixels that make up
/// the rectangular bounding box. Any threshold can be set to 0 to
/// disable its application altogether.
///
/// @param self        Object self-reference.
/// @param mask        32-bit integer mask cube.
/// @param radius_x    Merging radius in x.
/// @param radius_y    Merging radius in y.
/// @param radius_z    Merging radius in z.
/// @param min_size_x  Minimum size requirement for objects in x.
/// @param min_size_y  Minimum size requirement for objects in y.
/// @param min_size_z  Minimum size requirement for objects in z.
/// @param min_npix    Minimum required number of pixels.
/// @param min_fill    Minimum required filling factor of object.
/// @param max_size_x  Maximum size requirement for objects in x.
/// @param max_size_y  Maximum size requirement for objects in y.
/// @param max_size_z  Maximum size requirement for objects in z.
/// @param max_npix    Maximum required number of pixels.
/// @param max_fill    Maximum required filling factor of object.
/// @param pos_pix     If `true`, negative pixels will be discarded.
/// @param pos_src     If `true`, negative sources will be discarded.
/// @param rms         Global rms value by which all flux values will
///                    be normalised. 1 = no normalisation.


PUBLIC LinkerPar *DataCube_run_linker(const DataCube *self, DataCube *mask, const size_t radius_x, const size_t radius_y, const size_t radius_z, const size_t min_size_x, const size_t min_size_y, const size_t min_size_z, const size_t min_npix, const double min_fill, const size_t max_size_x, const size_t max_size_y, const size_t max_size_z, const size_t max_npix, const double max_fill, const bool pos_pix, const bool pos_src, const double rms)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	ensure(mask->data_type == 32, ERR_USER_INPUT, "Linker will only accept 32-bit integer masks.");
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Data cube must be of floating-point type for linking.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	if(min_fill < 0.0) warning("Ignoring minimum filling factor; value is negative.");
	if(max_fill < 0.0) warning("Ignoring maximum filling factor; value is negative.");
	
	// Print some information
	message("Linker settings:");
	message(" - Merging radii:   %zu, %zu, %zu", radius_x, radius_y, radius_z);
	message(" - Minimum size:    %zu x %zu x %zu", min_size_x, min_size_y, min_size_z);
	if(max_size_x || max_size_y || max_size_z) message(" - Maximum size:    %zu x %zu x %zu", max_size_x, max_size_y, max_size_z);
	if(min_npix || max_npix) message(" - Min/max pixels:  %zu, %zu", min_npix, max_npix);
	if(min_fill > 0.0 || max_fill > 0.0) message(" - Min/max fill:    %.1f%%, %.1f%%", 100.0 * min_fill, 100.0 * max_fill);
	message(" - Keep negative:   %s\n", pos_src ? "no" : "yes");
	
	// Create empty linker parameter object
	LinkerPar *lpar = LinkerPar_new(self->verbosity);
	
	// Define a few parameters
	const size_t nx = mask->axis_size[0];
	const size_t ny = mask->axis_size[1];
	const size_t nz = mask->axis_size[2];
	const size_t max_x = mask->axis_size[0] ? mask->axis_size[0] - 1 : 0;
	const size_t max_y = mask->axis_size[1] ? mask->axis_size[1] - 1 : 0;
	const size_t max_z = mask->axis_size[2] ? mask->axis_size[2] - 1 : 0;
	const double rms_inv = 1.0 / rms;
	int32_t label = 1;
	int32_t *ptr_mask = (int32_t *)(mask->data);
	const size_t cadence = (nz / 100) ? nz / 100 : 1;  // Only used for updating progress bar
	
	// Link pixels into sources
	for(size_t z = nz; z--;)
	{
		if(z % cadence == 0) progress_bar("Progress: ", nz - z, nz);
		
		for(size_t y = ny; y--;)
		{
			for(size_t x = nx; x--;)
			{
				// Obtain index
				const size_t index = DataCube_get_index(mask, x, y, z);
				
				// Check if pixel is detected
				if(ptr_mask[index] < 0)
				{
					// Get flux value and check for NaN, Inf and positivity
					const double flux = DataCube_get_data_flt(self, x, y, z);
					if(!isfinite(flux) || (pos_pix && flux < 0.0))
					{
						ptr_mask[index] = 0;
						continue;
					}
					
					// Set pixel to label
					ptr_mask[index] = label;
					
					// Set quality flag
					unsigned char flag = 0;
					if(x == 0 || x == max_x || y == 0 || y == max_y) flag |= 1;
					if(z == 0 || z == max_z) flag |= 2;
					
					// Create a new linker parameter entry
					LinkerPar_push(lpar, label, x, y, z, flux * rms_inv, flag);
					
					// Recursively process neighbouring pixels
					Stack *stack = Stack_new();
					Stack_push(stack, index);
					DataCube_process_stack(self, mask, stack, radius_x, radius_y, radius_z, label, lpar, rms_inv, pos_pix);
					Stack_delete(stack);
					
					// Check if new source outside size (and other) requirements
					if(LinkerPar_get_obj_size(lpar, label, 0) < min_size_x
					|| LinkerPar_get_obj_size(lpar, label, 1) < min_size_y
					|| LinkerPar_get_obj_size(lpar, label, 2) < min_size_z
					|| (max_size_x && LinkerPar_get_obj_size(lpar, label, 0) > max_size_x)
					|| (max_size_y && LinkerPar_get_obj_size(lpar, label, 1) > max_size_y)
					|| (max_size_z && LinkerPar_get_obj_size(lpar, label, 2) > max_size_z)
					|| LinkerPar_get_npix(lpar, label) < min_npix
					|| (max_npix && LinkerPar_get_npix(lpar, label) > max_npix)
					|| (min_fill > 0.0 && (double)LinkerPar_get_npix(lpar, label) / (double)(LinkerPar_get_obj_size(lpar, label, 0) * LinkerPar_get_obj_size(lpar, label, 1) * LinkerPar_get_obj_size(lpar, label, 2)) < min_fill)
					|| (max_fill > 0.0 && (double)LinkerPar_get_npix(lpar, label) / (double)(LinkerPar_get_obj_size(lpar, label, 0) * LinkerPar_get_obj_size(lpar, label, 1) * LinkerPar_get_obj_size(lpar, label, 2)) > max_fill)
					|| (pos_src && LinkerPar_get_flux(lpar, label) < 0.0))
					{
						// Yes, it is -> discard source
						// Get source bounding box
						size_t x_min, x_max, y_min, y_max, z_min, z_max;
						LinkerPar_get_bbox(lpar, label, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
						
						// Set all source pixels to 0 in mask
						for(size_t zz = z_min; zz <= z_max; ++zz)
						{
							for(size_t yy = y_min; yy <= y_max; ++yy)
							{
								for(size_t xx = x_min; xx <= x_max; ++xx)
								{
									// Get index and mask value of pixel
									int32_t *ptr2 = ptr_mask + DataCube_get_index(mask, xx, yy, zz);
									if(*ptr2 == label) *ptr2 = 0;
								}
							}
						}
						
						// Discard source entry
						LinkerPar_pop(lpar);
					}
					else
					{
						// No it isn't -> retain source
						// Set flags as necessary
						size_t x_min, x_max, y_min, y_max, z_min, z_max;
						LinkerPar_get_bbox(lpar, label, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
						if(x_min == 0 || x_max == max_x || y_min == 0 || y_max == max_y) LinkerPar_update_flag(lpar, flag |= 1);
						if(z_min == 0 || z_max == max_z) LinkerPar_update_flag(lpar, flag |= 2);
						
						// Increment label
						ensure(++label > 0, ERR_INT_OVERFLOW, "Too many sources for 32-bit signed integer mask.");
					}
				}
			}
		}
	}
	
	// Print information
	LinkerPar_print_info(lpar);
	
	// Return LinkerPar object
	return lpar;
}



/// @brief Recursive method for labelling neighbouring pixels
///
/// Private method for checking whether any neighbouring pixels of
/// the object on the stack within the specified merging radii are
/// detected by the source finder (value of < 0). If so, their value
/// will be set to the same label as the object on the stack and the
/// LinkerPar object will be updated to include the new pixel. If
/// `pos_pix` is set to `true`, then only positive pixels will be
/// accepted and negative ones discarded.
///
/// The method will then process the neighbours of each neighbour
/// recursively by using an internal stack rather than recursive
/// function calls, which makes stack overflows controllable and
/// ensures that the stack is implemented on the heap to allow its
/// size to be dynamically adjusted and take up as much memory as
/// needed.
///
/// @param self      Object self-reference.
/// @param mask      32-bit mask cube.
/// @param stack     Stack object to be processed.
/// @param radius_x  Merging radius in x.
/// @param radius_y  Merging radius in y.
/// @param radius_z  Merging radius in z.
/// @param label     Label to be assigned to detected neighbours.
///                  Must be > 1, as 1 means not yet labelled.
/// @param lpar      Pointer to LinkerPar object containing the
///                  recorded object parameters. This will be updated
///                  whenever a new pixel is assigned to the same
///                  object currently getting linked.
/// @param rms_inv   Inverse of the global rms value by which all
///                  flux values will multiplied. If set to 1, no
///                  normalisation will occur.
/// @param pos_pix   If `true`, negative pixels will be discarded.

PRIVATE void DataCube_process_stack(const DataCube *self, DataCube *mask, Stack *stack, const size_t radius_x, const size_t radius_y, const size_t radius_z, const int32_t label, LinkerPar *lpar, const double rms_inv, const bool pos_pix)
{
	// Set up a few parameters
	size_t x, y, z;
	const size_t radius_x_squ   = radius_x * radius_x;
	const size_t radius_y_squ   = radius_y * radius_y;
	const size_t radius_z_squ   = radius_z * radius_z;
	const size_t radius_xy_squ  = radius_x_squ * radius_y_squ;
	const size_t radius_xz_squ  = radius_x_squ * radius_z_squ;
	const size_t radius_yz_squ  = radius_y_squ * radius_z_squ;
	const size_t radius_xyz_squ = radius_x_squ * radius_yz_squ;
	const float *data_flt = (float *)(self->data);
	const double *data_dbl = (double *)(self->data);
	int32_t *ptr_mask = (int32_t *)(mask->data);
	unsigned char flag = 0;
	
	// Loop until the stack is empty
	while(Stack_get_size(stack))
	{
		// Pop last element from stack and get its x, y and z coordinates
		DataCube_get_xyz(mask, Stack_pop(stack), &x, &y, &z);
		
		// Determine bounding box within which to search for neighbours
		const size_t x1 = (x > radius_x) ? (x - radius_x) : 0;
		const size_t y1 = (y > radius_y) ? (y - radius_y) : 0;
		const size_t z1 = (z > radius_z) ? (z - radius_z) : 0;
		const size_t x2 = (x + radius_x + 1 < mask->axis_size[0]) ? (x + radius_x) : (mask->axis_size[0] - 1);
		const size_t y2 = (y + radius_y + 1 < mask->axis_size[1]) ? (y + radius_y) : (mask->axis_size[1] - 1);
		const size_t z2 = (z + radius_z + 1 < mask->axis_size[2]) ? (z + radius_z) : (mask->axis_size[2] - 1);
		
		// Loop over entire bounding box
		for(size_t zz = z1; zz <= z2; ++zz)
		{
			const size_t dz_squ = zz > z ? (zz - z) * (zz - z) * radius_xy_squ : (z - zz) * (z - zz) * radius_xy_squ;
			
			for(size_t yy = y1; yy <= y2; ++yy)
			{
				const size_t dy_squ = yy > y ? (yy - y) * (yy - y) * radius_xz_squ : (y - yy) * (y - yy) * radius_xz_squ;
				
				for(size_t xx = x1; xx <= x2; ++xx)
				{
					const size_t dx_squ = xx > x ? (xx - x) * (xx - x) * radius_yz_squ : (x - xx) * (x - xx) * radius_yz_squ;
					
					// Check merging radius, assuming ellipsoid (with dx^2 / rx^2 + dy^2 / ry^2 + dz^2 / rz^2 = 1)
					if(dx_squ + dy_squ + dz_squ > radius_xyz_squ) continue;
					
					// Get index, mask value and flux of neighbour
					const size_t index = DataCube_get_index(mask, xx, yy, zz);
					int32_t *ptr = ptr_mask + index;
					
					// WARNING: The following implicitly assumes that data are of floating-point type!
					const double flux = (self->data_type == -32) ? data_flt[index] : data_dbl[index];
					
					// Check if NaN or Inf
					if(!isfinite(flux))
					{
						*ptr = 0;                                // unmask pixel
						LinkerPar_update_flag(lpar, flag |= 4);  // update flag
						continue;
					}
					
					// Check for positivity
					if(pos_pix && flux < 0.0)
					{
						*ptr = 0;  // unmask pixel
						continue;
					}
					
					// If detected, but not yet labelled
					if(*ptr < 0)
					{
						*ptr = label;                                              // label pixel
						LinkerPar_update(lpar, xx, yy, zz, flux * rms_inv, flag);  // update linker parameter object
						Stack_push(stack, index);                                  // push pixel onto stack
					}
				}
			}
		}
	}
	
	return;
}



/// @brief Source parameterisation
///
/// Public method for measuring advanced parameters of all sources
/// contained in the specified catalogue. The mask cube must be of
/// 32-bit integer type and must have the same dimensions as the
/// data cube. All sources found in the catalogue must also be
/// recorded in the mask with their catalogued source ID number.
/// All parameters derived by this method will be appended at the
/// end of the catalogue or updated if it already exists.
///
/// If `use_wcs` is set to `true`, the method will attempt to convert
/// certain parameters to WCS and append those to the catalogue in
/// addition to their pixel-based equivalents.
///
/// @param self      Object self-reference.
/// @param mask      32-bit mask cube.
/// @param cat       Catalogue of sources to be parameterised.
/// @param use_wcs   If `true`, attempt to convert the position of
///                  the source to WCS.
/// @param physical  If `true`, convert relevant parameters to
///                  physical units using information from the header.
///                  If `false`, native pixel units will be used.
/// @param prefix    Prefix to be used in source names. Defaults to
///                  `SoFiA` if set to `NULL`.

PUBLIC void DataCube_parameterise(const DataCube *self, const DataCube *mask, Catalog *cat, bool use_wcs, bool physical, const char *prefix)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	check_null(cat);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Parameterisation only possible with floating-point data.");
	ensure(mask->data_type > 0, ERR_USER_INPUT, "Mask must be of integer type.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	
	// Establish catalogue size
	const size_t cat_size = Catalog_get_size(cat);
	ensure(cat_size, ERR_USER_INPUT, "No sources in catalogue; nothing to parameterise.");
	message("Found %zu source%s in need of parameterisation.", cat_size, (cat_size > 1 ? "s" : ""));
	
	// Relevant header information
	String *unit_flux_dens = NULL;
	String *unit_flux = NULL;
	String *label_lon = NULL;
	String *label_lat = NULL;
	String *label_spec = NULL;
	String *ucd_lon = NULL;
	String *ucd_lat = NULL;
	String *ucd_spec = NULL;
	String *unit_lon = NULL;
	String *unit_lat = NULL;
	String *unit_spec = NULL;
	double beam_area = 0.0;
	double chan_size = 0.0;
	
	// Extract relevant header information
	DataCube_get_wcs_info(self, &unit_flux_dens, &unit_flux, &label_lon, &label_lat, &label_spec, &ucd_lon, &ucd_lat, &ucd_spec, &unit_lon, &unit_lat, &unit_spec, &beam_area, &chan_size);
	String *label_lon_peak = String_new(String_get(label_lon));
	String *label_lat_peak = String_new(String_get(label_lat));
	String *label_spec_peak = String_new(String_get(label_spec));
	String_append(label_lon_peak,  "_peak");
	String_append(label_lat_peak,  "_peak");
	String_append(label_spec_peak, "_peak");
	
	// Extract WCS object if requested
	WCS *wcs = NULL;
	use_wcs = use_wcs ? (wcs = DataCube_extract_wcs(self)) != NULL : use_wcs;
	
	// Establish if physical parameters can be calculated
	// (only supported if BUNIT is Jy/beam)
	if((physical = physical ? String_compare(unit_flux_dens, "Jy/beam") : physical)) message("Attempting to measure parameters in physical units.");
	
	// Create string holding source name
	String *source_name = String_new("");
	
	// Loop over all sources in catalogue
	for(size_t i = 0; i < cat_size; ++i)
	{
		// Extract source
		Source *src = Catalog_get_source(cat, i);
		
		// Extract source ID
		const size_t src_id = Source_get_par_by_name_int(src, "id");
		ensure(src_id, ERR_USER_INPUT, "Source ID missing from catalogue; cannot parameterise.");
		progress_bar("Progress: ", i + 1, cat_size);
		
		// Extract number of detected pixels
		const size_t n_pix = Source_get_par_by_name_int(src, "n_pix");
		
		// Extract source bounding box
		const size_t x_min = Source_get_par_by_name_int(src, "x_min");
		const size_t x_max = Source_get_par_by_name_int(src, "x_max");
		const size_t y_min = Source_get_par_by_name_int(src, "y_min");
		const size_t y_max = Source_get_par_by_name_int(src, "y_max");
		const size_t z_min = Source_get_par_by_name_int(src, "z_min");
		const size_t z_max = Source_get_par_by_name_int(src, "z_max");
		ensure(x_min <= x_max && y_min <= y_max && z_min <= z_max, ERR_INDEX_RANGE, "Illegal source bounding box: min > max!");
		ensure(x_max < self->axis_size[0] && y_max < self->axis_size[1] && z_max < self->axis_size[2], ERR_INDEX_RANGE, "Source bounding box outside data cube boundaries.");
		
		const size_t nx = x_max - x_min + 1;
		const size_t ny = y_max - y_min + 1;
		const size_t nz = z_max - z_min + 1;
		
		// Check if source has negative flux
		const bool is_negative = (Source_get_par_by_name_flt(src, "f_sum") < 0.0);
		
		// Initialise source parameters
		double rms = 0.0;
		double pos_x = 0.0;
		double pos_y = 0.0;
		double pos_z = 0.0;
		size_t pos_x_peak = 0.0;
		size_t pos_y_peak = 0.0;
		size_t pos_z_peak = 0.0;
		double f_sum = 0.0;
		double f_min = INFINITY;
		double f_max = -INFINITY;
		double w50 = 0.0;
		double w20 = 0.0;
		double wm50 = 0.0;
		double err_x = 0.0;
		double err_y = 0.0;
		double err_z = 0.0;
		double err_f_sum = 0.0;
		double longitude = 0.0;
		double latitude = 0.0;
		double spectral = 0.0;
		double longitude_peak = 0.0;
		double latitude_peak = 0.0;
		double spectral_peak = 0.0;
		double ell_maj = 0.0;
		double ell_min = 0.0;
		double ell_pa = 0.0;
		double ell3s_maj = 0.0;
		double ell3s_min = 0.0;
		double ell3s_pa = 0.0;
		double kin_pa = 0.0;
		
		// Auxiliary parameters and storage
		double *kpa_cenX = (double *)memory(MALLOC, nz, sizeof(double));
		double *kpa_cenY = (double *)memory(MALLOC, nz, sizeof(double));
		double *kpa_sum  = (double *)memory(MALLOC, nz, sizeof(double));
		size_t kpa_first = z_max - z_min;
		size_t kpa_last  = 0;
		size_t kpa_counter = 0;
		
		Array_dbl *array_rms = Array_dbl_new(0);
		double *spectrum   = (double *)memory(CALLOC, nz, sizeof(double));
		double *moment_map = (double *)memory(CALLOC, nx * ny,   sizeof(double));
		size_t *count_map  = (size_t *)memory(CALLOC, nx * ny,   sizeof(size_t));
		
		double sum_pos = 0.0;
		
		// First pass
		for(size_t z = z_min; z <= z_max; ++z)
		{
			for(size_t y = y_min; y <= y_max; ++y)
			{
				for(size_t x = x_min; x <= x_max; ++x)
				{
					const size_t id    = DataCube_get_data_int(mask, x, y, z);
					const double value = is_negative ? -DataCube_get_data_flt(self, x, y, z) : DataCube_get_data_flt(self, x, y, z);
					
					if(id == src_id)
					{
						// ALL PIXELS
						// Flux
						f_sum += value;
						if(f_min > value) f_min = value;
						if(f_max < value) f_max = value;
						
						// Moment map for ellipse fitting and peak position
						moment_map[x - x_min + nx * (y - y_min)] += value;
						count_map [x - x_min + nx * (y - y_min)] += 1;
						
						// Spectrum for line width
						spectrum[z - z_min] += value;
						
						// POSITIVE PIXELS ONLY
						if(value > 0.0)
						{
							// Centroid position
							pos_x += value * x;
							pos_y += value * y;
							pos_z += value * z;
							sum_pos += value;
						}
					}
					else if(id == 0)
					{
						// Retain non-source pixels for RMS measurement
						Array_dbl_push(array_rms, value);
					}
				}
			}
		}
		
		// Finalise centroid
		pos_x /= sum_pos;
		pos_y /= sum_pos;
		pos_z /= sum_pos;
		
		// Measure position of peak in moment map and spectrum
		double value_peak = 0.0;
		for(size_t y = 0; y < ny; ++y)
		{
			for(size_t x = 0; x < nx; ++x)
			{
				if(moment_map[x + nx * y] > value_peak)
				{
					value_peak = moment_map[x + nx * y];
					pos_x_peak = x + x_min;
					pos_y_peak = y + y_min;
				}
			}
		}
		
		value_peak = 0.0;
		for(size_t z = 0; z < nz; ++z)
		{
			if(spectrum[z] > value_peak)
			{
				value_peak = spectrum[z];
				pos_z_peak = z + z_min;
			}
		}
		
		// Measure local RMS
		if(Array_dbl_get_size(array_rms)) rms = MAD_TO_STD * mad_val_dbl(Array_dbl_get_ptr(array_rms), Array_dbl_get_size(array_rms), 0.0, 1, 0);
		else warning_verb(self->verbosity, "Failed to measure local noise level for source %zu.", src_id);
		
		// Second pass
		// (for parameters dependent on first-pass output)
		for(size_t z = z_min; z <= z_max; ++z)
		{
			kpa_cenX[z - z_min] = 0.0;
			kpa_cenY[z - z_min] = 0.0;
			kpa_sum[z - z_min] = 0.0;
			
			for(size_t y = y_min; y <= y_max; ++y)
			{
				for(size_t x = x_min; x <= x_max; ++x)
				{
					const size_t id    = DataCube_get_data_int(mask, x, y, z);
					const double value = is_negative ? -DataCube_get_data_flt(self, x, y, z) : DataCube_get_data_flt(self, x, y, z);
					
					if(id == src_id)
					{
						// POSITIVE PIXELS
						if(value > 0.0)
						{
							err_x += ((double)x - pos_x) * ((double)x - pos_x);
							err_y += ((double)y - pos_y) * ((double)y - pos_y);
							err_z += ((double)z - pos_z) * ((double)z - pos_z);
						}
						
						// PIXELS > 3 SIGMA
						if(value > 3.0 * rms)
						{
							// Centroids for kinematic major axis
							kpa_cenX[z - z_min] += value * (double)x;
							kpa_cenY[z - z_min] += value * (double)y;
							kpa_sum [z - z_min] += value;
						}
					}
				}
			}
			
			// Determine centroid in each channel
			if(kpa_sum[z - z_min] > 0.0)
			{
				kpa_cenX[z - z_min] /= kpa_sum[z - z_min];
				kpa_cenY[z - z_min] /= kpa_sum[z - z_min];
				++kpa_counter;
				
				if(kpa_first > z - z_min) kpa_first = z - z_min;
				if(kpa_last  < z - z_min) kpa_last  = z - z_min;
			}
		}
		
		// Measure kinematic major axis
		if(kpa_counter < 2)
		{
			warning_verb(self->verbosity, "Failed to determine kinematic major axis for source %zu.\n         Emission is too faint.", i);
			kin_pa = -1.0;
		}
		else
		{
			if(kpa_counter == 2) warning_verb(self->verbosity, "Kinematic major axis for source %zu based on just 2 data points.", i);
			kin_pa = kin_maj_axis_dbl(kpa_cenX, kpa_cenY, kpa_sum, nz, kpa_first, kpa_last);
		}
		
		// Ellipse fit to moment-0 map
		moment_ellipse_fit_dbl(moment_map, count_map, nx, ny, pos_x - x_min, pos_y - y_min, rms, &ell_maj, &ell_min, &ell_pa, &ell3s_maj, &ell3s_min, &ell3s_pa);
		
		// Determine w20 and w50 from spectrum (moving inwards)
		spectral_line_width_dbl(spectrum, nz, &w20, &w50);
		
		// Determine wm50
		wm50 = wm50_line_width_dbl(spectrum, nz);
		
		// Determine uncertainties
		err_x = sqrt(err_x) * rms / sum_pos;
		err_y = sqrt(err_y) * rms / sum_pos;
		err_z = sqrt(err_z) * rms / sum_pos;
		err_f_sum = rms * sqrt(n_pix);
		
		// Carry out WCS conversion if requested
		if(use_wcs)
		{
			WCS_convertToWorld(wcs, pos_x, pos_y, pos_z, &longitude, &latitude, &spectral);
			WCS_convertToWorld(wcs, (double)pos_x_peak, (double)pos_y_peak, (double)pos_z_peak, &longitude_peak, &latitude_peak, &spectral_peak);
			DataCube_create_src_name(self, &source_name, prefix, longitude, latitude, label_lon);
		}
		else
		{
			String_set(source_name, prefix ? prefix : "SoFiA");
			String_append_int(source_name, "-%04zu", src_id);
		}
		
		// Invert flux-related parameters of negative sources
		if(is_negative)
		{
			swap(&f_min, &f_max);
			f_min = -f_min;
			f_max = -f_max;
			f_sum = -f_sum;
		}
		
		// Update catalogue entries
		Source_set_identifier(src, String_get(source_name));
		Source_set_par_flt(src, "x",      pos_x,      "pix",                      "pos.cartesian.x");
		Source_set_par_flt(src, "y",      pos_y,      "pix",                      "pos.cartesian.y");
		Source_set_par_flt(src, "z",      pos_z,      "pix",                      "pos.cartesian.z");
		Source_set_par_flt(src, "rms",    rms,        String_get(unit_flux_dens), "instr.det.noise");
		Source_set_par_flt(src, "f_min",  f_min,      String_get(unit_flux_dens), "phot.flux.density;stat.min");
		Source_set_par_flt(src, "f_max",  f_max,      String_get(unit_flux_dens), "phot.flux.density;stat.max");
		
		if(physical)
		{
			Source_set_par_flt(src, "f_sum", f_sum * chan_size / beam_area, String_get(unit_flux), "phot.flux");
			Source_set_par_flt(src, "w20",   w20 * chan_size,               String_get(unit_spec), "spect.line.width");
			Source_set_par_flt(src, "w50",   w50 * chan_size,               String_get(unit_spec), "spect.line.width");
			Source_set_par_flt(src, "wm50",  wm50 * chan_size,              String_get(unit_spec), "spect.line.width");
		}
		else
		{
			Source_set_par_flt(src, "f_sum", f_sum, String_get(unit_flux_dens), "phot.flux");
			Source_set_par_flt(src, "w20",   w20,   "pix",                      "spect.line.width");
			Source_set_par_flt(src, "w50",   w50,   "pix",                      "spect.line.width");
			Source_set_par_flt(src, "wm50",  wm50,  "pix",                      "spect.line.width");
		}
		
		Source_set_par_flt(src, "ell_maj",   ell_maj,   "pix", "phys.angSize");
		Source_set_par_flt(src, "ell_min",   ell_min,   "pix", "phys.angSize");
		Source_set_par_flt(src, "ell_pa",    ell_pa,    "deg", "pos.posAng");
		Source_set_par_flt(src, "ell3s_maj", ell3s_maj, "pix", "phys.angSize");
		Source_set_par_flt(src, "ell3s_min", ell3s_min, "pix", "phys.angSize");
		Source_set_par_flt(src, "ell3s_pa",  ell3s_pa,  "deg", "pos.posAng");
		Source_set_par_flt(src, "kin_pa",    kin_pa,    "deg", "pos.posAng");
		
		if(physical)
		{
			Source_set_par_flt(src, "err_x",     err_x * sqrt(beam_area),                 "pix",                 "stat.error;pos.cartesian.x");
			Source_set_par_flt(src, "err_y",     err_y * sqrt(beam_area),                 "pix",                 "stat.error;pos.cartesian.y");
			Source_set_par_flt(src, "err_z",     err_z * sqrt(beam_area),                 "pix",                 "stat.error;pos.cartesian.z");
			Source_set_par_flt(src, "err_f_sum", err_f_sum * chan_size / sqrt(beam_area), String_get(unit_flux), "stat.error;phot.flux");
		}
		else
		{
			Source_set_par_flt(src, "err_x",     err_x,     "pix",                      "stat.error;pos.cartesian.x");
			Source_set_par_flt(src, "err_y",     err_y,     "pix",                      "stat.error;pos.cartesian.y");
			Source_set_par_flt(src, "err_z",     err_z,     "pix",                      "stat.error;pos.cartesian.z");
			Source_set_par_flt(src, "err_f_sum", err_f_sum, String_get(unit_flux_dens), "stat.error;phot.flux");
		}
		
		if(use_wcs)
		{
			// Centroid
			Source_set_par_flt(src, String_get(label_lon),  longitude, String_get(unit_lon),  String_get(ucd_lon));
			Source_set_par_flt(src, String_get(label_lat),  latitude,  String_get(unit_lat),  String_get(ucd_lat));
			Source_set_par_flt(src, String_get(label_spec), spectral,  String_get(unit_spec), String_get(ucd_spec));
		}
		
		Source_set_par_int(src, "x_peak", pos_x_peak, "pix", "pos.cartesian.x");
		Source_set_par_int(src, "y_peak", pos_y_peak, "pix", "pos.cartesian.y");
		Source_set_par_int(src, "z_peak", pos_z_peak, "pix", "pos.cartesian.z");
		
		if(use_wcs)
		{
			// Peak
			Source_set_par_flt(src, String_get(label_lon_peak),  longitude_peak, String_get(unit_lon),  String_get(ucd_lon));
			Source_set_par_flt(src, String_get(label_lat_peak),  latitude_peak,  String_get(unit_lat),  String_get(ucd_lat));
			Source_set_par_flt(src, String_get(label_spec_peak), spectral_peak,  String_get(unit_spec), String_get(ucd_spec));
		}
		
		// Clean up (per source)
		Array_dbl_delete(array_rms);
		free(spectrum);
		free(moment_map);
		free(count_map);
		
		free(kpa_cenX);
		free(kpa_cenY);
		free(kpa_sum);
	}
	
	// Clean up (globally)
	WCS_delete(wcs);
	String_delete(unit_flux_dens);
	String_delete(unit_flux);
	String_delete(unit_lon);
	String_delete(unit_lat);
	String_delete(unit_spec);
	String_delete(label_lon);
	String_delete(label_lat);
	String_delete(label_spec);
	String_delete(label_lon_peak);
	String_delete(label_lat_peak);
	String_delete(label_spec_peak);
	String_delete(ucd_lon);
	String_delete(ucd_lat);
	String_delete(ucd_spec);
	String_delete(source_name);
	
	return;
}



/// @brief Extract WCS-related keywords from header
///
/// Private method for extracting header information related to the
/// World Coordinate System (WCS) of the data cube. All information
/// will be written to the variable pointers specified by the user.
///
/// @param self            Object self-reference.
/// @param unit_flux_dens  String to hold flux density unit.
/// @param unit_flux       String to hold flux unit.
/// @param label_lon       String to hold longitude axis name.
/// @param label_lat       String to hold latitude axis name.
/// @param label_spec      String to hold spectral axis name.
/// @param ucd_lon         String to hold UCD for longitude axis.
/// @param ucd_lat         String to hold UCD for latitude axis.
/// @param ucd_spec        String to hold UCD for spectral axis.
/// @param unit_lon        String to hold longitude axis unit.
/// @param unit_lat        String to hold latitude axis unit.
/// @param unit_spec       String to hold spectral axis unit.
/// @param beam_area       Variable to hold beam solid angle in pixels.
/// @param chan_size       Variable to hold spectral channel width.

PRIVATE void DataCube_get_wcs_info(const DataCube *self, String **unit_flux_dens, String **unit_flux, String **label_lon, String **label_lat, String **label_spec, String **ucd_lon, String **ucd_lat, String **ucd_spec, String **unit_lon, String **unit_lat, String **unit_spec, double *beam_area, double *chan_size)
{
	// Extract flux density unit from header
	*unit_flux_dens = Header_get_string(self->header, "BUNIT");
	if(String_size(*unit_flux_dens) == 0)
	{
		warning_verb(self->verbosity, "No flux unit (\'BUNIT\') defined in header.");
		String_set(*unit_flux_dens, "???");
	}
	else String_trim(*unit_flux_dens);
	
	// Fix commonly encountered misspellings
	if(String_compare(*unit_flux_dens, "JY/BEAM") || String_compare(*unit_flux_dens, "Jy/Beam")) String_set(*unit_flux_dens, "Jy/beam");
	
	// Make flux unit the same as flux density unit (might get updated later)
	*unit_flux = String_copy(*unit_flux_dens);
	
	// Determine axis types, units and UCDs
	*label_lon  = String_new("lon");
	*label_lat  = String_new("lat");
	*label_spec = String_new("spec");
	*ucd_lon    = String_new("");
	*ucd_lat    = String_new("");
	*ucd_spec   = String_new("");
	*unit_lon   = String_trim(Header_get_string(self->header, "CUNIT1"));
	*unit_lat   = String_trim(Header_get_string(self->header, "CUNIT2"));
	*unit_spec  = String_trim(Header_get_string(self->header, "CUNIT3"));
	if(String_size(*unit_lon) == 0) String_set(*unit_lon, "deg");
	if(String_size(*unit_lat) == 0) String_set(*unit_lat, "deg");
	
	// Longitude axis
	if(DataCube_cmphd(self, "CTYPE1", "RA--", 4))
	{
		String_set(*label_lon, "ra");
		String_set(*ucd_lon, "pos.eq.ra");
	}
	else if(DataCube_cmphd(self, "CTYPE1", "GLON", 4))
	{
		String_set(*label_lon, "l");
		String_set(*ucd_lon, "pos.galactic.lon");
	}
	else warning("Unsupported CTYPE1 value. Supported: RA, GLON.");
	
	// Latitude axis
	if(DataCube_cmphd(self, "CTYPE2", "DEC-", 4))
	{
		String_set(*label_lat, "dec");
		String_set(*ucd_lat, "pos.eq.dec");
	}
	else if(DataCube_cmphd(self, "CTYPE2", "GLAT", 4))
	{
		String_set(*label_lat, "b");
		String_set(*ucd_lat, "pos.galactic.lat");
	}
	else warning("Unsupported CTYPE2 value. Supported: DEC, GLAT.");
	
	// Spectral axis
	if(DataCube_cmphd(self, "CTYPE3", "FREQ", 4))
	{
		String_set(*label_spec, "freq");
		String_set(*ucd_spec, "em.freq");
		if(String_size(*unit_spec) == 0) String_set(*unit_spec, "Hz");
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VRAD", 4))
	{
		String_set(*label_spec, "v_rad");
		String_set(*ucd_spec, "spect.dopplerVeloc.radio");
		if(String_size(*unit_spec) == 0) String_set(*unit_spec, "m/s");
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VOPT", 4) || DataCube_cmphd(self, "CTYPE3", "FELO", 4))
	{
		String_set(*label_spec, "v_opt");
		String_set(*ucd_spec, "spect.dopplerVeloc.opt");
		if(String_size(*unit_spec) == 0) String_set(*unit_spec, "m/s");
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VELO", 4))
	{
		String_set(*label_spec, "v_app");
		String_set(*ucd_spec, "spect.dopplerVeloc");
		if(String_size(*unit_spec) == 0) String_set(*unit_spec, "m/s");
	}
	else
	{
		warning("Unsupported CTYPE3 value. Supported: FREQ, VRAD, VOPT, VELO.");
		if(String_size(*unit_spec) == 0) String_set(*unit_spec, "???");
	}
	
	// Extract spectral channel width
	*chan_size = fabs(Header_get_flt(self->header, "CDELT3"));
	
	if(IS_NAN(*chan_size))
	{
		warning("Header keyword \'CDELT3\' not found; assuming value of 1.");
		*chan_size = 1.0;
	}
	
	// Extract beam solid angle in pixels
	*beam_area = DataCube_get_beam_area(self);
	
	if(IS_NAN(*beam_area))
	{
		*beam_area = 1.0;
		String_append(*unit_flux, "*");
		String_append(*unit_flux, String_get(*unit_spec));
	}
	else
	{
		String_set(*unit_flux, "Jy*");  // WARNING: Flux density unit 'Jy' currently hard-coded here!!!
		String_append(*unit_flux, String_get(*unit_spec));
	}
	
	return;
}



/// @brief Generate source name from WCS information
///
/// Private method for generating a source name based on the
/// coordinates and WCS information specified by the user. The name
/// will consist of a prefix (`SoFiA` by default) followed by a
/// space followed by the coordinate part of the source, the format
/// of which will depend on the prefix and coordinate type. If the
/// prefix is `WALLABY`, then the official WALLABY source naming
/// convention will be used:
///
///  * `WALLABY Jhhmmss-ddmmss`
///
/// In all other cases the source name will be
///
///  * `prefix (J/B)hhmmss.ss-ddmmss.s`  for equatorial coordinates,
///  * `prefix Glll.llll-dd.dddd`        for Galactic coordinates, and
///  * `prefix lll.llll-dd.dddd`         otherwise.
///
/// The final source name will be written to the String pointer
/// specified by the user.
///
/// @param self         Object self-reference.
/// @param source_name  String to hold source name.
/// @param prefix       C string specifying the desired prefix.
/// @param longitude    Longitude of the source.
/// @param latitude     Latitude of the source.
/// @param label_lon    Longitude axis name.

PRIVATE void DataCube_create_src_name(const DataCube *self, String **source_name, const char *prefix, const double longitude, const double latitude, const String *label_lon)
{
	// Create source name
	String_set(*source_name, prefix ? prefix : "SoFiA");
	String_append(*source_name, " ");
	
	if(String_compare(label_lon, "ra"))
	{
		// Equatorial coordinates; try to figure out equinox
		double equinox = Header_get_flt(self->header, "EQUINOX");
		if(IS_NAN(equinox)) equinox = Header_get_flt(self->header, "EPOCH");
		
		if(equinox < 2000.0) String_append(*source_name, "B");  // assume Besselian equinox
		else String_append(*source_name, "J");                  // assume Julian equinox as the default
		
		// Determine coordinate part
		const double ra  = longitude / 15.0;  // WARNING: Assuming degrees!
		const double rah = floor(ra);
		const double ram = floor(60.0 * (ra - rah));
		const double ras = 3600.0 * (ra - rah - ram / 60.0);
		
		String_append_int(*source_name, "%02d", (int)rah);
		String_append_int(*source_name, "%02d", (int)ram);
		if(strcmp(prefix, "WALLABY")) String_append_flt(*source_name, "%05.2f", ras);  // round to 2 decimals in general
		else String_append_int(*source_name, "%02d", (int)ras);                        // round down if 'WALLABY'
		
		const double de  = fabs(latitude);    // WARNING: Assuming degrees!
		const double ded = floor(de);
		const double dem = floor(60.0 * (de - ded));
		const double des = 3600.0 * (de - ded - dem / 60.0);
		
		String_append(*source_name, latitude < 0.0 ? "-" : "+");
		String_append_int(*source_name, "%02d", (int)ded);
		String_append_int(*source_name, "%02d", (int)dem);
		if(strcmp(prefix, "WALLABY")) String_append_flt(*source_name, "%04.1f", des);  // round to 1 decimal in general
		else String_append_int(*source_name, "%02d", (int)des);                        // round down if 'WALLABY'
	}
	else if(String_compare(label_lon, "glon"))
	{
		// Galactic coordinates
		String_append(*source_name, "G");
		String_append_flt(*source_name, "%08.4f", longitude);
		String_append(*source_name, longitude < 0.0 ? "-" : "+");
		String_append_flt(*source_name, "%07.4f", latitude);
	}
	else
	{
		// Unknown coordinates
		String_append_flt(*source_name, "%08.4f", longitude);
		String_append(*source_name, longitude < 0.0 ? "-" : "+");
		String_append_flt(*source_name, "%07.4f", latitude);
	}
	
	return;
}



/// @brief Generate moment maps from data cube
///
/// Public method for generating spectral moment maps from the
/// specified data cube for all pixels that are `!= 0` in the mask.
/// The generated maps will be pointed to by the `mom0`, `mom1` and
/// `mom2` pointers provided in the function call. Note that these must
/// be uninitialised pointers to a DataCube object, i.e. they must **not**
/// be pointing to any valid DataCube object before being passed on
/// to the function. It is the user's responsibility to call the
/// destructor on each of the moment maps once they are no longer
/// required.
///
/// If positive is set to true, then only pixels with positive flux
/// will contribute to the calculation of the first and second
/// moment maps. This can be useful to prevent large negative signals
/// from affecting the moment calculation.
///
/// @param  self       Object self-reference.
/// @param  mask       32-bit mask cube.
/// @param  mom0       Pointer to a data cube object that will be
///                    pointing to the generated moment 0 map.
/// @param  mom1       Pointer to a data cube object that will be
///                    pointing to the generated moment 1 map.
/// @param  mom2       Pointer to a data cube object that will be
///                    pointing to the generated moment 2 map.
/// @param  chan       Pointer to a data cube object that will be
///                    pointing to the generated map containing the
///                    number of channels per pixel in the moment 0
///                    map.
/// @param  snr        Pointer to a data cube object that will be
///                    pointing to the generated map containing the
///                    signal-to-noise ratio (SNR) per pixel in the
///                    moment 0 map.
/// @param  obj_name   Name of the object for `OBJECT` header entry.
///                    If `NULL`, no `OBJECT` entry will be created.
/// @param  use_wcs    If `true`, convert channel numbers to WCS.
/// @param  threshold  Flux threshold to be used in the calculation
///                    of moment 1 and 2. Note that moment 0 will
///                    always include all channels to avoid a
///                    positive flux bias.
/// @param  rms        If positive, then convert the number of
///                    channels map to a proper SNR map. Set to 0
///                    to disable this conversion.

PUBLIC void DataCube_create_moments(const DataCube *self, const DataCube *mask, DataCube **mom0, DataCube **mom1, DataCube **mom2, DataCube **chan, DataCube **snr, const char *obj_name, bool use_wcs, const double threshold, const double rms)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Moment maps only possible with floating-point data.");
	ensure(mask->data_type > 0, ERR_USER_INPUT, "Mask must be of integer type.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	
	// Is data cube a 2-D image?
	const bool is_3d = DataCube_get_axis_size(self, 2) > 1;
	if(!is_3d) warning("Image is not 3D; moments 1 and 2 will not be created.");
	
	// Extract WCS information if requested
	WCS *wcs = NULL;
	use_wcs = (use_wcs && is_3d) ? (wcs = DataCube_extract_wcs(self)) != NULL : false;  // Ensure that WCS information is valid and cube is 3D
	
	// Extract spectral unit
	String *unit_spec  = Header_get_string(self->header, "CUNIT3");
	String_trim(unit_spec);
	
	if(String_size(unit_spec) == 0 && is_3d)
	{
		if(DataCube_cmphd(self, "CTYPE3", "FREQ", 4)) String_set(unit_spec,  "Hz");
		else if(DataCube_cmphd(self, "CTYPE3", "VRAD", 4) || DataCube_cmphd(self, "CTYPE3", "VOPT", 4) || DataCube_cmphd(self, "CTYPE3", "VELO", 4) || DataCube_cmphd(self, "CTYPE3", "FELO", 4)) String_set(unit_spec,  "m/s");
		else warning("Unsupported CTYPE3 value. Supported: FREQ, VRAD, VOPT, VELO.");
	}
	
	// Extract flux unit
	String *unit_flux_dens = Header_get_string(self->header, "BUNIT");
	String_trim(unit_flux_dens);
	
	// Fix commonly encountered misspellings
	if(String_compare(unit_flux_dens, "JY/BEAM") || String_compare(unit_flux_dens, "Jy/Beam")) String_set(unit_flux_dens, "Jy/beam");
	
	// Multiply flux unit by spectral unit
	if(use_wcs)
	{
		if(String_size(unit_flux_dens) == 0) warning_verb(self->verbosity, "No flux unit (\'BUNIT\') defined in header.");
		else String_append(unit_flux_dens, "*");
		String_append(unit_flux_dens, String_get(unit_spec));
	}
	
	// Initialise everything with NULL
	*mom0 = NULL;
	*mom1 = NULL;
	*mom2 = NULL;
	*chan = NULL;
	*snr  = NULL;
	
	// Create empty moment 0 map
	*mom0 = DataCube_blank(self->axis_size[0], self->axis_size[1], 1, -32, self->verbosity);
	
	// Create additional map containing summed flux for moment 1 and 2 calculations
	DataCube *sum_pos = NULL;
	
	// Copy WCS and other header elements from data cube to moment map
	Header_copy_wcs(self->header, (*mom0)->header);
	Header_copy_misc(self->header, (*mom0)->header, true, true);
	if(use_wcs) Header_set_str((*mom0)->header, "BUNIT", String_get(unit_flux_dens));
	if(obj_name != NULL) Header_set_str((*mom0)->header, "OBJECT", obj_name);
	
	if(is_3d)
	{
		// 3-D cube; create empty moment 1 and 2 maps (by copying empty moment 0 map)
		*mom1 = DataCube_copy(*mom0);
		*mom2 = DataCube_copy(*mom0);
		
		// Create empty sum map (by copying empty moment 0 map)
		sum_pos = DataCube_copy(*mom0);
		
		// Create empty channel map of 32-bit integer type
		*chan = DataCube_blank(self->axis_size[0], self->axis_size[1], 1, 32, self->verbosity);
		Header_copy_wcs(self->header, (*chan)->header);
		Header_copy_misc(self->header, (*chan)->header, false, true);
		if(obj_name != NULL) Header_set_str((*chan)->header, "OBJECT", obj_name);
		
		// Set BUNIT keyword in moments 1 and 2 and channel map
		Header_set_str((*mom1)->header, "BUNIT", use_wcs ? String_get(unit_spec) : " ");
		Header_set_str((*mom2)->header, "BUNIT", use_wcs ? String_get(unit_spec) : " ");
		Header_set_str((*chan)->header, "BUNIT", " ");
	}
	
	// Determine moments 0 and 1
	// NOTE: This should not be multi-threaded, as otherwise the moment-0 map would
	//       no longer be deterministic due to the arbitrary order of the summation
	//       which would result in slightly different rounding errors. While those
	//       differences are negligible, the moment maps from different runs would
	//       no longer be binary-identical, making unit testing impossible.
	for(size_t z = self->axis_size[2]; z--;)
	{
		double spectral = z;
		if(use_wcs) WCS_convertToWorld(wcs, 0, 0, z, NULL, NULL, &spectral);
		
		for(size_t y = self->axis_size[1]; y--;)
		{
			for(size_t x = self->axis_size[0]; x--;)
			{
				if(DataCube_get_data_int(mask, x, y, z))
				{
					const double flux = DataCube_get_data_flt(self, x, y, z);
					DataCube_add_data_flt(*mom0, x, y, 0, flux);
					
					if(is_3d)
					{
						DataCube_add_data_int(*chan, x, y, 0, 1);
						
						if(flux > threshold)
						{
							DataCube_add_data_flt(*mom1, x, y, 0, flux * spectral);
							DataCube_add_data_flt(sum_pos, x, y, 0, flux);
						}
					}
				}
			}
		}
	}
	
	// If image is 2-D then return, as only mom0 needed and nothing else left to do
	if(!is_3d) return;
	
	// Convert channel map to SNR map if requested
	if(rms > 0.0)
	{
		// Create empty SNR map
		*snr = DataCube_blank(self->axis_size[0], self->axis_size[1], 1, -32, self->verbosity);
		Header_copy_wcs(self->header, (*snr)->header);
		Header_copy_misc(self->header, (*snr)->header, false, true);
		Header_set_str((*snr)->header, "BUNIT", " ");
		if(obj_name != NULL) Header_set_str((*snr)->header, "OBJECT", obj_name);
		
		// Calculate SNR values
		#pragma omp parallel for collapse(2) schedule(static)
		for(size_t y = 0; y < self->axis_size[1]; ++y)
		{
			for(size_t x = 0; x < self->axis_size[0]; ++x)
			{
				const long int value_chan = DataCube_get_data_int(*chan, x, y, 0);
				const double   value_flux = DataCube_get_data_flt(*mom0, x, y, 0);
				if(value_chan > 0) DataCube_set_data_flt(*snr, x, y, 0, value_flux / (rms * sqrt((double)value_chan)));
				else DataCube_set_data_flt(*snr, x, y, 0, NAN);
			}
		}
	}
	
	// Continue with creation of moments 1 and 2
	// Divide moment 1 by moment 0
	#pragma omp parallel for collapse(2) schedule(static)
	for(size_t y = 0; y < self->axis_size[1]; ++y)
	{
		for(size_t x = 0; x < self->axis_size[0]; ++x)
		{
			const double flux = DataCube_get_data_flt(sum_pos, x, y, 0);
			if(flux > 0.0) DataCube_set_data_flt(*mom1, x, y, 0, DataCube_get_data_flt(*mom1, x, y, 0) / flux);
			else DataCube_set_data_flt(*mom1, x, y, 0, NAN);
		}
	}
	
	// Determine moment 2
	for(size_t z = self->axis_size[2]; z--;)
	{
		double spectral = z;
		if(use_wcs) WCS_convertToWorld(wcs, 0, 0, z, NULL, NULL, &spectral);
		
		for(size_t y = self->axis_size[1]; y--;)
		{
			for(size_t x = self->axis_size[0]; x--;)
			{
				if(DataCube_get_data_int(mask, x, y, z))
				{
					const double flux = DataCube_get_data_flt(self, x, y, z);
					
					if(flux > threshold)
					{
						const double velo = DataCube_get_data_flt(*mom1, x, y, 0) - spectral;
						DataCube_add_data_flt(*mom2, x, y, 0, velo * velo * flux);
					}
				}
			}
		}
	}
	
	// Divide moment 2 by summed flux density and take square root.
	#pragma omp parallel for collapse(2) schedule(static)
	for(size_t y = 0; y < self->axis_size[1]; ++y)
	{
		for(size_t x = 0; x < self->axis_size[0]; ++x)
		{
			// Moment 2
			const double flux = DataCube_get_data_flt(sum_pos, x, y, 0);
			const double sigma = DataCube_get_data_flt(*mom2, x, y, 0);
			if(flux > 0.0 && sigma > 0.0) DataCube_set_data_flt(*mom2, x, y, 0, sqrt(sigma / flux));
			else DataCube_set_data_flt(*mom2, x, y, 0, NAN);
		}
	}
	
	// Multiply moment 0 by CDELT3 if requested.
	if(use_wcs) DataCube_multiply_const(*mom0, fabs(Header_get_flt(self->header, "CDELT3")));
	
	// Clean up
	DataCube_delete(sum_pos);
	WCS_delete(wcs);
	String_delete(unit_flux_dens);
	String_delete(unit_spec);
	
	return;
}



/// @brief Create position-velocity (PV) diagram
///
/// Public method for generating a 2D position-velocity (PV) map from the
/// specified data cube through the spatial pixel position (`x0`, `y0`) with
/// a position angle specified by `angle`. The step size along the line
/// defined by the centre and position angle is specified by `step_size`
/// (in pixels). Bi-linear interpolation is used to interpolate the flux
/// values along the specified line. A DataCube object containing the PV
/// diagram will be returned.
///
/// @param self       Object self-reference.
/// @param x0         X coordinate of the centre position of the PV diagram
///                   in pixels.
/// @param y0         Y coordinate of the centre position of the PV diagram
///                   in pixels.
/// @param angle      Position angle (in radians) of the line through the
///                   centre along which the PV diagram is extracted. The
///                   angle increases anti-clockwise, with 0 pointing up
///                   (i.e. in increasing direction of the 2nd axis).
/// @param step_size  Spatial step size (in pixels) of the PV diagram. This
///                   does not need to be an integer number, but can be a
///                   fraction of the original pixel size (e.g. 0.5 for
///                   half-pixel increments).
/// @param obj_name   Name of the object to be written into the `OBJECT`
///                   header card of the output FITS file. If set to `NULL`
///                   then no `OBJECT` card will be written.
///
/// @return  Pointer to `DataCube` object containing the PV diagram.

PUBLIC DataCube *DataCube_create_pv(const DataCube *self, const double x0, const double y0, const double angle, const double step_size, const char *obj_name)
{
	// Sanity checks
	check_null(self);
	
	// Get current date and time
	char current_time_string[32];
	time_t current_time = time(NULL);
	strftime(current_time_string, 32, "%Y-%m-%dT%H:%M:%S", gmtime(&current_time));
	
	// Sort out geometry
	const size_t nx = DataCube_get_axis_size(self, 0);
	const size_t ny = DataCube_get_axis_size(self, 1);
	const size_t nz = DataCube_get_axis_size(self, 2);
	const size_t steps = (double)(nx > ny ? nx : ny) / (2.0 * step_size);
	
	// Create empty PV diagram
	DataCube *pv = DataCube_blank(2 * steps + 1, nz, 1, -32, self->verbosity);
	
	// Create/copy WCS header entries
	char value[FITS_HEADER_VALUE_SIZE + 1];
	DataCube_puthd_str(pv, "CTYPE1", "ANGLE   ");
	if(DataCube_chkhd(self, "CUNIT2"))
	{
		DataCube_gethd_str(self, "CUNIT2", value);
		DataCube_puthd_str(pv, "CUNIT1", value);
	}
	if(DataCube_chkhd(self, "CDELT2")) DataCube_puthd_flt(pv, "CDELT1", DataCube_gethd_flt(self, "CDELT2") * step_size);
	DataCube_puthd_flt(pv, "CRVAL1", 0.0);
	DataCube_puthd_flt(pv, "CRPIX1", steps + 1);  // NOTE: FITS is 1-based, hence +1 here!
	if(DataCube_chkhd(self, "CTYPE3"))
	{
		DataCube_gethd_str(self, "CTYPE3", value);
		DataCube_puthd_str(pv, "CTYPE2", value);
	}
	if(DataCube_chkhd(self, "CUNIT3"))
	{
		DataCube_gethd_str(self, "CUNIT3", value);
		DataCube_puthd_str(pv, "CUNIT2", value);
	}
	if(DataCube_chkhd(self, "CDELT3")) DataCube_puthd_flt(pv, "CDELT2", DataCube_gethd_flt(self, "CDELT3"));
	if(DataCube_chkhd(self, "CRVAL3")) DataCube_puthd_flt(pv, "CRVAL2", DataCube_gethd_flt(self, "CRVAL3"));
	if(DataCube_chkhd(self, "CRPIX3")) DataCube_puthd_flt(pv, "CRPIX2", DataCube_gethd_flt(self, "CRPIX3"));
	DataCube_puthd_str(pv, "ORIGIN", SOFIA_VERSION_FULL " (" SOFIA_CREATION_DATE ")");
	DataCube_puthd_str(pv, "DATE", current_time_string);
	if(DataCube_chkhd(self, "BUNIT"))
	{
		DataCube_gethd_str(self, "BUNIT", value);
		DataCube_puthd_str(pv, "BUNIT", value);
	}
	if(obj_name != NULL) DataCube_puthd_str(pv, "OBJECT", obj_name);
	if(DataCube_chkhd(self, "RESTFREQ")) DataCube_puthd_flt(pv, "RESTFREQ", DataCube_gethd_flt(self, "RESTFREQ"));
	if(DataCube_chkhd(self, "RESTFRQ"))  DataCube_puthd_flt(pv, "RESTFRQ",  DataCube_gethd_flt(self, "RESTFRQ"));
	if(DataCube_chkhd(self, "VELREF"))  DataCube_puthd_int(pv, "VELREF",  DataCube_gethd_int(self, "VELREF"));
	if(DataCube_chkhd(self, "SPECSYS"))
	{
		DataCube_gethd_str(self, "SPECSYS", value);
		DataCube_puthd_str(pv, "SPECSYS", value);
	}
	
	// Extract PV diagram
	for(size_t x = 0; x <= 2 * steps; ++x)
	{
		// Work out position
		const double dr = step_size * ((double)x - (double)steps);
		const double x_new = x0 - dr * sin(angle);
		const double y_new = y0 + dr * cos(angle);
		
		// Determine 4 nearest pixels
		const size_t x1 = (size_t)floor(x_new);
		const size_t x2 = (size_t)ceil(x_new);
		const size_t y1 = (size_t)floor(y_new);
		const size_t y2 = (size_t)ceil(y_new);
		
		// Blank if pixels are beyond axis range
		if(x1 >= nx || x2 >= nx || y1 >= ny || y2 >= ny || x2 <= x1 || y2 <= y1)
		{
			for(size_t z = 0; z < nz; ++z) DataCube_set_data_flt(pv, x, z, 0, NAN);
			continue;
		}
		
		for(size_t z = 0; z < nz; ++z)
		{
			// Bi-linear interpolation
			const double f1 = (x2 - x_new) * DataCube_get_data_flt(self, x1, y1, z) + (x_new - x1) * DataCube_get_data_flt(self, x2, y1, z);
			const double f2 = (x2 - x_new) * DataCube_get_data_flt(self, x1, y2, z) + (x_new - x1) * DataCube_get_data_flt(self, x2, y2, z);
			DataCube_set_data_flt(pv, x, z, 0, (y2 - y_new) * f1 + (y_new - y1) * f2);
		}
	}
	
	return pv;
}



/// @brief Create cubelets and other source products
///
/// Public method for generating cubelets and other data products
/// for each source in the specified mask and catalogue. The method
/// will generate cut-outs of the data cube and mask cube around
/// each source and also generate moment maps (0-2) and integrate
/// spectra. All data products will be saved to disc and then
/// deleted again.
///
/// @param  self       Object self-reference (data cube).
/// @param  mask       Mask cube.
/// @param  cat        Source catalogue.
/// @param  basename   Base name to be used for output files.
/// @param  overwrite  Replace existing files (`true`) or not (`false`)?
/// @param  use_wcs    Try to convert channel numbers to WCS?
/// @param  physical   If `true`, correct flux for beam solid angle.
/// @param  margin     Margin in pixels to be added around each
///                    source. If 0, sources will be cut out exactly.
/// @param  threshold  Flux threshold to be used for moment 1 and 2.
/// @param  offset_z   Offset to be added to channel numbers in spectrum.
/// @param  par        SoFiA parameter settings; these will be added
///                    to the output FITS file history in the header.
///                    If `NULL`, then no history will be written.

PUBLIC void DataCube_create_cubelets(const DataCube *self, const DataCube *mask, const Catalog *cat, const char *basename, const bool overwrite, bool use_wcs, bool physical, const size_t margin, const double threshold, const size_t offset_z, const Parameter *par)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(mask);
	check_null(mask->data);
	check_null(cat);
	ensure(self->data_type == -32 || self->data_type == -64, ERR_USER_INPUT, "Cubelets only possible with floating-point data.");
	ensure(mask->data_type > 0, ERR_USER_INPUT, "Mask must be of integer type.");
	ensure(self->axis_size[0] == mask->axis_size[0] && self->axis_size[1] == mask->axis_size[1] && self->axis_size[2] == mask->axis_size[2], ERR_USER_INPUT, "Data cube and mask cube have different sizes.");
	ensure(Catalog_get_size(cat), ERR_USER_INPUT, "Empty source catalogue provided.");
	
	// Create string for file names
	String *filename_template = String_append(String_new(basename), "_");
	String *filename = String_new("");
	
	// Extract flux unit from header
	String *unit_flux_dens = Header_get_string(self->header, "BUNIT");
	if(String_size(unit_flux_dens) == 0)
	{
		warning_verb(self->verbosity, "No flux unit (\'BUNIT\') defined in header.");
		String_set(unit_flux_dens, "???");
	}
	else String_trim(unit_flux_dens);
	
	// Fix commonly encountered misspellings
	if(String_compare(unit_flux_dens, "JY/BEAM") || String_compare(unit_flux_dens, "Jy/Beam")) String_set(unit_flux_dens, "Jy/beam");
	
	// Make flux unit the same as flux density for now (may get updated later)
	String *unit_flux = String_copy(unit_flux_dens);
	
	// Check if physical parameters can be derived
	// (only supported if BUNIT is Jy/beam)
	physical = physical ? String_compare(unit_flux_dens, "Jy/beam") : physical;
	
	// Extract WCS information if requested
	WCS *wcs = NULL;
	use_wcs = use_wcs ? (wcs = DataCube_extract_wcs(self)) != NULL : false;
	
	// Extract spectral unit from header
	String *label_spec = String_trim(Header_get_string(self->header, "CTYPE3"));
	String *unit_spec  = String_trim(Header_get_string(self->header, "CUNIT3"));
	
	if(DataCube_cmphd(self, "CTYPE3", "FREQ", 4))
	{
		String_set(label_spec, "freq");
		if(String_size(unit_spec) == 0) String_set(unit_spec,  "Hz");  // FITS default
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VRAD", 4))
	{
		String_set(label_spec, "v_rad");
		if(String_size(unit_spec) == 0) String_set(unit_spec,  "m/s");  // FITS default
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VOPT", 4) || DataCube_cmphd(self, "CTYPE3", "FELO", 4))
	{
		String_set(label_spec, "v_opt");
		if(String_size(unit_spec) == 0) String_set(unit_spec,  "m/s");  // FITS default
	}
	else if(DataCube_cmphd(self, "CTYPE3", "VELO", 4))
	{
		String_set(label_spec, "v_app");
		if(String_size(unit_spec) == 0) String_set(unit_spec,  "m/s");  // FITS default
	}
	else
	{
		warning("Unsupported CTYPE3 value. Supported: FREQ, VRAD, VOPT, VELO.");
		if(String_size(unit_spec) == 0) String_set(unit_spec, "???");  // Unit undefined
	}
	
	// Extract beam solid angle in pixels from header
	double beam_area = 1.0;
	
	if(physical)
	{
		// Get beam area
		beam_area = DataCube_get_beam_area(self);
		if(IS_NAN(beam_area)) beam_area = 1.0;
		else String_set(unit_flux, "Jy");
	}
	
	// Loop over all sources in the catalogue
	for(size_t i = 0; i < Catalog_get_size(cat); ++i)
	{
		const Source *src = Catalog_get_source(cat, i);
		
		// Get source ID
		const size_t src_id = Source_get_par_by_name_int(src, "id");
		ensure(src_id, ERR_USER_INPUT, "Source ID missing from catalogue; cannot create cubelets.");
		
		// Get local RMS
		double rms = Source_get_par_by_name_flt(src, "rms");
		if IS_NAN(rms) rms = 0.0;
		
		// Get source bounding box
		size_t x_min = Source_get_par_by_name_int(src, "x_min");
		size_t x_max = Source_get_par_by_name_int(src, "x_max");
		size_t y_min = Source_get_par_by_name_int(src, "y_min");
		size_t y_max = Source_get_par_by_name_int(src, "y_max");
		size_t z_min = Source_get_par_by_name_int(src, "z_min");
		size_t z_max = Source_get_par_by_name_int(src, "z_max");
		ensure(x_min <= x_max && y_min <= y_max && z_min <= z_max, ERR_INDEX_RANGE, "Illegal source bounding box: min > max!");
		ensure(x_max < self->axis_size[0] && y_max < self->axis_size[1] && z_max < self->axis_size[2], ERR_INDEX_RANGE, "Source bounding box outside data cube boundaries.");
		
		// Add margin if requested
		if(margin)
		{
			x_min = margin > x_min ? 0 : x_min - margin;
			y_min = margin > y_min ? 0 : y_min - margin;
			z_min = margin > z_min ? 0 : z_min - margin;
			x_max = x_max + margin < self->axis_size[0] ? x_max + margin : self->axis_size[0] - 1;
			y_max = y_max + margin < self->axis_size[1] ? y_max + margin : self->axis_size[1] - 1;
			z_max = z_max + margin < self->axis_size[2] ? z_max + margin : self->axis_size[2] - 1;
		}
		
		const size_t nx = x_max - x_min + 1;
		const size_t ny = y_max - y_min + 1;
		const size_t nz = z_max - z_min + 1;
		
		// Create empty cubelet
		DataCube *cubelet = DataCube_blank(nx, ny, nz, self->data_type, self->verbosity);
		
		// Copy and adjust header information
		Header_copy_wcs(self->header, cubelet->header);
		Header_adjust_wcs_to_subregion(cubelet->header, x_min, x_max, y_min, y_max, z_min, z_max);
		Header_copy_misc(self->header, cubelet->header, true, true);
		Header_set_str(cubelet->header, "OBJECT", Source_get_identifier(src));
		
		// Create empty masklet
		DataCube *masklet = DataCube_blank(nx, ny, nz, 8, self->verbosity);
		
		// Copy and adjust header information
		Header_copy_wcs(self->header, masklet->header);
		Header_adjust_wcs_to_subregion(masklet->header, x_min, x_max, y_min, y_max, z_min, z_max);
		Header_set_str(masklet->header, "BUNIT", " ");
		Header_set_str(masklet->header, "OBJECT", Source_get_identifier(src));
		
		// Create data array for spectrum
		double *spectrum = (double *)memory(CALLOC, nz, sizeof(double));
		size_t *pixcount = (size_t *)memory(CALLOC, nz, sizeof(size_t));
		
		// Copy data into cubelet, masklet and spectrum
		for(size_t z = z_min; z <= z_max; ++z)
		{
			for(size_t y = y_min; y <= y_max; ++y)
			{
				for(size_t x = x_min; x <= x_max; ++x)
				{
					// Cubelet
					DataCube_set_data_flt(cubelet, x - x_min, y - y_min, z - z_min, DataCube_get_data_flt(self, x, y, z));
					
					// Masklet
					const size_t id = DataCube_get_data_int(mask, x, y, z);
					if(id == src_id)
					{
						DataCube_set_data_int(masklet, x - x_min, y - y_min, z - z_min, 1);
						spectrum[z - z_min] += DataCube_get_data_flt(self, x, y, z);
						pixcount[z - z_min] += 1;
					}
					else DataCube_set_data_int(masklet, x - x_min, y - y_min, z - z_min, 0);
				}
			}
		}
		
		// Create moment maps etc.
		DataCube *mom0;
		DataCube *mom1;
		DataCube *mom2;
		DataCube *chan;
		DataCube *snr;
		DataCube_create_moments(cubelet, masklet, &mom0, &mom1, &mom2, &chan, &snr, Source_get_identifier(src), use_wcs, threshold * rms, rms);
		
		// Create PV diagram
		DataCube *pv = DataCube_create_pv(cubelet, Source_get_par_by_name_flt(src, "x") - x_min, Source_get_par_by_name_flt(src, "y") - y_min, Source_get_par_by_name_flt(src, "kin_pa") * M_PI / 180.0, 1.0, Source_get_identifier(src));
		
		// Save output products...
		// ...cubelet
		String_set(filename, String_get(filename_template));
		String_append_int(filename, "%ld", src_id);
		String_append(filename, "_cube.fits");
		DataCube_add_history(cubelet, par);
		DataCube_save(cubelet, String_get(filename), overwrite, DESTROY);
		
		// ...masklet
		String_set(filename, String_get(filename_template));
		String_append_int(filename, "%ld", src_id);
		String_append(filename, "_mask.fits");
		DataCube_add_history(masklet, par);
		DataCube_save(masklet, String_get(filename), overwrite, DESTROY);
		
		// ...moment maps
		if(mom0 != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_mom0.fits");
			DataCube_add_history(mom0, par);
			DataCube_save(mom0, String_get(filename), overwrite, DESTROY);
		}
		
		if(mom1 != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_mom1.fits");
			DataCube_add_history(mom1, par);
			DataCube_save(mom1, String_get(filename), overwrite, DESTROY);
		}
		
		if(mom2 != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_mom2.fits");
			DataCube_add_history(mom2, par);
			DataCube_save(mom2, String_get(filename), overwrite, DESTROY);
		}
		
		if(chan != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_chan.fits");
			DataCube_add_history(chan, par);
			DataCube_save(chan, String_get(filename), overwrite, DESTROY);
		}
		
		if(snr != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_snr.fits");
			DataCube_add_history(snr, par);
			DataCube_save(snr, String_get(filename), overwrite, DESTROY);
		}
		
		if(pv != NULL)
		{
			String_set(filename, String_get(filename_template));
			String_append_int(filename, "%ld", src_id);
			String_append(filename, "_pv.fits");
			DataCube_add_history(pv, par);
			DataCube_save(pv, String_get(filename), overwrite, DESTROY);
		}
		
		// ...spectrum
		String_set(filename, String_get(filename_template));
		String_append_int(filename, "%ld", src_id);
		String_append(filename, "_spec.txt");
		message("Creating text file: %s", strrchr(String_get(filename), '/') == NULL ? String_get(filename) : strrchr(String_get(filename), '/') + 1);
		
		// Get current date and time
		char current_time_string[64];
		time_t current_time = time(NULL);
		strftime(current_time_string, 64, "%a, %d %b %Y, %H:%M:%S", localtime(&current_time));
		
		FILE *fp;
		if(overwrite) fp = fopen(String_get(filename), "wb");
		else fp = fopen(String_get(filename), "wxb");
		ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", String_get(filename));
		
		fprintf(fp, "# Integrated source spectrum\n");
		fprintf(fp, "# Creator: %s (%s)\n", SOFIA_VERSION_FULL, SOFIA_CREATION_DATE);
		fprintf(fp, "# Time:    %s\n", current_time_string);
		fprintf(fp, "#\n");
		fprintf(fp, "# Description of parameters:\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# - chan    Spectral channel number (zero-based).\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# - v_opt / v_rad / v_app\n");
		fprintf(fp, "#           Radial velocity corresponding to the channel number as\n");
		fprintf(fp, "#           described by the WCS information in the header. The suf-\n");
		fprintf(fp, "#           fix denotes optical, radio or apparent radial velocity,\n");
		fprintf(fp, "#           respectively.\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# - freq    Frequency corresponding to the channel number as described\n");
		fprintf(fp, "#           by the WCS information in the header.\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# - f_sum   Sum of flux density values of all spatial pixels covered\n");
		fprintf(fp, "#           by the source in that channel. If the unit is Jy, then\n");
		fprintf(fp, "#           the flux density has already been corrected for the solid\n");
		fprintf(fp, "#           angle of the beam. If instead the unit is Jy/beam, you\n");
		fprintf(fp, "#           will need to manually divide by the beam area which, for\n");
		fprintf(fp, "#           Gaussian beams, will be\n");
		fprintf(fp, "#\n");
		fprintf(fp, "#             pi * a * b / (4 * ln(2))\n");
		fprintf(fp, "#\n");
		fprintf(fp, "#           where a and b are the major and minor axis of the beam in\n");
		fprintf(fp, "#           units of pixels.\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# - n_pix   Number of spatial pixels covered by the source in that\n");
		fprintf(fp, "#           channel. This can be used to determine the statistical\n");
		fprintf(fp, "#           uncertainty of the summed flux value. Again, this has\n");
		fprintf(fp, "#           not yet been corrected for any potential spatial correla-\n");
		fprintf(fp, "#           tion of pixels due to the beam solid angle.\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# Note that a WCS-related column will only be present if WCS conversion\n");
		fprintf(fp, "# was explicitly requested when running the pipeline.\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# Header rows:\n");
		fprintf(fp, "#   1 = column number\n");
		fprintf(fp, "#   2 = parameter name\n");
		fprintf(fp, "#   3 = parameter unit\n");
		fprintf(fp, "#\n");
		if(use_wcs)
		{
			// 4 columns (chan, wcs, f_sum, n_pix)
			fprintf(fp, "#%*s%*s%*s%*s\n", 9, "1",    18, "2",                    18, "3",                   10, "4");
			fprintf(fp, "#%*s%*s%*s%*s\n", 9, "chan", 18, String_get(label_spec), 18, "f_sum",               10, "n_pix");
			fprintf(fp, "#%*s%*s%*s%*s\n", 9, "-",    18, String_get(unit_spec),  18, String_get(unit_flux), 10, "-");
		}
		else
		{
			// 3 columns (chan, f_sum, n_pix)
			fprintf(fp, "#%*s%*s%*s\n", 9, "1",    18, "2",                   10, "3");
			fprintf(fp, "#%*s%*s%*s\n", 9, "chan", 18, "f_sum",               10, "n_pix");
			fprintf(fp, "#%*s%*s%*s\n", 9, "-",    18, String_get(unit_flux), 10, "-");
		}
		fprintf(fp, "#\n");
		
		for(size_t j = 0; j < nz; ++j)
		{
			// Convert z to WCS if requested and possible
			if(use_wcs)
			{
				double spectral = 0.0;
				WCS_convertToWorld(wcs, 0, 0, j + z_min, NULL, NULL, &spectral);
				fprintf(fp, "%*zu%*.7e%*.7e%*zu\n", 10, j + z_min + offset_z, 18, spectral, 18, spectrum[j] / beam_area, 10, pixcount[j]);
			}
			else fprintf(fp, "%*zu%*.7e%*zu\n", 10, j + z_min + offset_z, 18, spectrum[j] / beam_area, 10, pixcount[j]);
		}
		
		fclose(fp);
		
		// Delete output products again
		DataCube_delete(cubelet);
		DataCube_delete(masklet);
		DataCube_delete(mom0);
		DataCube_delete(mom1);
		DataCube_delete(mom2);
		DataCube_delete(chan);
		DataCube_delete(snr);
		DataCube_delete(pv);
		free(spectrum);
		free(pixcount);
	}
	
	// Clean up
	String_delete(filename_template);
	String_delete(filename);
	String_delete(unit_flux_dens);
	String_delete(unit_flux);
	String_delete(unit_spec);
	String_delete(label_spec);
	WCS_delete(wcs);
	
	return;
}



/// @brief Extract beam solid angle from header
///
/// Private method for extracting the solid angle of the beam from
/// the header, reading the `BMAJ`, `BMIN` and `CDELT2` keywords. This
/// assumes that the beam is Gaussian, and the solid angle will be
/// calculated in units of pixels under the assumption that the
/// units of `BMAJ`, `BMIN` and `CDELT2` are the same. If the beam
/// cannot be determined, `NaN` will instead be returned.
///
/// @param self  Object self-reference.
///
/// @return Solid angle of the beam in pixels.

PRIVATE double DataCube_get_beam_area(const DataCube *self)
{
	// Extract beam information
	// WARNING: We assume here that BMAJ, BMIN and CDELT2 have the same unit!
	const double beam_maj   = Header_get_flt(self->header, "BMAJ");
	const double beam_min   = Header_get_flt(self->header, "BMIN");
	const double pixel_size = Header_get_flt(self->header, "CDELT2");
	
	if(IS_NAN(beam_maj) || IS_NAN(beam_min) || IS_NAN(pixel_size) || beam_maj == 0.0 || beam_min == 0.0 || pixel_size == 0.0)
	{
		warning("Failed to determine beam size from header.");
		return NAN;
	}
	
	message("Assuming beam size of %.1f x %.1f pixels.\n", beam_maj / pixel_size, beam_min / pixel_size);
	return M_PI * beam_maj * beam_min / (4.0 * log(2.0) * pixel_size * pixel_size);
}




/// @brief Extract WCS information from header
///
/// Public method for extracting WCS information from the header
/// of the data cube pointed to by `self`. If valid WCS information
/// was found, a WCS object will be returned; otherwise, the method
/// will return a `NULL` pointer.
///
/// @param self  Object self-reference.
///
/// @return Pointer to WCS object if valid, NULL otherwise.
///
/// @note It is the responsibility of the user to call the destructor
///       on the returned WCS object once it is no longer required in
///       order to release its memory.

PUBLIC WCS *DataCube_extract_wcs(const DataCube *self)
{
	WCS *wcs = NULL;
	int *dim_axes = (int *)memory(MALLOC, self->dimension, sizeof(int));  // NOTE: WCSlib requires int!
	
	for(size_t i = 0; i < self->dimension; ++i) dim_axes[i] = (i < 4 && self->axis_size[i]) ? self->axis_size[i] : 1;
	wcs = WCS_new(Header_get(self->header), Header_get_size(self->header) / FITS_HEADER_LINE_SIZE, self->dimension, dim_axes);
	
	free(dim_axes);
	
	return WCS_is_valid(wcs) ? wcs : NULL;
}



/// @brief Swap byte order of data array
///
/// Private method for swapping the byte order of the data array
/// stored in the object referred to by `self`. The function will
/// check if byte order swapping is necessary and, if so, loop over
/// the entire array and call the corresponding swapping function
/// swap_byte_order() defined in common.c on each array element.
///
/// @param self  Object self-reference.

PRIVATE void DataCube_swap_byte_order(const DataCube *self)
{
	if(is_little_endian() && self->word_size > 1)
	{
		#pragma omp parallel for schedule(static)
		for(size_t i = 0; i < self->data_size * self->word_size; i += self->word_size)
		{
			swap_byte_order(self->data + i, self->word_size);
		}
	}
	
	return;
}
