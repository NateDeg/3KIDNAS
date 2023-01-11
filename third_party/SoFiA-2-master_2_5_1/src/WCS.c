// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (WCS.c) - Source Finding Application                     //
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

/// @file   WCS.c
/// @author Tobias Westmeier
/// @date   25/11/2021
/// @brief  Class for handling and converting world coordinates.


#include <stdlib.h>

#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>

#include "WCS.h"



/// @brief Class for handling and converting world coordinates.
///
/// The purpose of this class is to provide a support for World
/// Coordinate System (WCS) conversions in the form of a wrapper around
/// the wcslib package. The class provides methods for setting up WCS
/// information from a FITS file header and converting between pixel
/// coordinates (x, y, z) and world coordinates (lon, lat, spec) in
/// three dimensions.

CLASS WCS
{
	bool valid;               ///< WCS correctly initialised?
	struct wcsprm *wcs_pars;  ///< wcsprm struct from wcslib.
	int  n_wcs_rep;           ///< Number of WCS representations (required by wcslib).
};



/// @brief Standard constructor
///
/// Standard constructor. Will create a new WCS object and return a
/// pointer to the newly created object. The object will be set up
/// from the information in the FITS header supplied by the user.
/// If initialisation of the WCS information fails, the property
/// `valid` will be set to `false`. Note that the destructor will
/// need to be called explicitly once the object is no longer
/// required to release any memory allocated during the lifetime
/// of the object.
///
/// @param header    String containing the raw header information of
///                  the FITS data cube (not null-terminated).
/// @param n_keys    Number of header keywords. Alternatively, it is
///                  also possible to provide the number of header
///                  lines here, as anything beyond the `END` keyword
///                  will be ignored by wcslib anyway.
/// @param n_axes    Number of WCS axes in the data cube.
/// @param dim_axes  Array holding the size of each axis.
///
/// @return Pointer to newly created WCS object.

PUBLIC WCS *WCS_new(const char *header, const int n_keys, const int n_axes, const int *dim_axes)
{
	// Sanity checks
	check_null(header);
	check_null(dim_axes);
	ensure(n_axes, ERR_USER_INPUT, "Failed to set up WCS; FITS header has no WCS axes.");
	
	// Create new WCS object
	WCS *self = (WCS *)memory(MALLOC, 1, sizeof(WCS));
	
	// Initialise properties
	self->valid = false;
	self->wcs_pars = NULL;
	self->n_wcs_rep = 0;
	
	// Set up WCS object from header information
	WCS_setup(self, header, n_keys, n_axes, dim_axes);
	
	// Return new WCS object
	return self;
}



/// @brief Destructor
///
/// Destructor. Note that the destructor must be called explicitly
/// if the object is no longer required. This will release the
/// memory occupied by the object.
///
/// @param self  Object self-reference.

PUBLIC void WCS_delete(WCS *self)
{
	if(WCS_is_valid(self))
	{
		wcsvfree(&self->n_wcs_rep, &self->wcs_pars);
		wcsfree(self->wcs_pars);
		free(self->wcs_pars);
	}
	
	free(self);
	
	return;
}



/// @brief Check if WCS information valid
///
/// Public method for checking if the WCS object has been correctly
/// set up and contains valid WCS information. This should be used
/// before any coordinate conversion is attempted to avoid termination
/// of the process with an error message.
///
/// @param self  Object self-reference.
///
/// @return `true` if WCS is valid, `false` otherwise.

PUBLIC bool WCS_is_valid(const WCS *self)
{
	return self != NULL && self->wcs_pars != NULL && self->valid;
}



/// @brief Set up WCS information from FITS header
///
/// Private method for setting up the WCS object with information
/// from the supplied FITS header. If extraction of WCS information
/// fails for some reason, the object will be left in an invalid
/// state; the property `valid` will be set to `false` in this case,
/// which can later be checked by calling WCS_is_valid().
///
/// @param self      Object self-reference.
/// @param header    C string containing the raw header information of
///                  the FITS data cube (not null-terminated).
/// @param n_keys    Number of header keywords. Alternatively, it is
///                  also possible to provide the number of header
///                  lines here, as anything beyond the `END` keyword
///                  will be ignored by wcslib.
/// @param n_axes    Number of WCS axes in the data cube.
/// @param dim_axes  Array holding the size of each axis.

PRIVATE void WCS_setup(WCS *self, const char *header, const int n_keys, const int n_axes, const int *dim_axes)
{
	// Some variables needed by wcslib
	int n_rejected = 0;
	int stat[NWCSFIX];
	
	// Allocate memory for wcsprm structure
	self->wcs_pars = (struct wcsprm *)memory(CALLOC, 1, sizeof(struct wcsprm));
	self->wcs_pars->flag = -1;
	
	// Initialise wcsprm structure
	int status = wcsini(true, n_axes, self->wcs_pars);
	
	// Parse the FITS header to fill in the wcsprm structure
	if(!status) status = wcspih((char *)header, n_keys, WCSHDR_all, 0, &n_rejected, &self->n_wcs_rep, &self->wcs_pars);
	// NOTE: The (char *) cast is necessary as wcspih would actually
	//       manipulate the header if the 4th argument was negative!
	
	// Apply all necessary corrections to wcsprm structure
	// (missing cards, non-standard units or spectral types, etc.)
	if(!status) status = wcsfix(1, dim_axes, self->wcs_pars, stat);
	
	// Set up additional parameters in wcsprm structure derived from imported data
	if(!status) status = wcsset(self->wcs_pars);
	
	// Redo the corrections to account for things like NCP projections
	if(!status) status = wcsfix(1, dim_axes, self->wcs_pars, stat);
	
	if(status)
	{
		warning("wcslib error %d: %s\n         Failed to parse WCS information.", status, wcs_errmsg[status]);
		self->valid = false;
	}
	else
	{
		message("WCS setup successful.\n");
		self->valid = true;
	}
	
	return;
}



/// @brief Convert from pixel to world coordinates
///
/// Public method for converting the pixel coordinates (x, y, z) to
/// world coordinates (longitude, latitude, spectral). Note that
/// the implicit assumption is made that the first up-to-three axes
/// of the cube are in the aforementioned order. Pixel coordinates
/// must be zero-based; world coordinates will be in the native
/// units of the data cube. Longitude, latitude or spectral can be
/// `NULL`, in which case they are not updated. If invalid input
/// coordinates are supplied by the user, then a warning message will
/// be printed and the output coordinate variables will be left
/// unchanged.
///
/// @param self       Object self-reference.
/// @param x          x coordinate (0-based).
/// @param y          y coordinate (0-based).
/// @param z          z coordinate (0-based).
/// @param longitude  Pointer for holding longitude coordinate.
/// @param latitude   Pointer for holding latitude coordinate.
/// @param spectral   Pointer for holding spectral coordinate.

PUBLIC void WCS_convertToWorld(const WCS *self, const double x, const double y, const double z, double *longitude, double *latitude, double *spectral)
{
	// Sanity checks
	ensure(WCS_is_valid(self), ERR_USER_INPUT, "Failed to convert coordinates; no valid WCS definition found.");
	
	// Determine number of WCS axes
	const size_t n_axes = self->wcs_pars->naxis;
	ensure(n_axes, ERR_USER_INPUT, "Failed to convert coordinates; no valid WCS axes found.");
	
	// Allocate memory for coordinate arrays
	double *coord_pixel = (double *)memory(MALLOC, n_axes, sizeof(double));
	double *coord_world = (double *)memory(MALLOC, n_axes, sizeof(double));
	double *tmp_world   = (double *)memory(MALLOC, n_axes, sizeof(double));
	
	// Initialise pixel coordinates
	for(size_t i = 0; i < n_axes; ++i)
	{
		// NOTE: WCS pixel arrays are 1-based!!!
		if(i == 0)      coord_pixel[i] = 1.0 + x;
		else if(i == 1) coord_pixel[i] = 1.0 + y;
		else if(i == 2) coord_pixel[i] = 1.0 + z;
		else            coord_pixel[i] = 1.0;
	}
	
	// Declare a few variables
	double phi;
	double theta;
	int stat;
	
	// Call WCS conversion module
	int status = wcsp2s(self->wcs_pars, 1, n_axes, coord_pixel, tmp_world, &phi, &theta, coord_world, &stat);
	if(!status)
	{
		// Pass back world coordinates
		if(n_axes > 0 && longitude != NULL) *longitude = coord_world[0];
		if(n_axes > 1 &&  latitude != NULL) *latitude  = coord_world[1];
		if(n_axes > 2 &&  spectral != NULL) *spectral  = coord_world[2];
	}
	else warning("wcslib error %d: %s", status, wcs_errmsg[status]);
	
	// Clean up
	free(coord_pixel);
	free(coord_world);
	free(tmp_world);
	
	return;
}



/// @brief Convert from world to pixel coordinates
///
/// Public method for converting the world coordinates (longitude,
/// latitude, spectral) to world coordinates (x, y, z). Note that
/// the implicit assumption is made that the first up-to-three axes
/// of the cube are in the aforementioned order. Pixel coordinates
/// will be zero-based; world coordinates must be in the native
/// units of the data cube. X, y or z can be `NULL`, in which case
/// they are not updated. If invalid input coordinates are supplied
/// by the user, then a warning message will be printed and the
/// output coordinate variables will be left unchanged.
///
/// @param self       Object self-reference.
/// @param longitude  Longitude coordinate.
/// @param latitude   Latitude coordinate.
/// @param spectral   Spectral coordinate.
/// @param x          Pointer for holding x coordinate (0-based).
/// @param y          Pointer for holding y coordinate (0-based).
/// @param z          Pointer for holding z coordinate (0-based).

PUBLIC void WCS_convertToPixel(const WCS *self, const double longitude, const double latitude, const double spectral, double *x, double *y, double *z)
{
	// Sanity checks
	ensure(WCS_is_valid(self), ERR_USER_INPUT, "Failed to convert coordinates; no valid WCS definition found.");
	
	// Determine number of WCS axes
	const size_t n_axes = self->wcs_pars->naxis;
	ensure(n_axes, ERR_USER_INPUT, "Failed to convert coordinates; no valid WCS axes found.");
	
	// Allocate memory for coordinate arrays
	double *coord_pixel = (double *)memory(MALLOC, n_axes, sizeof(double));
	double *coord_world = (double *)memory(MALLOC, n_axes, sizeof(double));
	double *tmp_world   = (double *)memory(MALLOC, n_axes, sizeof(double));
	
	// Initialise pixel coordinates
	for(size_t i = 0; i < n_axes; ++i)
	{
		if(i == 0)      coord_world[i] = longitude;
		else if(i == 1) coord_world[i] = latitude;
		else if(i == 2) coord_world[i] = spectral;
		else            coord_world[i] = 0.0;
	}
	
	// Declare a few variables
	double phi;
	double theta;
	int stat;
	
	int status = wcss2p(self->wcs_pars, 1, n_axes, coord_world, &phi, &theta, tmp_world, coord_pixel, &stat);
	if(!status)
	{
		// Pass back pixel coordinates
		// NOTE: WCS pixel arrays are 1-based!!!
		if(n_axes > 0 && x != NULL) *x = coord_pixel[0] - 1.0;
		if(n_axes > 1 && y != NULL) *y = coord_pixel[1] - 1.0;
		if(n_axes > 2 && z != NULL) *z = coord_pixel[2] - 1.0;
	}
	else warning("wcslib error %d: %s", status, wcs_errmsg[status]);
	
	// Clean up
	free(coord_pixel);
	free(coord_world);
	free(tmp_world);
	
	return;
}
