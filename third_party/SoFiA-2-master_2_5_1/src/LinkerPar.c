// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (LinkerPar.c) - Source Finding Application               //
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

/// @file   LinkerPar.c
/// @author Tobias Westmeier
/// @date   03/08/2021
/// @brief  Class for managing the parameters of detections formed by the linker in SoFiA 2.


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "LinkerPar.h"
#include "String.h"
#include "statistics_dbl.h"

/// Set to 1 if measurement of the flux-weighted centroid is required.
#define MEASURE_CENTROID_POSITION 0

/// Set to 1 if additional statistical parameters are required in the output catalogue.
#define MEASURE_ADDITIONAL_STATS  1


/// @brief Class for handling detections created by the linker
///
/// The purpose of this class is to provide a structure for storing
/// and updating source parameters handled by the linker implemented
/// in the DataCube class.

CLASS LinkerPar
{
	size_t  size;          ///< Number of detections currently stored.
	int     verbosity;     ///< Verbosity level (0 or 1).
	size_t *label;         ///< Label of the detection.
	size_t *n_pix;         ///< Number of pixels contained in detection mask.
	size_t *x_min;         ///< Lower end of bounding box in x.
	size_t *x_max;         ///< Upper end of bounding box in x.
	size_t *y_min;         ///< Lower end of bounding box in y.
	size_t *y_max;         ///< Upper end of bounding box in y.
	size_t *z_min;         ///< Lower end of bounding box in z.
	size_t *z_max;         ///< Upper end of bounding box in z.
	#if MEASURE_CENTROID_POSITION
	double *x_ctr;         ///< Centroid position in x.
	double *y_ctr;         ///< Centroid position in y.
	double *z_ctr;         ///< Centroid position in z.
	#endif
	double *f_min;         ///< Minimum flux density within mask.
	double *f_max;         ///< Maximum flux density within mask.
	double *f_sum;         ///< Sum of flux densities within mask.
	double *rel;           ///< Reliability of detection.
	unsigned char *flags;  ///< Quality flags.
	double *fill;          ///< Filing factor of mask within bounding box.
	double *m1;            ///< 1st moment of flux densities within mask.
	double *m2;            ///< 2nd moment of flux densities within mask.
	double *m3;            ///< 3rd moment of flux densities within mask.
	double *m4;            ///< 4th moment of flux densities within mask.
};


/// @brief Standard constructor
///
/// Will create a new and empty LinkerPar object and return a pointer to the newly
/// created object. No memory will be allocated other than for the object itself.
/// Note that the destructor will need to be called explicitly once the object is
/// no longer required to release any memory allocated during the lifetime of the
/// object.
///
/// @param verbosity  Verbosity level of the new object.
///
/// @return Pointer to newly created LinkerPar object.

PUBLIC LinkerPar *LinkerPar_new(const bool verbosity)
{
	LinkerPar *self = (LinkerPar *)memory(MALLOC, 1, sizeof(LinkerPar));
	
	self->verbosity = verbosity;
	self->size = 0;
	
	self->label = NULL;
	self->n_pix = NULL;
	self->x_min = NULL;
	self->x_max = NULL;
	self->y_min = NULL;
	self->y_max = NULL;
	self->z_min = NULL;
	self->z_max = NULL;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr = NULL;
	self->y_ctr = NULL;
	self->z_ctr = NULL;
	#endif
	self->f_min = NULL;
	self->f_max = NULL;
	self->f_sum = NULL;
	self->rel   = NULL;
	self->flags = NULL;
	self->fill  = NULL;
	self->m1    = NULL;
	self->m2    = NULL;
	self->m3    = NULL;
	self->m4    = NULL;
	
	return self;
}


/// @brief Destructor
///
/// Note that the destructor must be called explicitly if the object is no
/// longer required. This will release the memory occupied by the object.
///
/// @param self  Object self-reference

PUBLIC void LinkerPar_delete(LinkerPar *self)
{
	if(self != NULL)
	{
		self->size = 0;
		LinkerPar_reallocate_memory(self);
		free(self);
	}
	
	return;
}


/// @brief Return number of sources in LinkerPar object
///
/// Public method for returning the size of the LinkerPar object,
/// i.e. the number of sources it currently contains.
///
/// @param self  Object self-reference.
///
/// @return Number of sources stored in LinkerPar object.

PUBLIC size_t LinkerPar_get_size(const LinkerPar *self)
{
	return self == NULL ? 0 : self->size;
}


/// @brief Insert new object at end of list
///
/// Public method for adding a new object to the end of the list pointed to by
/// `self`. The label will be set to `label`, the number of pixels to 1, and
/// the (x, y, z) position will be used as the initial `x_min`, `x_max`, etc.
/// The memory allocation of the object will automatically be expanded if necessary.
///
/// @param self   Object self-reference.
/// @param label  Label of the new object.
/// @param x      x-position of the new object.
/// @param y      y-position of the new object.
/// @param z      z-position of the new object.
/// @param flux   Flux value of the new object.
/// @param flag   Flag values to be set; 1 = spatial edge; 2 = spectral edge; 4 = blanked pixels; 8 = other sources.

PUBLIC void LinkerPar_push(LinkerPar *self, const size_t label, const size_t x, const size_t y, const size_t z, const double flux, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	
	// Increment size counter
	++self->size;
	
	// Allocate additional memory
	LinkerPar_reallocate_memory(self);
	
	// Insert new element at end
	self->label[self->size - 1] = label;
	self->n_pix[self->size - 1] = 1;
	self->x_min[self->size - 1] = x;
	self->x_max[self->size - 1] = x;
	self->y_min[self->size - 1] = y;
	self->y_max[self->size - 1] = y;
	self->z_min[self->size - 1] = z;
	self->z_max[self->size - 1] = z;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr[self->size - 1] = flux * x;
	self->y_ctr[self->size - 1] = flux * y;
	self->z_ctr[self->size - 1] = flux * z;
	#endif
	self->f_min[self->size - 1] = flux;
	self->f_max[self->size - 1] = flux;
	self->f_sum[self->size - 1] = flux;
	self->rel  [self->size - 1] = 0.0;  // NOTE: Must be 0 (default for neg. sources), as only pos. sources will be updated later!
	self->flags[self->size - 1] = flag;
	self->fill[self->size - 1]  = 1;
	self->m1[self->size - 1]    = flux;
	self->m2[self->size - 1]    = 0.0;
	self->m3[self->size - 1]    = 0.0;
	self->m4[self->size - 1]    = 0.0;
	
	return;
}


/// @brief Remove last object from list
///
/// Public method for removing the most recently added object from the list.
/// The process will be terminated if the original list is empty.
///
/// @param self  Object self-reference.

PUBLIC void LinkerPar_pop(LinkerPar *self)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_FAILURE, "Failed to pop element from empty LinkerPar object.");
	
	// Decrement size
	--self->size;
	
	// Reallocate memory
	LinkerPar_reallocate_memory(self);
	
	return;
}


/// @brief Add another pixel to last object in list
///
/// Public method for adding another pixel to the last object in the current linker list.
/// The object's `x_min`, `x_max`, `y_min`, etc. values will be checked against the
/// newly added pixel and updated if necessary. The programme will terminate if the list
/// is found to be empty.
///
/// @param self  Object self-reference.
/// @param x     x-position of the new pixel.
/// @param y     y-position of the new pixel.
/// @param z     z-position of the new pixel.
/// @param flux  Flux value of the new pixel.
/// @param flag  Flag values to be set; 1 = spatial edge; 2 = spectral edge; 4 = blanked
///              pixels; 8 = other sources.

PUBLIC void LinkerPar_update(LinkerPar *self, const size_t x, const size_t y, const size_t z, const double flux, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_USER_INPUT, "Failed to update LinkerPar object; list is currently empty.");
	
	// Get index
	const size_t index = self->size - 1;
	
	++self->n_pix[index];
	if(x < self->x_min[index]) self->x_min[index] = x;
	if(x > self->x_max[index]) self->x_max[index] = x;
	if(y < self->y_min[index]) self->y_min[index] = y;
	if(y > self->y_max[index]) self->y_max[index] = y;
	if(z < self->z_min[index]) self->z_min[index] = z;
	if(z > self->z_max[index]) self->z_max[index] = z;
	#if MEASURE_CENTROID_POSITION
	self->x_ctr[index] += flux * x;
	self->y_ctr[index] += flux * y;
	self->z_ctr[index] += flux * z;
	#endif
	if(flux > self->f_max[index]) self->f_max[index] = flux;
	if(flux < self->f_min[index]) self->f_min[index] = flux;
	self->f_sum[index] += flux;
	self->flags[index] |= flag;
	self->fill[index] = (double)(self->n_pix[index]) / (double)((self->x_max[index] - self->x_min[index] + 1) * (self->y_max[index] - self->y_min[index] + 1) * (self->z_max[index] - self->z_min[index] + 1));
	const size_t n = self->n_pix[index];
	const size_t n1 = self->n_pix[index] - 1;
	const double delta = flux - self->m1[index];
	const double delta_n = delta / n;
	const double delta_n2 = delta_n * delta_n;
	const double term1 = delta * delta_n * n1;
	self->m1[index] += delta_n;
	self->m4[index] += term1 * delta_n2 * (n * n - 3 * n + 3) + 6.0 * delta_n2 * self->m2[index] - 4.0 * delta_n * self->m3[index];
	self->m3[index] += term1 * delta_n * (n - 2) - 3.0 * delta_n * self->m2[index];
	self->m2[index] += term1;
	
	return;
}

/// @brief Update flag of last object in current list
///
/// Public method for updating the flag of the last object in the current linker list.
/// The programme will terminate if the list is found to be empty.
///
/// @param self  Object self-reference.
/// @param flag  Flag value to be set; 1 = spatial edge; 2 = spectral edge; 4 = blanked
///              pixels; 8 = other sources.

PUBLIC void LinkerPar_update_flag(LinkerPar *self, const unsigned char flag)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_USER_INPUT, "Failed to update LinkerPar object; list is currently empty.");
	
	// Set flag
	self->flags[self->size - 1] |= flag;
	
	return;
}


/// @brief Get size of object in x, y or z
///
/// Public method for returning the size of the specified object along the specified axis.
/// The programme will terminate if the `axis` or `label` are out of range.
///
/// @param self   Object self-reference.
/// @param label  Index of the object to be retrieved.
/// @param axis   Axis for which size should be returned; 0 = x, 1 = y, and 2 = z.
///
/// @return Size of the object in pixels along the specified axis.

PUBLIC size_t LinkerPar_get_obj_size(const LinkerPar *self, const size_t label, const int axis)
{
	// Sanity checks
	check_null(self);
	ensure(axis >= 0 && axis <= 2, ERR_USER_INPUT, "Invalid axis selection (%d) in LinkerPar object.", axis);
	
	// Determine index
	const size_t index = LinkerPar_get_index(self, label);
	
	if(axis == 0) return self->x_max[index] - self->x_min[index] + 1;
	if(axis == 1) return self->y_max[index] - self->y_min[index] + 1;
	return self->z_max[index] - self->z_min[index] + 1;
}


/// @brief Get number of pixels of object
///
/// Public method for returning the number of pixels that have been recorded for the specified object.
/// The programme will terminate if the label is out of range.
///
/// @param self   Object self-reference.
/// @param label  Index of the object to be retrieved.
///
/// @return Number of pixels of the specified object.

PUBLIC size_t LinkerPar_get_npix(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->n_pix[LinkerPar_get_index(self, label)];
}


/// @brief Get total flux of object
///
/// Public method for returning the total flux of the object specified by `label`.
/// The programme will terminate if the label is out of range.
///
/// @param self   Object self-reference.
/// @param label  Index of the object to be retrieved.
///
/// @return Total flux of the specified object.

PUBLIC double LinkerPar_get_flux(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->f_sum[LinkerPar_get_index(self, label)];
}


/// @brief Get reliability of object
///
/// Public method for returning the reliability of the object specified by `label`.
/// The programme will terminate if the label is out of range.
///
/// @param self   Object self-reference.
/// @param label  Index of the object to be retrieved.
///
/// @return Reliability of the specified object.

PUBLIC double LinkerPar_get_rel(const LinkerPar *self, const size_t label)
{
	check_null(self);
	return self->rel[LinkerPar_get_index(self, label)];
}


/// @brief Get label of object by index
///
/// Public method for returning the label of the object with the specified `index`.
/// The programme will terminate if the `index` is out of range.
///
/// @param self   Object self-reference.
/// @param index  Index of the object to be retrieved.
///
/// @return Reliability of the specified object.

PUBLIC size_t LinkerPar_get_label(const LinkerPar *self, const size_t index)
{
	// Sanity checks
	check_null(self);
	ensure(index < self->size, ERR_INDEX_RANGE, "Index out of range. Cannot retrieve label.");
	
	return self->label[index];
}


/// @brief Get bounding box of object
///
/// Public method for retrieving the bounding box of the object with the specified label.
/// The process will be terminated if the label does not exist.
///
/// @param self   Object self-reference.
/// @param label  Index of the object to be retrieved.
/// @param x_min  Pointer for holding the bounding box x_min value.
/// @param x_max  Pointer for holding the bounding box x_max value.
/// @param y_min  Pointer for holding the bounding box y_min value.
/// @param y_max  Pointer for holding the bounding box y_max value.
/// @param z_min  Pointer for holding the bounding box z_min value.
/// @param z_max  Pointer for holding the bounding box z_max value.

PUBLIC void LinkerPar_get_bbox(const LinkerPar *self, const size_t label, size_t *x_min, size_t *x_max, size_t *y_min, size_t *y_max, size_t *z_min, size_t *z_max)
{
	// Sanity checks
	check_null(self);
	
	// Determine index
	size_t index = LinkerPar_get_index(self, label);
	
	// Copy bounding box values
	*x_min = self->x_min[index];
	*x_max = self->x_max[index];
	*y_min = self->y_min[index];
	*y_max = self->y_max[index];
	*z_min = self->z_min[index];
	*z_max = self->z_max[index];
	
	return;
}


/// @brief Create source catalogue from LinkerPar object
///
/// Public method for generating a source catalogue from the specified LinkerPar
/// object. A pointer to the newly created catalogue will be returned. Note that
/// the user will assume ownership of the catalogue and will be responsible for
/// explicitly calling the destructor if the catalogue is no longer required.
/// Unreliable sources can be excluded from the catalogue by providing a non-empty
/// filter object that contains the old and new labels of all sources deemed reliable.
///
/// The old labels will be replaced by the new ones in the source ID column of the
/// catalogue. If filter is `NULL` or empty, all sources will be copied into the
/// catalogue without filtering.
///
/// @param self       Object self-reference.
/// @param filter     Map object containing old and new label pairs of only those
///                   sources considered as reliable.
/// @param flux_unit  String containing the flux unit of the data.
///
/// @return Pointer to newly created catalogue.

PUBLIC Catalog *LinkerPar_make_catalog(const LinkerPar *self, const Map *filter, const char *flux_unit)
{
	// Sanity checks
	check_null(self);
	check_null(flux_unit);
	
	// Check if reliability filtering requested
	const bool remove_unreliable = filter != NULL && Map_get_size(filter);
	
	// Create an empty source catalogue
	Catalog *cat = Catalog_new();
	
	// Create string for holding identifier
	String *identifier = String_new("");
	
	// Loop over all LinkerPar entries
	for(size_t i = 0; i < self->size; ++i)
	{
		const size_t new_label = remove_unreliable && Map_key_exists(filter, self->label[i]) ? Map_get_value(filter, self->label[i]) : self->label[i];
		
		if(!remove_unreliable || Map_key_exists(filter, self->label[i]))
		{
			// Create a new source
			Source *src = Source_new(self->verbosity);
			
			// Set the identifier to the current label
			String_set_int(identifier, "%zu", self->label[i]);
			Source_set_identifier(src, String_get(identifier));
			
			// Set other parameters
			Source_add_par_int(src, "id",    new_label,                       "",        "meta.id");
			#if MEASURE_CENTROID_POSITION
			Source_add_par_flt(src, "x",     self->x_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.x");
			Source_add_par_flt(src, "y",     self->y_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.y");
			Source_add_par_flt(src, "z",     self->z_ctr[i] / self->f_sum[i], "pix",     "pos.cartesian.z");
			#else
			// Insert placeholder if no centroid requested
			Source_add_par_flt(src, "x",     0.0,                             "pix",     "pos.cartesian.x");
			Source_add_par_flt(src, "y",     0.0,                             "pix",     "pos.cartesian.y");
			Source_add_par_flt(src, "z",     0.0,                             "pix",     "pos.cartesian.z");
			#endif
			Source_add_par_int(src, "x_min", self->x_min[i],                  "pix",     "pos.cartesian.x;stat.min");
			Source_add_par_int(src, "x_max", self->x_max[i],                  "pix",     "pos.cartesian.x;stat.max");
			Source_add_par_int(src, "y_min", self->y_min[i],                  "pix",     "pos.cartesian.y;stat.min");
			Source_add_par_int(src, "y_max", self->y_max[i],                  "pix",     "pos.cartesian.y;stat.max");
			Source_add_par_int(src, "z_min", self->z_min[i],                  "pix",     "pos.cartesian.z;stat.min");
			Source_add_par_int(src, "z_max", self->z_max[i],                  "pix",     "pos.cartesian.z;stat.max");
			Source_add_par_int(src, "n_pix", self->n_pix[i],                  "",        "meta.number;instr.pixel");
			Source_add_par_flt(src, "f_min", self->f_min[i],                  flux_unit, "phot.flux.density;stat.min");
			Source_add_par_flt(src, "f_max", self->f_max[i],                  flux_unit, "phot.flux.density;stat.max");
			Source_add_par_flt(src, "f_sum", self->f_sum[i],                  flux_unit, "phot.flux");
			Source_add_par_flt(src, "rel",   self->rel  [i],                  "",        "stat.probability");
			Source_add_par_int(src, "flag",  self->flags[i],                  "",        "meta.code.qual");
			#if MEASURE_ADDITIONAL_STATS
			const double n = (double)(self->n_pix[i]);
			const double variance = self->m2[i] / (n - 1.0);
			const double skewness = sqrt(n) * self->m3[i] / pow(self->m2[i], 1.5);
			const double kurtosis = n * self->m4[i] / (self->m2[i] * self->m2[i]) - 3.0;
			Source_add_par_flt(src, "fill",  self->fill[i],                   "",        "stat.filling");
			Source_add_par_flt(src, "mean",  self->m1[i],                     flux_unit, "phot.flux.density;stat.mean");
			Source_add_par_flt(src, "std",   sqrt(variance),                  flux_unit, "phot.flux.density;stat.stdev");
			Source_add_par_flt(src, "skew",  skewness,                        "",        "stat.param");
			Source_add_par_flt(src, "kurt",  kurtosis,                        "",        "stat.param");
			#endif
			
			// Add source to catalogue
			Catalog_add_source(cat, src);
		}
	}
	
	// Clean up
	String_delete(identifier);
	
	// Return catalogue
	return cat;
}


/// @brief Create reliability parameter catalogues from LinkerPar object
///
/// Public method for generating catalogues of relevant reliability
/// parameters of all negative and positive detections for debugging
/// purposes. These include the source ID, geometric centroid (x, y, z),
/// reliability and the three parameters used in the reliability
/// calculation: log(f_max), log(f_sum) and log(f_mean). Two pointers
/// to empty catalogues must be provided; they will be holding
/// the negative and positive detections, respectively.
///
/// @param self             Object self-reference.
/// @param flux_unit        String holding the flux unit of the data.
/// @param cat_rel_par_neg  Pointer to catalogue that will hold the
///                         reliability parameters of all negative
///                         detections.
/// @param cat_rel_par_pos  Pointer to catalogue that will hold the
///                         reliability parameters of all positive
///                         detections.

PUBLIC void LinkerPar_get_rel_cat(const LinkerPar *self, const char *flux_unit, Catalog **cat_rel_par_neg, Catalog **cat_rel_par_pos)
{
	// Sanity checks
	check_null(self);
	check_null(flux_unit);
	check_null(cat_rel_par_neg);
	check_null(cat_rel_par_pos);
	ensure(Catalog_get_size(*cat_rel_par_neg) == 0 && Catalog_get_size(*cat_rel_par_pos) == 0, ERR_USER_INPUT, "Non-empty reliability parameter catalogue provided.");
	
	// Create string for holding identifier
	String *identifier = String_new("");
	
	// Loop over all LinkerPar entries
	for(size_t i = 0; i < self->size; ++i)
	{
		// Check if flux is negative
		const bool is_neg = (self->f_sum[i] < 0.0);
		
		// Create a new source
		Source *src = Source_new(self->verbosity);
		
		// Set identifier according to current label
		String_set_int(identifier, is_neg ? "neg_%zu" : "pos_%zu", self->label[i]);
		Source_set_identifier(src, String_get(identifier));
		
		// Add basic parameters
		Source_add_par_int(src, "id",  self->label[i],                                  "",    "meta.id");
		Source_add_par_flt(src, "x",   (double)(self->x_max[i] + self->x_min[i]) / 2.0, "pix", "pos.cartesian.x");
		Source_add_par_flt(src, "y",   (double)(self->y_max[i] + self->y_min[i]) / 2.0, "pix", "pos.cartesian.y");
		Source_add_par_flt(src, "z",   (double)(self->z_max[i] + self->z_min[i]) / 2.0, "pix", "pos.cartesian.z");
		Source_add_par_flt(src, "nx",  self->x_max[i] - self->x_min[i] + 1,             "pix", "pos.cartesian.x;arith.diff");
		Source_add_par_flt(src, "ny",  self->y_max[i] - self->y_min[i] + 1,             "pix", "pos.cartesian.y;arith.diff");
		Source_add_par_flt(src, "nz",  self->z_max[i] - self->z_min[i] + 1,             "pix", "pos.cartesian.z;arith.diff");
		Source_add_par_flt(src, "rel", self->rel[i],                                    "",    "stat.probability");
		
		// Add relevant reliability parameters
		Source_add_par_flt(src, "log_f_peak", is_neg ? log10(-self->f_min[i]) : log10(self->f_max[i]), flux_unit, "phot.flux.density;stat.max");
		Source_add_par_flt(src, "log_f_sum", is_neg ? log10(-self->f_sum[i]) : log10(self->f_sum[i]), flux_unit, "phot.flux");
		Source_add_par_flt(src, "log_f_mean", is_neg ? log10(-self->f_sum[i] / self->n_pix[i]) : log10(self->f_sum[i] / self->n_pix[i]), flux_unit, "phot.flux;stat.mean");
		
		// Add extra statistical parameters
		const double n = (double)(self->n_pix[i]);
		const double variance = self->m2[i] / (n - 1.0);
		const double skewness = sqrt(n) * self->m3[i] / pow(self->m2[i], 1.5);
		const double kurtosis = n * self->m4[i] / (self->m2[i] * self->m2[i]) - 3.0;
		Source_add_par_flt(src, "fill",  self->fill[i],  "",        "stat.filling");
		Source_add_par_flt(src, "std",   sqrt(variance), flux_unit, "phot.flux.density;stat.stdev");
		Source_add_par_flt(src, "skew",  skewness,       "",        "stat.param");
		Source_add_par_flt(src, "kurt",  kurtosis,       "",        "stat.param");
		
		// Add source to appropriate catalogue
		Catalog_add_source(is_neg ? *cat_rel_par_neg : *cat_rel_par_pos, src);
	}
	
	// Clean up
	String_delete(identifier);
	
	return;
}


/// @brief Print basic information about LinkerPar object
///
/// Public method for printing some basic information on the size and memory
/// usage of the LinkerPar object pointed to by `self`.
///
/// @param self  Object self-reference.

PUBLIC void LinkerPar_print_info(const LinkerPar *self)
{
	// Sanity checks
	check_null(self);
	
	// Calculate memory usage
	#if MEASURE_CENTROID_POSITION
		const double memory_usage = (double)(self->size * (8 * sizeof(size_t) + 12 * sizeof(double) + 1 * sizeof(char)));
	#else
		const double memory_usage = (double)(self->size * (8 * sizeof(size_t) + 9 * sizeof(double) + 1 * sizeof(char)));
	#endif
	
	// Print size and memory information
	message("Linker status:");
	message(" - No. of objects:  %zu", self->size);
	if(memory_usage < MEGABYTE) message(" - Memory usage:    %.2f kB\n", memory_usage / KILOBYTE);
	else message(" - Memory usage:    %.2f MB\n", memory_usage / MEGABYTE);
	
	return;
}


/// @brief Return index of element by label
///
/// Private method for finding and returning the index of the element
/// with the specified label. The process will be terminated if the
/// requested label does not exist.
///
/// @param self   Object self-reference.
/// @param label  Label of the element the index of which is to be returned.
///
/// @return Index of the element identified with `label`

PRIVATE size_t LinkerPar_get_index(const LinkerPar *self, const size_t label)
{
	size_t index = 0;
	while(index < self->size && self->label[index] != label) ++index;
	ensure(self->size && self->label[index] == label, ERR_USER_INPUT, "Label not found.");
	return index;
}


/// @brief Reallocate memory for LinkerPar object
///
/// Private method for reallocating the memory requirements of the
/// specified LinkerPar object, e.g. as necessitated by a change in
/// size. If the new size is 0, all memory will be de-allocated and
/// the pointers will be set to `NULL`.
///
/// @param self  Object self-reference.

PRIVATE void LinkerPar_reallocate_memory(LinkerPar *self)
{
	if(self->size)
	{
		// Reallocate memory
		self->label = (size_t *)memory_realloc(self->label, self->size, sizeof(size_t));
		self->n_pix = (size_t *)memory_realloc(self->n_pix, self->size, sizeof(size_t));
		self->x_min = (size_t *)memory_realloc(self->x_min, self->size, sizeof(size_t));
		self->x_max = (size_t *)memory_realloc(self->x_max, self->size, sizeof(size_t));
		self->y_min = (size_t *)memory_realloc(self->y_min, self->size, sizeof(size_t));
		self->y_max = (size_t *)memory_realloc(self->y_max, self->size, sizeof(size_t));
		self->z_min = (size_t *)memory_realloc(self->z_min, self->size, sizeof(size_t));
		self->z_max = (size_t *)memory_realloc(self->z_max, self->size, sizeof(size_t));
		#if MEASURE_CENTROID_POSITION
		self->x_ctr = (double *)memory_realloc(self->x_ctr, self->size, sizeof(double));
		self->y_ctr = (double *)memory_realloc(self->y_ctr, self->size, sizeof(double));
		self->z_ctr = (double *)memory_realloc(self->z_ctr, self->size, sizeof(double));
		#endif
		self->f_min = (double *)memory_realloc(self->f_min, self->size, sizeof(double));
		self->f_max = (double *)memory_realloc(self->f_max, self->size, sizeof(double));
		self->f_sum = (double *)memory_realloc(self->f_sum, self->size, sizeof(double));
		self->rel   = (double *)memory_realloc(self->rel,   self->size, sizeof(double));
		self->flags = (unsigned char *)memory_realloc(self->flags, self->size, sizeof(unsigned char));
		self->fill  = (double *)memory_realloc(self->fill,  self->size, sizeof(double));
		self->m1    = (double *)memory_realloc(self->m1,    self->size, sizeof(double));
		self->m2    = (double *)memory_realloc(self->m2,    self->size, sizeof(double));
		self->m3    = (double *)memory_realloc(self->m3,    self->size, sizeof(double));
		self->m4    = (double *)memory_realloc(self->m4,    self->size, sizeof(double));
	}
	else
	{
		free(self->label);
		free(self->n_pix);
		free(self->x_min);
		free(self->x_max);
		free(self->y_min);
		free(self->y_max);
		free(self->z_min);
		free(self->z_max);
		#if MEASURE_CENTROID_POSITION
		free(self->x_ctr);
		free(self->y_ctr);
		free(self->z_ctr);
		#endif
		free(self->f_min);
		free(self->f_max);
		free(self->f_sum);
		free(self->rel);
		free(self->flags);
		free(self->fill);
		free(self->m1);
		free(self->m2);
		free(self->m3);
		free(self->m4);
		
		self->label = NULL;
		self->n_pix = NULL;
		self->x_min = NULL;
		self->x_max = NULL;
		self->y_min = NULL;
		self->y_max = NULL;
		self->z_min = NULL;
		self->z_max = NULL;
		#if MEASURE_CENTROID_POSITION
		self->x_ctr = NULL;
		self->y_ctr = NULL;
		self->z_ctr = NULL;
		#endif
		self->f_min = NULL;
		self->f_max = NULL;
		self->f_sum = NULL;
		self->rel   = NULL;
		self->flags = NULL;
		self->fill  = NULL;
		self->m1    = NULL;
		self->m2    = NULL;
		self->m3    = NULL;
		self->m4    = NULL;
	}
	
	return;
}


/// @brief Determine reliability of detections
///
/// Public method for measuring the reliability of all the sources in the specified
/// LinkerPar object. This will set the rel property of the object, but not yet
/// filter out unreliable sources. Reliability measurement works by comparing the
/// density of positive and negative detections in an N-dimensional parameter space.
/// For this purpose, the covariance matrix of the distribution of negative sources
/// in parameter space is first calculated. The covariance matrix is assumed to
/// describe the multivariate normal distribution of the Gaussian noise of the data.
/// Next, the sum of the probability density functions of all positive and negative
/// sources is evaluated at the location of each positive detection (multivariate
/// Gaussian kernel density estimation). From this, the reliability is estimated as
///
/// R = (P - N) / N,
///
/// where P is the sum of the PDFs of the positive sources, and N is the sum
/// of the PDFs of the negative sources. If N > P, R is set to 0 to
/// ensure that the resulting reliability is always in the range of 0 to 1. Note that
/// the reliability will only be determined for positive sources above the `fmin`
/// threshold, where `fmin` is the summed flux divided by the square root of the
/// number of pixels contributing to a source. In order to be able to exclude certain
/// negative artefacts from affecting the reliability calculation, the user has the
/// option of specifying a table of (x, y) pixel positions using the parameter
/// `rel_cat`. All negative detections the (x, y) bounding box of which contains one
/// of those positions will be excluded from the reliability calculation. `rel_cat`
/// must contain exactly two columns (x and y in pixels). If set to `NULL`, this
/// feature will be disabled altogether. This method can also create a Skellam array
/// to assist with the optimisation of the kernel scale. This can be controlled using
/// the `skellam` parameter. If set to `NULL`, no Skellam array will be generated.
///
/// It is also possible to employ an automatic kernel scaling algorithm which will
/// iteratively change the kernel scale factor until convergence is achieved. The 
/// auto-kernel algorithm will start at a low `scale_kernel = 0.1` and will
/// incrementally increase that value until the absolute value of the median of the
/// Skellam distribution drops below a user-defined tolerance or the maximum
/// number of iterations is exceeded. The `scale_kernel` increments will become
/// smaller as the algorithm approaches the threshold. If the algorithm fails to
/// converge, the original value of `scale_kernel` will instead be used. Automatic
/// kernel scaling can be enabled by setting `autokernel` to `true`.
///
/// @param self           Object self-reference.
/// @param rel_par_space  Array of parameters to be used to determine the
///                       reliability of detections. These must be integer
///                       values corresponding to the parameters defined
///                       in the header file.
/// @param scale_kernel   The size of the convolution kernel used in determining
///                       the density of positive and negative detections in
///                       parameter space will be scaled by this factor. If set to 1,
///                       the original covariance matrix derived from the distribution
///                       of negative sources is used. This parameter must be passed by
///                       reference and will hold the new kernel scale factor determined
///                       by the auto-kernel algorithm if enabled.
/// @param minpix         Minimum number of pixels for a source to be considered reliable.
/// @param fmin           Value of the `fmin` parameter, where `fmin = sum / sqrt(N)`.
/// @param rel_cat        Table of pixel coordinates on the sky. All negative detections
///                       with bounding boxes including those positions will be removed
///                       before reliability calculation. NULL can be used to disable
///                       this feature.
/// @param skellam        Pointer to an Array of type double to hold the Skellam array.
///                       Set NULL to disable.
/// @param autokernel     If set to `true` then automatic scaling of the kernel will be
///                       enabled. `scale_kernel` will be used as the fall-back option
///                       if the auto-kernel algorithm fails to converge.
/// @param iterations     Maximum number of iterations used by the auto-kernel algorithm.
///                       If the algorithm does not converge, `scale_kernel` will be used
///                       instead.
/// @param tolerance      Skellam parameter tolerance for convergence of the auto-kernel
///                       algorithm. The algorithm converges when the absolute value of
///                       the median of the Skellam distribution drops below this value.
///
/// @return Covariance matrix from the negative detections.

PUBLIC Matrix *LinkerPar_reliability(LinkerPar *self, const Array_siz *rel_par_space, double *scale_kernel, const double fmin, const size_t minpix, const Table *rel_cat, Array_dbl **skellam, const bool autokernel, const int iterations, const double tolerance)
{
	// Sanity checks
	check_null(self);
	ensure(self->size, ERR_NO_SRC_FOUND, "No sources left after linking. Cannot proceed.");
	ensure(skellam != NULL || !autokernel, ERR_USER_INPUT, "With kernel auto-scaling enabled, skellam must not be NULL.");
	if(*scale_kernel <= 0.0)
	{
		warning("Kernel scale factor is non-positive; using default of 0.4 instead.");
		*scale_kernel = 0.4;
	}
	
	// Dimensionality of parameter space
	const int dim = Array_siz_get_size(rel_par_space);
	const char *par_names[] = {"peak", "sum", "mean", "chan", "pix", "fill", "std", "skew", "kurt"};
	message("Using %zuD parameter space:", Array_siz_get_size(rel_par_space));
	for(size_t i = 0; i < Array_siz_get_size(rel_par_space); ++i) message(" - %s", par_names[Array_siz_get(rel_par_space, i)]);
	
	// Define a few parameters
	size_t n_neg = 0;
	size_t n_pos = 0;
	size_t counter_neg = 0;
	size_t counter_pos = 0;
	const size_t threshold_warning = 50;
	const double fmin_squared = fmin * fmin;
	
	// Determine number of positive and negative detections
	for(size_t i = self->size; i--;)
	{
		if(self->f_sum[i] < 0.0) ++n_neg;
		else if(self->f_sum[i] > 0.0) ++n_pos;
	}
	
	ensure(n_neg, ERR_FAILURE, "No negative sources found. Cannot proceed.");
	ensure(n_pos, ERR_FAILURE, "No positive sources found. Cannot proceed.");
	message("Found %zu positive and %zu negative sources.", n_pos, n_neg);
	if(n_neg < threshold_warning) warning("Only %zu negative %s found.\n         Reliability calculation may not be accurate.", n_neg, n_neg > 1 ? "detections" : "detection");
	
	// Extract relevant parameters
	// Arrays of parameters and indices for associated positive and negative detections
	double *par_pos = (double *)memory(MALLOC, dim * n_pos, sizeof(double));
	size_t *idx_pos = (size_t *)memory(MALLOC, n_pos, sizeof(size_t));
	double *par_neg = (double *)memory(MALLOC, dim * n_neg, sizeof(double));
	size_t *idx_neg = (size_t *)memory(MALLOC, n_neg, sizeof(size_t));
	
	// Loop over all detections
	for(size_t i = self->size; i--;)
	{
		if(self->f_sum[i] < 0.0)
		{
			// Negative detection
			bool include_source = true;
			
			if(rel_cat != NULL)
			{
				// Exclude positions from user-specified catalogue
				for(size_t row = 0; row < Table_rows(rel_cat); ++row)
				{
					const double cat_x = Table_get(rel_cat, row, 0);
					const double cat_y = Table_get(rel_cat, row, 1);
					
					if(cat_x >= (double)(self->x_min[i]) && cat_x <= (double)(self->x_max[i]) && cat_y >= (double)(self->y_min[i]) && cat_y <= (double)(self->y_max[i]))
					{
						include_source = false;
						break;
					}
				}
			}
			
			if(include_source)
			{
				ensure(self->f_min[i] < 0.0, ERR_FAILURE, "Non-negative minimum assigned to source with negative flux!");
				
				for(int j = 0; j < dim; ++j)
				{
					switch(Array_siz_get(rel_par_space, j))
					{
						case LINKERPAR_PEAK:
							// peak
							par_neg[dim * counter_neg + j] = log10(-self->f_min[i]);
							break;
						case LINKERPAR_SUM:
							// sum
							par_neg[dim * counter_neg + j] = log10(-self->f_sum[i]);
							break;
						case LINKERPAR_MEAN:
							// mean
							par_neg[dim * counter_neg + j] = log10(-self->f_sum[i] / self->n_pix[i]);
							break;
						case LINKERPAR_CHAN:
							// n_chan
							par_neg[dim * counter_neg + j] = self->z_max[i] - self->z_min[i] + 1;
							break;
						case LINKERPAR_PIX:
							// n_pix
							par_neg[dim * counter_neg + j] = log10(self->n_pix[i]);
							break;
						case LINKERPAR_FILL:
							// filling factor
							par_neg[dim * counter_neg + j] = log10(self->fill[i]);
							break;
						case LINKERPAR_STD:
							// std. dev.
							par_neg[dim * counter_neg + j] = sqrt(self->m2[i] / (self->n_pix[i] - 1));
							break;
						case LINKERPAR_SKEW:
							// skewness
							par_neg[dim * counter_neg + j] = sqrt(self->n_pix[i]) * self->m3[i] / pow(self->m2[i], 1.5);
							break;
						case LINKERPAR_KURT:
							// kurtosis
							par_neg[dim * counter_neg + j] = self->m4[i] * self->n_pix[i] / (self->m2[i] * self->m2[i]) - 3.0;
							break;
						default:
							ensure(false, ERR_USER_INPUT, "Unknown parameter requested for reliability measurement (%d).", Array_siz_get(rel_par_space, j));
					}
				}
				
				idx_neg[counter_neg] = i;
				++counter_neg;
			}
		}
		else if(self->f_sum[i] > 0.0)
		{
			// Positive detection
			ensure(self->f_max[i] > 0.0, ERR_FAILURE, "Non-positive maximum assigned to source with positive flux!");
			
			for(int j = 0; j < dim; ++j)
			{
				switch(Array_siz_get(rel_par_space, j))
				{
					case LINKERPAR_PEAK:
						// peak
						par_pos[dim * counter_pos + j] = log10(self->f_max[i]);
						break;
					case LINKERPAR_SUM:
						// sum
						par_pos[dim * counter_pos + j] = log10(self->f_sum[i]);
						break;
					case LINKERPAR_MEAN:
						// mean
						par_pos[dim * counter_pos + j] = log10(self->f_sum[i] / self->n_pix[i]);
						break;
					case LINKERPAR_CHAN:
						// n_chan
						par_pos[dim * counter_pos + j] = self->z_max[i] - self->z_min[i] + 1;
						break;
					case LINKERPAR_PIX:
						// n_pix
						par_pos[dim * counter_pos + j] = log10(self->n_pix[i]);
						break;
					case LINKERPAR_FILL:
						// filling factor
						par_pos[dim * counter_pos + j] = log10(self->fill[i]);
						break;
					case LINKERPAR_STD:
						// std. dev.
						par_pos[dim * counter_pos + j] = sqrt(self->m2[i] / (self->n_pix[i] - 1));
						break;
					case LINKERPAR_SKEW:
						// skewness
						par_pos[dim * counter_pos + j] = sqrt(self->n_pix[i]) * self->m3[i] / pow(self->m2[i], 1.5);
						break;
					case LINKERPAR_KURT:
						// skewness
						par_pos[dim * counter_pos + j] = self->m4[i] * self->n_pix[i] / (self->m2[i] * self->m2[i]) - 3.0;
						break;
					default:
						ensure(false, ERR_USER_INPUT, "Unknown parameter requested for reliability measurement (%d).", Array_siz_get(rel_par_space, j));
				}
			}
			
			idx_pos[counter_pos] = i;
			++counter_pos;
		}
	}
	
	// Adjust array sizes if necessary (some negative sources may have been removed)
	if(counter_neg < n_neg)
	{
		message("Excluding %zu out of %zu negative sources from reliability analysis.", n_neg - counter_neg, n_neg);
		n_neg = counter_neg;
		ensure(n_neg, ERR_FAILURE, "No negative sources found. Cannot proceed.");
		if(n_neg < threshold_warning) warning("Only %zu negative detections found.\n         Reliability calculation may not be accurate.", n_neg);
		par_neg = (double *)memory_realloc(par_neg, dim * n_neg, sizeof(double));
		idx_neg = (size_t *)memory_realloc(idx_neg, n_neg, sizeof(size_t));
	}
	else message("Retaining all negative detections.");
	
	// Determine covariance matrix from negative detections
	Matrix *covar = Matrix_covar(dim, n_neg, par_neg);
	Matrix *covar_inv = NULL;
	
	// Inverse of the square root of |2 * pi * covar| = (2 pi)^n |covar|
	// This is the scale factor needed to calculate the PDF of the multivariate normal distribution later on.
	//const double scal_fact = 1.0 / sqrt(Matrix_det(covar, 2.0 * M_PI));
	const double scal_fact = 1.0;
	// NOTE: This can be set to 1, as we don’t really care about the correct
	//       normalisation of the Gaussian kernel, so we might as well normalise
	//       the amplitude to 1 rather than the integral. The normalisation factor 
	//       does matter for the Skellam parameter, though.
	
	// Check if auto-kernel enabled
	if(!autokernel)
	{
		// No -> use fixed kernel scale factor
		Matrix_mul_scalar(covar, pow(*scale_kernel, 2));  // NOTE: Variance = sigma^2, hence scale_kernel^2 here.
		
		// Invert covariance matrix
		covar_inv = Matrix_invert(covar);
		ensure(covar_inv != NULL, ERR_FAILURE, "Covariance matrix is not invertible; cannot measure reliability.\n       Ensure that there are enough negative detections.");
		
		// Create Skellam array if requested
		if(skellam != NULL) LinkerPar_calculate_skellam(skellam, covar_inv, par_pos, par_neg, dim, n_pos, n_neg, scal_fact);
	}
	else
	{
		// Yes -> run auto-kernel by incrementally increasing scale
		message("Using auto-kernel feature.");
		
		int iter = 0;
		double scale = 0.1;
		double scale_old = 1.0;  // Must be initialised with 1!
		double skellam_med = 1e+5;
		const double scale_default = *scale_kernel;
		const double step_size = 0.02;
		
		while(iter < iterations && skellam_med > tolerance)
		{
			// Calculate skellam array
			Matrix_mul_scalar(covar, pow(scale / scale_old, 2));  // NOTE: Variance = sigma^2, hence scale_kernel^2 here.
			covar_inv = Matrix_invert(covar);
			ensure(covar_inv != NULL, ERR_FAILURE, "Covariance matrix is not invertible; cannot measure reliability.\n       Ensure that there are enough negative detections.");
			LinkerPar_calculate_skellam(skellam, covar_inv, par_pos, par_neg, dim, n_pos, n_neg, scal_fact);
			
			// Calculate new median
			skellam_med = fabs(median_dbl((double *)Array_dbl_get_ptr(*skellam), Array_dbl_get_size(*skellam), false));
			// NOTE: Casting constness away, as median needs to partially sort array.
			
			// Update with larger scale change far from target
			scale_old = scale;
			if(skellam_med < 10.0 * tolerance) scale += 2.0 * step_size;
			else if(skellam_med < 30.0 * tolerance) scale +=  5.0 * step_size;
			else if(skellam_med < 50.0 * tolerance) scale += 10.0 * step_size;
			else scale += step_size;
			
			++iter;
			message("  Iter. %*d: kernel = %.3f, median = %.3f", 2, iter, scale_old, skellam_med);
		}
		
		// Check if algorithm converged
		if(skellam_med <= tolerance)
		{
			*scale_kernel = scale_old;
			message("Converged to scale_kernel = %.3f after %d iterations.", scale_old, iter);
		}
		else
		{
			*scale_kernel = scale_default;
			
			Matrix_mul_scalar(covar, pow(*scale_kernel / scale_old, 2));
			covar_inv = Matrix_invert(covar);
			ensure(covar_inv != NULL, ERR_FAILURE, "Covariance matrix is not invertible; cannot measure reliability.\n       Ensure that there are enough negative detections.");
			if(skellam != NULL) LinkerPar_calculate_skellam(skellam, covar_inv, par_pos, par_neg, dim, n_pos, n_neg, scal_fact);
			warning("Auto-kernel failed to converge, defaulting to kernel scale of %.3f.", *scale_kernel);
		}
	}
	
	// Loop over all positive detections to measure their reliability
	const size_t cadence = (n_pos / 100) ? n_pos / 100 : 1;  // Only needed for progress bar
	size_t progress = 0;
	message("");
	
	#pragma omp parallel
	{
		Matrix *vector = Matrix_new(dim, 1);
		
		// Loop over all positive detections to calculate their reliability
		#pragma omp for schedule(static)
		for(size_t i = 0; i < n_pos; ++i)
		{
			#pragma omp critical
			if(++progress % cadence == 0 || progress == n_pos) progress_bar("Progress: ", progress, n_pos);
			
			// Only process sources above fmin and minpix
			if(self->f_sum[idx_pos[i]] * self->f_sum[idx_pos[i]] / self->n_pix[idx_pos[i]] > fmin_squared && self->n_pix[idx_pos[i]] > minpix)
			{
				// Multivariate kernel density estimation for negative detections
				double pdf_neg_sum = 0.0;
				
				for(double *ptr = par_neg; ptr < par_neg + n_neg * dim;)
				{
					// Set up relative position vector
					for(int j = 0; j < dim; ++j)
					{
						Matrix_set_value_nocheck(vector, j, 0, *ptr - par_pos[dim * i + j]);
						++ptr;
					}
					pdf_neg_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
				}
				
				// Multivariate kernel density estimation for positive detections
				double pdf_pos_sum = 0.0;
				
				for(double *ptr = par_pos; ptr < par_pos + n_pos * dim;)
				{
					// Set up relative position vector
					for(int j = 0; j < dim; ++j)
					{
						Matrix_set_value_nocheck(vector, j, 0, *ptr - par_pos[dim * i + j]);
						++ptr;
					}
					pdf_pos_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scal_fact);
				}
				
				// Determine reliability
				self->rel[idx_pos[i]] = pdf_pos_sum > pdf_neg_sum ? (pdf_pos_sum - pdf_neg_sum) / pdf_pos_sum : 0.0;
			}
		}
		
		Matrix_delete(vector);
	}
	
	// Release memory again
	Matrix_delete(covar_inv);
	free(par_pos);
	free(par_neg);
	free(idx_pos);
	free(idx_neg);
	
	return covar;
}


/// @brief Create reliability diagnostic plots
///
/// Public method for creating a diagnostic plot for the reliability measurement.
/// The method will generate an Encapsulated PostScript (EPS) file that shows the
/// distribution of negative and positive sources in 2-D projections of the parameter
/// space and highlight the ones that do or don't fulfil the reliability threshold or
/// `fmin` requirements. In addition, the location of `fmin` will be plotted, and
/// the Gaussian smoothing kernel used in the reliability measurement will be
/// overplotted as an ellipse based on the provided covariance matrix.
///
/// Note that if overwrite is set the false and the output file already exists, the
/// process will be terminated.
///
/// @param self           Object self-reference.
/// @param rel_par_space  Array of parameters used in determining the reliability of
///                       detections. These must be integer values corresponding to the
///                       parameters defined in the header file.
/// @param threshold      Reliability threshold.
/// @param fmin           Threshold for SNR filtering.
/// @param minSNR         Only needed for labelling plot, while `fmin` is the parameter
///                       used for drawing SNR line.
/// @param covar          Covariance matrix.
/// @param filename       Name of the output EPS file.
/// @param overwrite      If `true`, overwrite output file, otherwise do not overwrite.

PUBLIC void LinkerPar_rel_plots(const LinkerPar *self, const Array_siz *rel_par_space, const double threshold, const double fmin, const double minSNR, const Matrix *covar, const char *filename, const bool overwrite)
{
	// Sanity checks
	check_null(self);
	if(self->size == 0)
	{
		warning("No sources found; cannot generate reliability plots.");
		return;
	}
	ensure(filename != NULL && strlen(filename), ERR_USER_INPUT, "Empty file name for reliability plot provided.");
	
	// Dimensionality of parameter space
	const int dim = Array_siz_get_size(rel_par_space);
	const char *axis_labels[] = {"log\\(peak / rms\\)", "log\\(sum / rms\\)", "log\\(mean / rms\\)", "n_chan", "log\\(n_pix\\)", "log\\(fill_fact\\)", "std. dev.", "skewness", "kurtosis"};
	
	// Some settings
	const size_t plot_size_x = 300;  // pt
	const size_t plot_size_y = 300;  // pt
	
	const char *colour_neg = "1 0.4 0.4";
	const char *colour_pos = "0.4 0.4 1";
	const char *colour_rel = "0 0 0";
	const char *colour_kernel = "0.8 0.8 0.8";
	const char *colour_fmin = "0.5 0.5 0.5";
	const char *colour_axes = "0 0 0";
	
	// Create arrays for parameters
	double *data_x = (double *)memory(MALLOC, self->size, sizeof(double));
	double *data_y = (double *)memory(MALLOC, self->size, sizeof(double));
	
	// Open PS file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", filename);
	
	message("Creating postscript file: %s", strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/') + 1);
	
	// Print PS header
	String *bbox = String_new("0 ");
	String_append_int(bbox, "%d ", plot_size_y + 60);
	String_append_int(bbox, "%d ", 10 + (dim - 1) * (plot_size_x + 50));
	String_append_int(bbox, "%d", 10 + dim * (plot_size_y + 50));
	write_eps_header(fp, "SoFiA Reliability Plots", SOFIA_VERSION_FULL, String_get(bbox));
	String_delete(bbox);
	
	// Select font
	fprintf(fp, "roman\n");
	
	for(int p1 = 0; p1 < dim - 1; ++p1)
	{
		for(int p2 = p1 + 1; p2 < dim; ++p2)
		{
			// Read values and determine plotting range
			size_t plot_offset_x = 50 + p1 * (plot_size_x + 50);
			size_t plot_offset_y = 50 + p2 * (plot_size_y + 50);
			
			double data_min_x =  1e+30;
			double data_max_x = -1e+30;
			double data_min_y =  1e+30;
			double data_max_y = -1e+30;
			
			double radius_maj, radius_min, pa;
			Matrix_err_ellipse(covar, p1, p2, &radius_maj, &radius_min, &pa);
			
			for(size_t i = 0; i < self->size; ++i)
			{
				// Extract relevant parameters
				// x-axis
				switch(Array_siz_get(rel_par_space, p1))
				{
					case LINKERPAR_PEAK:
						if(self->f_sum[i] < 0.0) data_x[i] = log10(-self->f_min[i]);
						else data_x[i] = log10(self->f_max[i]);
						break;
					case LINKERPAR_SUM:
						if(self->f_sum[i] < 0.0) data_x[i] = log10(-self->f_sum[i]);
						else data_x[i] = log10(self->f_sum[i]);
						break;
					case LINKERPAR_MEAN:
						if(self->f_sum[i] < 0.0) data_x[i] = log10(-self->f_sum[i] / self->n_pix[i]);
						else data_x[i] = log10(self->f_sum[i] / self->n_pix[i]);
						break;
					case LINKERPAR_CHAN:
						data_x[i] = self->z_max[i] - self->z_min[i] + 1;
						break;
					case LINKERPAR_PIX:
						data_x[i] = log10(self->n_pix[i]);
						break;
					case LINKERPAR_FILL:
						data_x[i] = log10(self->fill[i]);
						break;
					case LINKERPAR_STD:
						data_x[i] = sqrt(self->m2[i] / (self->n_pix[i] - 1));
						break;
					case LINKERPAR_SKEW:
						data_x[i] = sqrt(self->n_pix[i]) * self->m3[i] / pow(self->m2[i], 1.5);
						break;
					case LINKERPAR_KURT:
						data_x[i] = self->m4[i] * self->n_pix[i] / (self->m2[i] * self->m2[i]) - 3.0;
						break;
					default:
						ensure(false, ERR_USER_INPUT, "Unknown parameter requested in reliability plot (%d).", Array_siz_get(rel_par_space, p1));
				}
				
				// y-axis
				switch(Array_siz_get(rel_par_space, p2))
				{
					case LINKERPAR_PEAK:
						if(self->f_sum[i] < 0.0) data_y[i] = log10(-self->f_min[i]);
						else data_y[i] = log10(self->f_max[i]);
						break;
					case LINKERPAR_SUM:
						if(self->f_sum[i] < 0.0) data_y[i] = log10(-self->f_sum[i]);
						else data_y[i] = log10(self->f_sum[i]);
						break;
					case LINKERPAR_MEAN:
						if(self->f_sum[i] < 0.0) data_y[i] = log10(-self->f_sum[i] / self->n_pix[i]);
						else data_y[i] = log10(self->f_sum[i] / self->n_pix[i]);
						break;
					case LINKERPAR_CHAN:
						data_y[i] = self->z_max[i] - self->z_min[i] + 1;
						break;
					case LINKERPAR_PIX:
						data_y[i] = log10(self->n_pix[i]);
						break;
					case LINKERPAR_FILL:
						data_y[i] = log10(self->fill[i]);
						break;
					case LINKERPAR_STD:
						data_y[i] = sqrt(self->m2[i] / (self->n_pix[i] - 1));
						break;
					case LINKERPAR_SKEW:
						data_y[i] = sqrt(self->n_pix[i]) * self->m3[i] / pow(self->m2[i], 1.5);
						break;
					case LINKERPAR_KURT:
						data_y[i] = self->m4[i] * self->n_pix[i] / (self->m2[i] * self->m2[i]) - 3.0;
						break;
					default:
						ensure(false, ERR_USER_INPUT, "Unknown parameter requested in reliability plot (%d).", Array_siz_get(rel_par_space, p2));
				}
				
				if(data_min_x > data_x[i]) data_min_x = data_x[i];
				if(data_max_x < data_x[i]) data_max_x = data_x[i];
				if(data_min_y > data_y[i]) data_min_y = data_y[i];
				if(data_max_y < data_y[i]) data_max_y = data_y[i];
			}
			
			double data_range_x = data_max_x - data_min_x;
			double data_range_y = data_max_y - data_min_y;
			
			// Add a little bit of margin
			data_min_x -= 0.05 * data_range_x;
			data_max_x += 0.05 * data_range_x;
			data_min_y -= 0.05 * data_range_y;
			data_max_y += 0.05 * data_range_y;
			
			// Determine optimal tick mark increments
			const double tick_inc_x = auto_tick(data_max_x - data_min_x, 4);
			const double tick_inc_y = auto_tick(data_max_y - data_min_y, 4);
			
			// Determine the mean of negative sources
			double mean_x = 0.0;
			double mean_y = 0.0;
			size_t counter = 0;
			
			for(size_t i = 0; i < self->size; ++i)
			{
				if(self->f_sum[i] < 0.0)
				{
					mean_x += data_x[i];
					mean_y += data_y[i];
					++counter;
				}
			}
			
			mean_x /= counter;
			mean_y /= counter;
			
			const double centre_x = (mean_x - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
			const double centre_y = (mean_y - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
			const double radius_x = radius_maj * plot_size_x / (data_max_x - data_min_x);
			const double radius_y = radius_min * plot_size_x / (data_max_x - data_min_x);
			
			// Determine scale factor of kernel ellipse
			const double scale_factor = data_range_x / data_range_y;
			
			// Plot negative sources
			fprintf(fp, "%s rgb\n", colour_neg);
			fprintf(fp, "0.5 lw\n");
			fprintf(fp, "np\n");
			
			for(size_t i = self->size; i--;)
			{
				if(self->f_sum[i] < 0.0)
				{
					const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
					const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
					
					fprintf(fp, "%.1f %.1f 1 0 360 af\n", plot_x, plot_y);
				}
			}
			
			// Plot unreliable positive sources
			fprintf(fp, "%s rgb\n", colour_pos);
			
			for(size_t i = self->size; i--;)
			{
				if(self->f_sum[i] > 0.0 && self->rel[i] < threshold)
				{
					const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
					const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
					
					fprintf(fp, "%.1f %.1f 1 0 360 af\n", plot_x, plot_y);
				}
			}
			
			// Plot reliable positive sources
			fprintf(fp, "%s rgb\n", colour_rel);
			
			for(size_t i = self->size; i--;)
			{
				if(self->f_sum[i] > 0.0 && self->rel[i] >= threshold)
				{
					const double plot_x = (data_x[i] - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
					const double plot_y = (data_y[i] - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
					
					if(self->f_sum[i] / sqrt(self->n_pix[i]) > fmin) fprintf(fp, "%.1f %.1f 2 0 360 af\n", plot_x, plot_y);
					else fprintf(fp, "%.1f %.1f 2 0 360 as\n", plot_x, plot_y);
				}
			}
			
			// Plot kernel ellipse
			fprintf(fp, "gsave\n");
			fprintf(fp, "%s rgb\n", colour_kernel);
			fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
			fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, radius_x, radius_y, 180.0 * pa / M_PI, scale_factor);
			fprintf(fp, "[2 2] 0 setdash\n");
			fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, 2.0 * radius_x, 2.0 * radius_y, 180.0 * pa / M_PI, scale_factor);
			fprintf(fp, "[0.5 1.5] 0 setdash\n");
			fprintf(fp, "%.2f %.2f %.2f %.2f %.2f %.2f ellipse\n", centre_x, centre_y, 3.0 * radius_x, 3.0 * radius_y, 180.0 * pa / M_PI, scale_factor);
			fprintf(fp, "grestore\n");
			
			// Plot fmin and npix lines if possible
			if(Array_siz_get(rel_par_space, p1) == LINKERPAR_SUM && Array_siz_get(rel_par_space, p2) == LINKERPAR_MEAN)
			{
				fprintf(fp, "gsave\n");
				
				// fmin
				double plot_x = plot_offset_x;
				double plot_y = (2.0 * log10(fmin) - data_min_x - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				
				fprintf(fp, "%s rgb\n", colour_fmin);
				fprintf(fp, "[3 3] 0 setdash\n");
				fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				
				plot_x = plot_offset_x + plot_size_x;
				plot_y = (2.0 * log10(fmin) - data_max_x - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				
				fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
				
				fprintf(fp, "np %zu %.zu m\n", plot_offset_x + 14, plot_offset_y + plot_size_y - 20);
				fprintf(fp, "%s rgb\n", colour_fmin);
				fprintf(fp, "(minSNR = %.1f) show\n", minSNR);
				
				// npix
				// 10
				fprintf(fp, "0.9 0.9 0.9 rgb\n");
				fprintf(fp, "[0.5 2] 0 setdash\n");
				fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
				
				plot_x = plot_offset_x;
				plot_y = (data_min_x - 1 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				
				plot_x = plot_offset_x + plot_size_x;
				plot_y = (data_max_x - 1 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
				
				plot_x = (data_min_y + 1 - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x + 10;
				plot_y = plot_offset_y + 10;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				fprintf(fp, "(10) show\n");
				
				// 100
				fprintf(fp, "0.8 0.8 0.8 rgb\n");
				fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
				
				plot_x = plot_offset_x;
				plot_y = (data_min_x - 2 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				
				plot_x = plot_offset_x + plot_size_x;
				plot_y = (data_max_x - 2 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
				
				plot_x = (data_min_y + 2 - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x + 10;
				plot_y = plot_offset_y + 10;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				fprintf(fp, "(100) show\n");
				
				// 1000
				fprintf(fp, "0.7 0.7 0.7 rgb\n");
				fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
				
				plot_x = plot_offset_x;
				plot_y = (data_min_x - 3 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				
				plot_x = plot_offset_x + plot_size_x;
				plot_y = (data_max_x - 3 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
				
				plot_x = (data_min_y + 3 - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x + 10;
				plot_y = plot_offset_y + 10;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				fprintf(fp, "(1000) show\n");
				
				// 10000
				fprintf(fp, "0.6 0.6 0.6 rgb\n");
				fprintf(fp, "np %zu %zu m %zu %zu l %zu %zu l %zu %zu l cp clip\n", plot_offset_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y, plot_offset_x + plot_size_x, plot_offset_y + plot_size_y, plot_offset_x, plot_offset_y + plot_size_y);
				
				plot_x = plot_offset_x;
				plot_y = (data_min_x - 4 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				
				plot_x = plot_offset_x + plot_size_x;
				plot_y = (data_max_x - 4 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "%.2f %.2f l s\n", plot_x, plot_y);
				
				plot_x = (data_min_y + 4 - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x + 10;
				plot_y = plot_offset_y + 10;
				fprintf(fp, "%.2f %.2f m\n", plot_x, plot_y);
				fprintf(fp, "(10000 px) show\n");
				
				fprintf(fp, "grestore\n");
			}
			
			// Plot frame
			fprintf(fp, "%s rgb\n", colour_axes);
			fprintf(fp, "[] 0 setdash\n");
			fprintf(fp, "np\n");
			fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
			fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
			fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
			fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
			fprintf(fp, "cp s\n");
			
			// Plot tick marks
			for(double tm = ceil(data_min_x / tick_inc_x) * tick_inc_x; tm <= data_max_x; tm += tick_inc_x)
			{
				if(fabs(tm) < 0.001) tm = 0.0;
				const double plot_x = (tm - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
				fprintf(fp, "np %.2f %zu m %.2f %zu l s\n", plot_x, plot_offset_y, plot_x, plot_offset_y + 5);
				fprintf(fp, "np %.2f %zu m (%.1f) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_x, plot_offset_y - 14, tm);
			}
			
			for(double tm = ceil(data_min_y / tick_inc_y) * tick_inc_y; tm <= data_max_y; tm += tick_inc_y)
			{
				if(fabs(tm) < 0.001) tm = 0.0;
				const double plot_y = (tm - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
				fprintf(fp, "np %zu %.2f m %zu %.2f l s\n", plot_offset_x, plot_y, plot_offset_x + 5, plot_y);
				fprintf(fp, "np %zu %.2f m (%.1f) dup stringwidth pop neg 0 rmoveto show\n", plot_offset_x - 4, plot_y - 4.0, tm);
			}
			
			// Print labels
			fprintf(fp, "np %zu %zu m (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_offset_x + plot_size_x / 2, plot_offset_y - 30, axis_labels[Array_siz_get(rel_par_space, p1)]);
			
			fprintf(fp, "np %zu %zu m gsave 90 rotate (%s) dup stringwidth pop 2 div neg 0 rmoveto show grestore\n", plot_offset_x - 34, plot_offset_y + plot_size_y / 2, axis_labels[Array_siz_get(rel_par_space, p2)]);
		}
	}
	
	// Print EPS footer
	write_eps_footer(fp);
	
	// Close output file
	fclose(fp);
	
	// Clean up
	free(data_x);
	free(data_y);
	
	return;
}


/// @brief Create Skellam diagnostic plot
///
/// Public function for generating a Skellam diagnostic plot showing the cumulative
/// distribution of values in the Array called `skellam` which must contain values
/// of (P - N) / sqrt(P + N) generated by the reliability module.The resulting plot
/// will be written to an EPS file with the specified file name.
///
/// @param skellam      Array of (P - N) / sqrt(P + N) values.
/// @param filename     Output file name for plot.
/// @param overwrite    If true, overwrite existing file.
/// @param kernelScale  Kernel scale factor, for labelling only.

PUBLIC void LinkerPar_skellam_plot(Array_dbl *skellam, const char *filename, const bool overwrite, const double kernelScale)
{
	// Sanity checks
	check_null(skellam);
	const size_t size = Array_dbl_get_size(skellam);
	ensure(size, ERR_USER_INPUT, "Failed to create Skellam plot; no valid data found.");
	
	// Sort the Skellam array
	Array_dbl_sort(skellam);
	
	// Scale all Skellam values by standard deviation
	// NOTE: This is necessary to ensure that the standard deviation of the
	//       Skellam parameter values is 1 such that their distribution can
	//       be readily compared to a standard Gaussian.
	//const double skel_mean = mean_dbl(Array_dbl_get_ptr(skellam), size);
	// NOTE: Using median instead of mean, as it is more robust
	const double skel_mean = IS_ODD(size) ? Array_dbl_get(skellam, size / 2) : 0.5 * (Array_dbl_get(skellam, size / 2 - 1) + Array_dbl_get(skellam, size / 2));
	const double skel_std  = std_dev_val_dbl(Array_dbl_get_ptr(skellam), size, skel_mean, 1, 0);
	for(size_t i = 0; i < size; ++i) Array_dbl_set(skellam, i, (Array_dbl_get(skellam, i) - skel_mean) / skel_std + skel_mean);
	
	// Determine plotting range
	//double data_min_x = Array_dbl_get(skellam, 0);
	//double data_max_x = Array_dbl_get(skellam, size - 1);
	const double data_min_x = -4.0;
	const double data_max_x = 4.0;
	const double data_min_y = 0.0;
	const double data_max_y = 1.0;
	
	// Plot geometry
	const size_t plot_size_x = 410;
	const size_t plot_size_y = 310;
	const size_t plot_offset_x = 50;
	const size_t plot_offset_y = 40;
	
	// Determine optimal tick mark increments
	const double tick_inc_x = auto_tick(data_max_x - data_min_x, 5);
	const double tick_inc_y = auto_tick(data_max_y - data_min_y, 5);
	
	// Colours
	const char *colour_data = "1 0 0";
	const char *colour_erf  = "0.3 0.3 0.3";
	const char *colour_zero = "0.7 0.7 0.7";
	const char *colour_axes = "0 0 0";
	
	// Labels
	const char *label_x_placeholder = "\\(P - N\\) / sqrt\\(P + N\\) normalised to s = 1";
	const char *label_x = "\\(P - N\\) / sqrt\\(P + N\\) normalised to ";
	const char *label_y = "Cumulative fraction";
	
	// Open output file
	FILE *fp;
	if(overwrite) fp = fopen(filename, "wb");
	else fp = fopen(filename, "wxb");
	ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", filename);
	
	message("Creating postscript file: %s", strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/') + 1);
	
	// Print PS header
	write_eps_header(fp, "SoFiA Skellam Plot", SOFIA_VERSION_FULL, "0 0 480 360");
	
	// Select font
	fprintf(fp, "roman\n");
	
	// Set clip path
	fprintf(fp, "gsave\n");
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
	fprintf(fp, "cp clip\n");
	
	// Plot vertical line at 0
	fprintf(fp, "np\n");
	fprintf(fp, "%.1f %.1f m\n", -data_min_x * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, -data_min_y * plot_size_y / (data_max_y - data_min_y) + plot_offset_y);
	fprintf(fp, "%.1f %.1f l\n", -data_min_x * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, (1.0 - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y);
	fprintf(fp, "%s rgb\n", colour_zero);
	fprintf(fp, "[5 3] 0 setdash\n");
	fprintf(fp, "s\n");
	
	// Plot vertical line at mean
	fprintf(fp, "np\n");
	fprintf(fp, "%.2f %zu m\n", (skel_mean - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, plot_offset_y);
	fprintf(fp, "%.2f %zu l\n", (skel_mean - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x, plot_offset_y + plot_size_y);
	fprintf(fp, "%s rgb\n", colour_data);
	fprintf(fp, "[1 3] 0 setdash\n");
	fprintf(fp, "s\n");
	
	// Plot error function for Gaussian with sigma = 1
	fprintf(fp, "np\n");
	for(int i = -100; i <= 100; ++i)
	{
		const double x = 0.05 * i;
		const double plot_x = (x - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		const double plot_y = (0.5 - 0.5 * erf(-x / sqrt(2.0)) - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		
		fprintf(fp, "%.2f %.2f %s\n", plot_x, plot_y, i == -100 ? "m" : "l");
	}
	fprintf(fp, "%s rgb\n", colour_erf);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	
	// Plot data
	fprintf(fp, "np\n");
	for(size_t i = 0; i < size; ++i)
	{
		const double plot_x = (Array_dbl_get(skellam, i) - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		const double plot_y = (((double)(i) / (double)(size - 1)) - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		
		fprintf(fp, "%.1f %.1f %s\n", plot_x, plot_y, i ? "l" : "m");
	}
	fprintf(fp, "%s rgb\n", colour_data);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "grestore\n");
	
	// Plot frame
	fprintf(fp, "%s rgb\n", colour_axes);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + plot_size_x, plot_offset_y + plot_size_y);
	fprintf(fp, "%zu %zu l\n", plot_offset_x, plot_offset_y + plot_size_y);
	fprintf(fp, "cp s\n");
	
	// Plot tick marks
	for(double tm = ceil(data_min_x / tick_inc_x) * tick_inc_x; tm <= data_max_x; tm += tick_inc_x)
	{
		if(fabs(tm) < 0.001) tm = 0.0;
		const double plot_x = (tm - data_min_x) * plot_size_x / (data_max_x - data_min_x) + plot_offset_x;
		fprintf(fp, "np %.2f %zu m %.2f %zu l s\n", plot_x, plot_offset_y, plot_x, plot_offset_y + 5);
		fprintf(fp, "np %.2f %zu m (%.1f) dup stringwidth pop 2 div neg 0 rmoveto show\n", plot_x, plot_offset_y - 14, tm);
	}
	
	for(double tm = ceil(data_min_y / tick_inc_y) * tick_inc_y; tm <= data_max_y; tm += tick_inc_y)
	{
		if(fabs(tm) < 0.001) tm = 0.0;
		const double plot_y = (tm - data_min_y) * plot_size_y / (data_max_y - data_min_y) + plot_offset_y;
		fprintf(fp, "np %zu %.2f m %zu %.2f l s\n", plot_offset_x, plot_y, plot_offset_x + 5, plot_y);
		fprintf(fp, "np %zu %.2f m (%.1f) dup stringwidth pop neg 0 rmoveto show\n", plot_offset_x - 4, plot_y - 4.0, tm);
	}
	
	// Print labels
	fprintf(fp, "np\n");
	fprintf(fp, "%zu 10 m (%s) stringwidth pop 2 div neg 0 rmoveto (%s) show\n", plot_offset_x + plot_size_x / 2, label_x_placeholder, label_x);
	fprintf(fp, "greek\n");
	fprintf(fp, "(s) show\n");
	fprintf(fp, "roman\n");
	fprintf(fp, "( = 1) show\n");
	fprintf(fp, "np %zu %zu m gsave 90 rotate (%s) dup stringwidth pop 2 div neg 0 rmoveto show grestore\n", plot_offset_x - 34, plot_offset_y + plot_size_y / 2, label_y);
	
	// Plot legend
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 20, plot_offset_y + plot_size_y - 20);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + 40, plot_offset_y + plot_size_y - 20);
	fprintf(fp, "%s rgb\n", colour_data);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 46, plot_offset_y + plot_size_y - 24);
	fprintf(fp, "(Data \\() show\n");
	fprintf(fp, "greek\n");
	fprintf(fp, "(m) show\n");
	fprintf(fp, "roman\n");
	fprintf(fp, "( = %.3f,) show\n", skel_mean);
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 46, plot_offset_y + plot_size_y - 39);
	fprintf(fp, "(kernel = %.2f\\)) show\n", kernelScale);
	
	fprintf(fp, "np\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 20, plot_offset_y + plot_size_y - 50);
	fprintf(fp, "%zu %zu l\n", plot_offset_x + 40, plot_offset_y + plot_size_y - 50);
	fprintf(fp, "%s rgb\n", colour_erf);
	fprintf(fp, "[] 0 setdash\n");
	fprintf(fp, "s\n");
	fprintf(fp, "%zu %zu m\n", plot_offset_x + 46, plot_offset_y + plot_size_y - 54);
	fprintf(fp, "(Gaussian \\() show\n");
	fprintf(fp, "greek\n");
	fprintf(fp, "(s) show\n");
	fprintf(fp, "roman\n");
	fprintf(fp, "( = 1\\)) show\n");
	
	// Print EPS footer
	write_eps_footer(fp);
	
	// Close output file
	fclose(fp);
	
	return;
}


/// @brief Calculate array of normalised Skellam values
///
/// Function for computing the Skellam array from positive and negative detections.
/// The Skellam value is given by
///
///   S = (P - N) / SQRT(P + N)
///
/// where P and N are kernel density estimations for positive and negative detections,
/// respectively. This function updates the provided `skellam` array with calculated
/// values.
///
/// The parameters for the positive and negative detections must be provided using the
/// `pos` and `neg` arrays which are expected to be flattened arrays of length
/// `dim * n_pos` and `dim * `n_neg`, respectively. Here, `dim` is the
/// number of parameters (i.e. dimensionality of the parameter space), while `n_pos`
/// and `n_neg` is the total number of positive and negative detections, respectively.
///
/// The size of the Gaussian density estimator kernel used in the calculation is
/// determined from the specified covariance matrix. An additional `scale` factor
/// can be specified to normalise the integral under the Gaussian kernel.
///
/// @param skellam    Pointer to skellam array.
/// @param covar_inv  Pointer to inverse of covariance matrix. Used to set the size of
///                   the Gaussian kernel. Can be calculated with Matrix_invert().
/// @param pos        Array of parameters for positive detections. Must be of length
///                   `dim * n_pos`.
/// @param neg        Array of parameters for negative detections. Must be of length
///                   `dim * n_neg`.
/// @param dim        Dimensionality of parameter space.
/// @param n_pos      Number of positive detections.
/// @param n_neg      Number of negative detections.
/// @param scale      Scale factor for normalisation.

PUBLIC void LinkerPar_calculate_skellam(Array_dbl **skellam, const Matrix *covar_inv, double *pos, double *neg, const int dim, const size_t n_pos, const size_t n_neg, const double scale)
{
	// Calculate skellam
	*skellam = Array_dbl_new(n_neg);
	
	#pragma omp parallel
	{
		Matrix *vector = Matrix_new(dim, 1);
		
		// Loop over all negative sources to derive Skellam distribution
		#pragma omp for schedule(static)
		for(size_t i = 0; i < n_neg; ++i)
		{
			// Multivariate kernel density estimation for negative detections
			double pdf_neg_sum = 0.0;
			
			for(double *ptr = neg; ptr < neg + n_neg * dim;)
			{
				// Set up relative position vector
				for(int j = 0; j < dim; ++j)
				{
					Matrix_set_value_nocheck(vector, j, 0, *ptr - neg[dim * i + j]);
					++ptr;
				}
				pdf_neg_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scale);
			}
			
			// Multivariate kernel density estimation for positive detections
			double pdf_pos_sum = 0.0;
			
			for(double *ptr = pos; ptr < pos + n_pos * dim;)
			{
				// Set up relative position vector
				for(int j = 0; j < dim; ++j)
				{
					Matrix_set_value_nocheck(vector, j, 0, *ptr - neg[dim * i + j]);
					++ptr;
				}
				pdf_pos_sum += Matrix_prob_dens_nocheck(covar_inv, vector, scale);
			}
			
			// Determine normalised Skellam parameter S = (P - N) / SQRT(P + N)
			Array_dbl_set(*skellam, i, (pdf_pos_sum - pdf_neg_sum) / sqrt(pdf_pos_sum + pdf_neg_sum));
		}
		
		Matrix_delete(vector);
	}
	return;
}
