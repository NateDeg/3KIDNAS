// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (statistics_flt.h) - Source Finding Application          //
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

/// @file   statistics_flt.h
/// @author Tobias Westmeier
/// @date   10/12/2021
/// @brief  Function templates for basic statistical analysis of data of type `float` (header).


// WARNING: This is a template that needs to be instantiated before use.
//          Do not edit template instances, as they are auto-generated
//          and will be overwritten during instantiation!


#ifndef STATISTICS_flt_H
#define STATISTICS_flt_H

#include <stdbool.h>
#include "common.h"

// -------------------------- //
// Settings for boxcar filter //
// -------------------------- //
#define BOXCAR_MIN_ITER 3  ///< Minimum number of iterations required for boxcar approximation of Gaussian smoothing kernel.
#define BOXCAR_MAX_ITER 6  ///< Maximum number of iterations allowed for boxcar approximation of Gaussian smoothing kernel.



// -------------------- //
// Statistics functions //
// -------------------- //

// Check if array contains NaN
bool contains_nan_flt(const float *data, const size_t size);
bool contains_inf_flt(float *data, const size_t size, const bool flag_inf);

// Maximum and minimum
void max_min_flt(const float *data, const size_t size, float *value_max, float *value_min);
float max_flt(const float *data, const size_t size);
float min_flt(const float *data, const size_t size);

// Sum and mean
double summation_flt(const float *data, const size_t size, const bool mean);
double sum_flt(const float *data, const size_t size);
double mean_flt(const float *data, const size_t size);

// N-th moment
double moment_flt(const float *data, const size_t size, unsigned int order, const double value);
void moments_flt(const float *data, const size_t size, const double value, double *m2, double *m3, double *m4);

// Standard deviation
double std_dev_flt(const float *data, const size_t size);
double std_dev_val_flt(const float *data, const size_t size, const double value, const size_t cadence, const int range);

// N-th-smallest element
float nth_element_flt(float *data, const size_t size, const size_t n);

// Median and MAD
float median_flt(float *data, const size_t size, const bool fast);
float median_safe_flt(const float *data, const size_t size, const bool fast);
float mad_flt(float *data, const size_t size);
float mad_val_flt(const float *data, const size_t size, const float value, const size_t cadence, const int range);

// Robust and fast noise measurement
float robust_noise_flt(const float *data, const size_t size);
float robust_noise_2_flt(const float *data, const size_t size);
float robust_noise_in_region_flt(const float *data, const size_t nx, const size_t ny, const size_t x1, const size_t x2, const size_t y1, const size_t y2, const size_t z1, const size_t z2);

// Gaussian fit to histogram
size_t *create_histogram_flt(const float *data, const size_t size, const size_t n_bins, const float data_min, const float data_max, const size_t cadence);
float gaufit_flt(const float *data, const size_t size, const size_t cadence, const int range);

// Skewness and kurtosis
void skew_kurt_flt(const float *data, const size_t size, double *skew, double *kurt);
double skewness_flt(const float *data, const size_t size);
double kurtosis_flt(const float *data, const size_t size);

// 1D boxcar filter
void filter_boxcar_1d_flt(float *data, float *data_copy, const size_t size, const size_t filter_radius);

// 2D Gaussian filter
void filter_gauss_2d_flt(float *data, float *data_copy, float *data_row, float *data_col, const size_t size_x, const size_t size_y, const size_t n_iter, const size_t filter_radius);

// Polynomial fitting
void shift_and_subtract_flt(float *data, const size_t size, const size_t shift);

// Source parameterisation
void   moment_ellipse_fit_flt  (const float *moment_map, const size_t *count_map, const size_t size_x, const size_t size_y, const float centroid_x, const float centroid_y, const float rms, float *ell_maj, float *ell_min, float *ell_pa, float *ell3s_maj, float *ell3s_min, float *ell3s_pa);
void   spectral_line_width_flt (const float *spectrum, const size_t size, float *w20, float *w50);
double wm50_line_width_flt     (const float *spectrum, const size_t size);
float kin_maj_axis_flt        (const float *centroidX, const float *centroidY, const float *sum, const size_t size, const size_t first, const size_t last);


// ----------------- //
// Utility functions //
// ----------------- //

// Filter size and iterations for Gaussian function
void optimal_filter_size_flt(const double sigma, size_t *filter_radius, size_t *n_iter);

#endif
