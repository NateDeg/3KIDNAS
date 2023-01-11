// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (statistics_flt.c) - Source Finding Application          //
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

/// @file   statistics_flt.c
/// @author Tobias Westmeier
/// @date   26/11/2021
/// @brief  Function templates for basic statistical analysis of data of type `float`.


// WARNING: This is a template that needs to be instantiated before use.
//          Do not edit template instances, as they are auto-generated
//          and will be overwritten during instantiation!


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "statistics_flt.h"



/// @brief Check if data contains `NaN`
///
/// Checks if the data array contains values of Not a
/// Number (`NaN`).
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return True if `NaN` found, false otherwise.

bool contains_nan_flt(const float *data, const size_t size)
{
	const float *ptr = data + size;
	while(ptr --> data) if(IS_NAN(*ptr)) return true;
	return false;
}



/// @brief Check if data contains `Inf`
///
/// Checks if the data array contains values of infinity
/// (`Inf`). If `flag_inf` is set to `true`, all values of
/// `Inf` will be replaced with `NaN`.
///
/// @param data      Pointer to the input data array.
/// @param size      Size of the input data array.
/// @param flag_inf  If true, set `Inf` values to `NaN`.
///
/// @return `true` if `Inf` found, `false` otherwise.

bool contains_inf_flt(float *data, const size_t size, const bool flag_inf)
{
	bool inf_found = false;
	float *ptr = data + size;
	
	while(ptr --> data)
	{
		if(isinf(*ptr))
		{
			if(!flag_inf) return true;
			inf_found = true;
			*ptr = NAN;
		}
	}
	
	return inf_found;
}



/// @brief Maximum and minimum
///
/// Determines the maximum and minimum value in the input
/// data array simultaneously. This is faster that calling
/// the max() and min() functions separately in situations
/// where both values are required. If the array contains
/// only `NaN`, `value_min` and `value_max` will both be
/// set to `NaN`.
///
/// @param data       Pointer to the input data array.
/// @param size       Size of the input data array.
/// @param value_max  Pointer to variable for holding maximum.
/// @param value_min  Pointer to variable for holding minimum.

void max_min_flt(const float *data, const size_t size, float *value_max, float *value_min)
{
	// Find the last non-NaN element
	const float *ptr = data + size - 1;
	while(IS_NAN(*ptr) && ptr > data) --ptr;
	
	// Copy values across
	*value_min = *ptr;
	*value_max = *ptr;
	
	// Check remaining elements
	for(const float *ptr2 = ptr; ptr2 --> data;)
	{
		if(*ptr2 < *value_min) *value_min = *ptr2;
		else if(*ptr2 > *value_max) *value_max = *ptr2;
	}
	
	return;
}



/// @brief Maximum
///
/// Determines the maximum value in the input data array.
/// If the array contains only `NaN`, then `NaN` is returned.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return Maximum value of the input array. `NaN` will be returned
/// for `NaN`-only input arrays.

float max_flt(const float *data, const size_t size)
{
	const float *result = data + size - 1;
	while(IS_NAN(*result) && data <-- result);
	const float *ptr = result;
	
	while(ptr --> data) if(*ptr > *result) result = ptr;
	
	return *result;
}



/// @brief Minimum
///
/// Determines the minimum value in the input data array.
/// If the array contains only `NaN`, then `NaN` is returned.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return Minimum value of the input array. `NaN` will be returned
/// for `NaN`-only input arrays.

float min_flt(const float *data, const size_t size)
{
	const float *result = data + size - 1;
	while(IS_NAN(*result) && data <-- result);
	const float *ptr = result;
	
	while(ptr --> data) if(*ptr < *result) result = ptr;
	
	return *result;
}



/// @brief Sum and mean
///
/// Depending on the value of mean, either the sum or the
/// mean of the non-`NaN` elements of the input array is
/// returned. Two convenient wrapper functions, sum() and
/// mean(), have been provided to explicitly return the sum
/// and the mean of the array, respectively.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
/// @param mean  If `true`, return the mean, otherwise the sum.
///
/// @return Sum or mean of the data array. `NaN` will be returned for
///         `NaN`-only input arrays.

double summation_flt(const float *data, const size_t size, const bool mean)
{
	double result = 0.0;
	size_t counter = 0;
	const float *ptr = data + size;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			result += *ptr;
			++counter;
		}
	}
	
	if(counter) return mean ? result / counter : result;
	return NAN;
}



/// @brief Sum
///
/// Wrapper function for summation_flt() that will return the
/// sum of the non-`NaN` elements of the input array.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return Sum of the data array values. `NaN` will be returned
///         for `NaN`-only input arrays.

double sum_flt(const float *data, const size_t size) { return summation_flt(data, size, false); }



/// @brief Mean
///
/// Wrapper function for summation_flt() that will return the
/// mean of the non-`NaN` elements of the input array.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return Mean of the data array values. `NaN` will be returned
///         for `NaN`-only input arrays.

double mean_flt(const float *data, const size_t size) { return summation_flt(data, size, true); }



/// @brief N-th moment about value
///
/// Calculates the N-th moment of the input data array. The
/// N-th moment is defined as
///
///   sum((x - value)^N) / n
///
/// where the summation is over all n data values, x, of
/// the input array that are not `NaN`. The most sensible
/// use would be to calculate the moment about the mean.
///
/// @param data   Pointer to the input data array.
/// @param size   Size of the input data array.
/// @param order  Order of the moment to be calculated.
/// @param value  Value about which to calculate the moment.
///
/// @return N-th moment of the data array. `NaN` will be returned for
///         `NaN`-only input arrays.

double moment_flt(const float *data, const size_t size, unsigned int order, const double value)
{
	if(order == 0) return 1.0;
	
	double result = 0.0;
	size_t counter = 0;
	double tmp;
	unsigned int i;
	const float *ptr = data + size;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			tmp = *ptr - value;
			for(i = order; --i;) tmp *= (*ptr - value);
			result += tmp;
			++counter;
		}
	}
	
	return counter ? result / counter : NAN;
}



/// @brief Moments 2, 3 and 4 about value
///
/// Calculates the 2nd, 3rd and 4th moment of a data array
/// simultaneously. This will enable faster calculation of
/// skewness and kurtosis, both of which depend on multiple
/// moments of the same data array. Moments are defined in
/// the same way as in function moment().
///
/// @param data   Pointer to the input data array.
/// @param size   Size of the input data array.
/// @param value  Value about which to calculate the moments.
/// @param m2     Pointer to a variable for storing moment 2.
/// @param m3     Pointer to a variable for storing moment 3.
/// @param m4     Pointer to a variable for storing moment 4.

void moments_flt(const float *data, const size_t size, const double value, double *m2, double *m3, double *m4)
{
	const float *ptr = data + size;
	double result2 = 0.0;
	double result3 = 0.0;
	double result4 = 0.0;
	size_t counter = 0;
	double tmp, tmp2;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			tmp = *ptr - value;
			tmp2 = tmp * tmp;
			result2 += tmp2;
			result3 += tmp2 * tmp;
			result4 += tmp2 * tmp2;
			++counter;
		}
	}
	
	if(counter)
	{
		*m2 = result2 / counter;
		*m3 = result3 / counter;
		*m4 = result4 / counter;
	}
	else
	{
		*m2 = NAN;
		*m3 = NAN;
		*m4 = NAN;
	}
	
	return;
}



/// @brief Standard deviation about value
///
/// Calculates the standard deviation about a user-specified
/// value of all elements in the input data array that are
/// not `NaN`. The standard deviation is defined as
///
///   sqrt[sum(x - value)^2 / n]
///
/// where the summation is over all n non-`NaN` elements, x,
/// of the input array.
///
/// @param data     Pointer to the input data array.
/// @param size     Size of the input data array.
/// @param value    Value about which to calculate the
///                 standard deviation.
/// @param cadence  Can be set to > 1 to speed up algorithm.
/// @param range    Flux range to be used. Can be Negative
///                 (-1), Full (0) or Positive (1).
///
/// @return Standard deviation about the specified value of
///         the elements in the input data array. `NaN` will
///         be returned for `NaN`-only input arrays.

double std_dev_val_flt(const float *data, const size_t size, const double value, const size_t cadence, const int range)
{
	const float *ptr = data + size;
	const float *ptr2 = data + cadence - 1;
	double result = 0.0;
	size_t counter = 0;
	double tmp;
	
	while(ptr > ptr2)
	{
		ptr -= cadence;
		
		if((range == 0 && IS_NOT_NAN(*ptr)) || (range < 0 && *ptr < 0.0) || (range > 0 && *ptr > 0.0))
		{
			tmp = (*ptr - value);
			result += tmp * tmp;
			++counter;
		}
	}
	
	return counter ? sqrt(result / counter) : NAN;
}



/// @brief Standard deviation about mean
///
/// Calculates the standard deviation about the mean of all
/// elements in the input data array that are not `NaN`. The
/// standard deviation is defined as
///
///   sqrt{sum[x - mean(x)]^2 / n}
///
/// where the summation is over all n non-`NaN` elements, x,
/// of the input array.
///
/// @param data  Pointer to the input data array.
/// @param size  Size of the input data array.
///
/// @return Standard deviation about the mean of the elements
///         in the input data array. `NaN` will be returned for
///         `NaN`-only input arrays.

double std_dev_flt(const float *data, const size_t size)
{
	return std_dev_val_flt(data, size, mean_flt(data, size), 1, 0);
}



/// @brief N-th smallest element in array
///
/// Partially sorts the input data array until the
/// n-th-smallest element is at the n-th position within
/// the array. After sorting, all elements below position
/// `n` will be smaller than or equal to the n-th-smallest
/// element, while all elements above position `n` will be
/// greater than or equal to the n-th-smallest element. The
/// value of the n-th-smallest element is returned.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
/// @param n     n-th smallest value will be returned.
///
/// @return Value of the n-th smallest array element.
///
/// @note This function is not `NaN`-safe and will modify
/// the original data array.

float nth_element_flt(float *data, const size_t size, const size_t n)
{
	float *l = data;
	float *m = data + size - 1;
	float *ptr = data + n;
	
	while(l < m)
	{
		float value = *ptr;
		float *i = l;
		float *j = m;
		
		do
		{
			while(*i < value) ++i;
			while(value < *j) --j;
			
			if(i <= j)
			{
				float tmp = *i;
				*i = *j;
				*j = tmp;
				++i;
				--j;
			}
		} while(i <= j);
		
		if(j < ptr) l = i;
		if(ptr < i) m = j;
	}
	
	return *ptr;
}



/// @brief Median
///
/// Calculates the median of the input data array.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
/// @param fast  If `true`, will return approximate median for
///              even-sized arrays. If `false`, the exact
///              median will be returned, which takes longer.
///
/// @return Median of the data array values.
///
/// @note This function is **not** `NaN`-safe and will modify
/// the original data array.

float median_flt(float *data, const size_t size, const bool fast)
{
	const float value = nth_element_flt(data, size, size / 2);
	return IS_ODD(size) || fast ? value : (value + max_flt(data, size / 2)) / 2.0;
}

/// @brief Median (safe)
///
/// Calculates the median of the input data array. Unlike median_flt(),
/// this function makes a copy of the data and does not alter the original
/// data array.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
/// @param fast  If `true`, will return approximate median for
///              even-sized arrays. If `false`, the exact
///              median will be returned, which takes longer.
///
/// @return Median of the data array values.
///
/// @note This function **is** `NaN`-safe.

float median_safe_flt(const float *data, const size_t size, const bool fast)
{
	float *data_copy = (float *)memory(MALLOC, size, sizeof(float));
	//memcpy(data_copy, data, size * sizeof(float));
	size_t size_copy = 0;
	for(const float *ptr = data + size; ptr --> data;) if(IS_NOT_NAN(*ptr)) data_copy[size_copy++] = *ptr;
	const float result = median_flt(data_copy, size_copy, fast);
	free(data_copy);
	return result;
}



/// @brief Median absolute deviation from value
///
/// Calculates the median absolute deviation (MAD) of the
/// data array values from a user-specified value. The MAD
/// is defined as
///
///   median(|x - value|)
///
/// where x denotes the data values from the input array.
/// Note that this function **is** `NaN`-safe and will **not**
/// modify the original data array.
///
/// @param data     Pointer to the data array to be sorted.
/// @param size     Size of the input array.
/// @param value    Value about which to calculate the MAD.
/// @param cadence  Can be set to > 1 to speed up algorithm.
/// @param range    Flux range to be used. Can be Negative
///                 (-1), Full (0) or Positive (1).
///
/// @return Median absolute deviation of the array values from
/// the specified data value.

float mad_val_flt(const float *data, const size_t size, const float value, const size_t cadence, const int range)
{
	// Create copy of data array with specified range and cadence
	const size_t data_copy_size = (range == 0) ? (size / cadence) : (size / (2 * cadence));
	float *data_copy = (float *)memory(MALLOC, data_copy_size, sizeof(float));
	
	// Some settings
	const float *ptr = data + size;
	float *ptr_copy = data_copy;
	size_t counter = 0;
	
	// Copy |*ptr - value| into array copy
	while((ptr -= cadence) > data && counter < data_copy_size)
	{
		if((range < 0 && *ptr < 0.0) || (range == 0 && IS_NOT_NAN(*ptr)) || (range > 0 && *ptr > 0.0))
		{
			ptr_copy[counter] = fabs(*ptr - value);
			++counter;
		}
	}
	
	// Determine median
	const float result = median_flt(data_copy, counter, false);
	
	// Clean up
	free(data_copy);
	
	return result;
}



/// @brief Median absolute deviation
///
/// Calculates the median absolute deviation (MAD) of the
/// data array values. The MAD is defined as
///
///   median(|x - median(x)|)
///
/// where x denotes the data values from the input array.
/// In the case of normally distributed random data values
/// the standard deviation of the data values about the
/// mean can be deduced by multiplying the MAD with the
/// constant `MAD_TO_STD`.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
///
/// @return Median absolute deviation of the array values.
///
/// @note This function is not `NaN`-safe and will modify
/// the original data array.

float mad_flt(float *data, const size_t size)
{
	return mad_val_flt(data, size, median_flt(data, size, false), 1, 0);
}



/// @brief Robust noise measurement for contiguous data array
///
/// Uses a robust way of measuring the noise within the
/// specified data array based on the median absolute
/// deviation (MAD). The following assumptions and conditions
/// apply:
///
///  * The noise is Gaussian and centred on zero.
///  * Only **negative** data points will be used in the noise
///    measurement.
///  * The algorithm is `NAN`-safe and discards all `NAN`s.
///
/// The result of the median of the absolute values of the
/// negative elements in the array will be returned. If no
/// valid data are found, `NaN` will instead be returned.
///
/// @param data Pointer to the data array.
/// @param size Size of the array.
///
/// @return Robust noise measurement of the data.

float robust_noise_flt(const float *data, const size_t size)
{
	// Allocate memory for 1D data copy
	float *data_copy = (float *)memory(MALLOC, size, sizeof(float));
	float *ptr_copy = data_copy;
	
	// Copy values of negative elements
	for(const float *ptr = data + size; ptr --> data;) if(*ptr < 0.0) *ptr_copy++ = *ptr;
	
	// Calculate pseudo-median
	const size_t size_copy = ptr_copy - data_copy;
	const float result = size_copy ? -MAD_TO_STD * nth_element_flt(data_copy, size_copy, size_copy / 2) : NAN;
	// NOTE: Multiplication by -1 because values are all negative.
	
	// Clean up and return result
	free(data_copy);
	return result;
}



/// @brief Robust noise measurement for contiguous data array
///
/// Uses a robust way of measuring the noise within the
/// specified data array based on the median absolute
/// deviation (MAD). The following assumptions and conditions
/// apply:
///
///  * The noise is Gaussian and centred on zero.
///  * Both **positive and negative** data points will be used
///    in the noise measurement.
///  * The algorithm is `NAN`-safe and discards all `NAN`s.
///
/// The result of the median of the absolute values of all
/// positive and negative elements in the array will be returned.
/// If no valid data are found, `NaN` will instead be returned.
///
/// @param data Pointer to the data array.
/// @param size Size of the array.
///
/// @return Robust noise measurement of the data.

float robust_noise_2_flt(const float *data, const size_t size)
{
	// Allocate memory for 1D data copy
	float *data_copy = (float *)memory(MALLOC, size, sizeof(float));
	float *ptr_copy = data_copy;
	
	// Copy values of all non-NaN elements
	for(const float *ptr = data + size; ptr --> data;)
	{
		if(*ptr < 0.0) *ptr_copy++ = -(*ptr);
		else if(*ptr >= 0.0) *ptr_copy++ = *ptr;
	}
	
	// Calculate pseudo-median
	const size_t size_copy = ptr_copy - data_copy;
	const float result = size_copy ? MAD_TO_STD * nth_element_flt(data_copy, size_copy, size_copy / 2) : NAN;
	
	// Clean up and return result
	free(data_copy);
	return result;
}



/// @brief Robust noise measurement in region of 3D array
///
/// Uses a robust way of measuring the noise within the
/// specified region of a 3D data array based on the median
/// absolute deviation (MAD). The following assumptions
/// and conditions apply:
///
///  * The array is contiguous in x and least contiguous
///    along the z axis.
///  * The noise is Gaussian and centred on zero.
///  * Only negative data points will be used in the noise
///    measurement.
///  * The algorithm is `NaN`-safe and discards all `NaN`s.
///
/// The result of the median of the absolute values of the
/// negative elements in the region will be returned. If no
/// valid data are found, `NaN` will instead be returned.
///
/// @param data  Pointer to the data array.
/// @param nx    Size of the data array in x.
/// @param ny    Size of the data array in y.
/// @param x1    Lower boundary of region in x.
/// @param x2    Upper boundary of region in x.
/// @param y1    Lower boundary of region in y.
/// @param y2    Upper boundary of region in y.
/// @param z1    Lower boundary of region in z.
/// @param z2    Upper boundary of region in z.
///
/// @return Robust noise measurement within the specified region.

float robust_noise_in_region_flt(const float *data, const size_t nx, const size_t ny, const size_t x1, const size_t x2, const size_t y1, const size_t y2, const size_t z1, const size_t z2)
{
	// Allocate memory for 1D data copy
	float *data_copy = (float *)memory(MALLOC, (x2 - x1 + 1) * (y2 - y1 + 1) * (z2 - z1 + 1), sizeof(float));
	float *ptr = data_copy;
	float value;
	
	// Copy absolute values of negative elements
	for(size_t z = z1; z <= z2; ++z)
	{
		for(size_t y = y1; y <= y2; ++y)
		{
			for(size_t x = x1; x <= x2; ++x)
			{
				value = data[x + nx * (y + ny * z)];
				if(value < 0.0) *ptr++ = -value;
			}
		}
	}
	
	// Calculate median
	const size_t size = ptr - data_copy;
	const float result = size ? MAD_TO_STD * nth_element_flt(data_copy, size, size / 2) : NAN;
	
	// Clean up
	free(data_copy);
	return result;
}



/// @brief Create histogram from data array
///
/// Generates a histogram of the data values in `data` with
/// the parameters as specified by the user. A pointer to
/// the generated histogram will be returned.
///
/// @param data      Pointer to the data array to be sorted.
/// @param size      Size of the input array.
/// @param n_bins    Number of bins of the histogram.
/// @param data_min  Lower flux limit of the histogram.
/// @param data_max  Upper flux limit of the histogram.
/// @param cadence   Cadence for generation of histogram. A
///                  cadence of N means that only every N-th
///                  data point will be used.
///
/// @return Pointer to the generated histogram.
///
/// @note It is the responsibility of the user to de-allocate
/// the memory occupied by the histogram once it is no longer
/// required (using `free()`) to ensure that there are no
/// memory leaks. This function is `NaN`-safe, and `NaN` values
/// will simply be ignored (due to `>=` and `<=` comparison).

size_t *create_histogram_flt(const float *data, const size_t size, const size_t n_bins, const float data_min, const float data_max, const size_t cadence)
{
	// Allocate memory
	size_t *histogram = (size_t *)memory(CALLOC, n_bins, sizeof(size_t));
	
	// Basic setup
	const float slope = (float)(n_bins - 1) / (data_max - data_min);
	const float *ptr  = data + size;
	const float tmp   = 0.5 - slope * data_min;  // The 0.5 is needed for correct rounding later on.
	
	// Generate histogram
	while((ptr -= cadence) > data)
	{
		if(*ptr >= data_min && *ptr <= data_max)
		{
			size_t bin = (size_t)(slope * *ptr + tmp);
			++histogram[bin];
		}
	}
	
	return histogram;
}



/// @brief Gaussian fit to histogram
///
/// Generates a histogram from the data values in the input
/// array and fits a Gaussian function to that histogram to
/// determine the standard deviation. This is useful for
/// measuring the noise level in the data cube.
///
/// The user can choose which pixels from the data cube
/// contribute to the histogram (positive, negative, or all
/// pixels and every N-th pixel as controlled by `cadence`).
/// For the purpose of fitting, the histogram range will be
/// set such that the second moment fills an optimal fraction
/// of the histogram width. A linear regression to the
/// logarithmic histogram, ln(h(x)) = a x^2 + b, will then
/// be performed to derive the standard deviation from the
/// slope of the fit. The advantage of this technique is
/// that the regression can be solved analytically, making
/// the algorithm extremely fast. The disadvantage is that
/// the weighting of histogram bins will be higher near the
/// centre of the flux distribution and lower in the wings,
/// although this can actually be an advantage, as the
/// wings are often affected by non-Gaussian effects such
/// as flux from actual objects.
///
/// @param data     Pointer to the data array to be sorted.
/// @param size     Size of the input array.
/// @param cadence  Cadence for generation of histogram. A
///                 cadence of N means that only every N-th
///                 data point will be used.
/// @param range    Flux range to use in histogram. Can be
///                 -1, 0 or +1 to use only negative pixels,
///                 all pixels, or only positive pixels,
///                 respectively.
///
/// @return Standard deviation from Gaussian fit.

float gaufit_flt(const float *data, const size_t size, const size_t cadence, const int range)
{
	// Determine maximum and minimum
	float data_max = 0.0;
	float data_min = 0.0;
	max_min_flt(data, size, &data_max, &data_min);
	
	if(data_min >= 0.0 || data_max <= 0.0)
	{
		warning("Maximum is not greater than minimum.");
		return NAN;
	}
	
	// Determine data range
	if(range < 0)
	{
		if(data_min >= 0.0)
		{
			warning("Minimum is not less than zero.");
			return NAN;
		}
		data_max = 0.0;
	}
	else if(range > 0)
	{
		if(data_max <= 0.0)
		{
			warning("Maximum is not greater than zero.");
			return NAN;
		}
		data_min = 0.0;
	}
	else
	{
		const float limit = fabs(data_min) < fabs(data_max) ? fabs(data_min) : fabs(data_max);
		data_min = -limit;
		data_max = limit;
	}
	
	// Create initial histogram
	const size_t n_bins = 101;
	size_t origin = n_bins / 2;
	if(range < 0) origin = n_bins - 1;
	else if(range > 0) origin = 0;
	const float inv_optimal_mom2 = 5.0 / n_bins;  // Require standard deviation to cover 1/5th of the histogram for optimal results
	size_t *histogram = create_histogram_flt(data, size, n_bins, data_min, data_max, cadence);
	
	// Calculate second moment
	float mom0 = 0.0;
	float mom1 = 0.0;
	float mom2 = 0.0;
	
	for(size_t i = n_bins; i--;)
	{
		mom0 += histogram[i];
		mom1 += histogram[i] * (float)i;
	}
	mom1 /= mom0;
	
	for(size_t i = n_bins; i--;) mom2 += histogram[i] * (mom1 - i) * (mom1 - i);
	mom2 = sqrt(mom2 / mom0);
	
	// Ensure that 2nd moment is equal to optimal_mom2 for optimal fitting
	if(range < 0) data_min *= mom2 * inv_optimal_mom2;
	else if(range > 0) data_max *= mom2 * inv_optimal_mom2;
	else
	{
		data_min *= mom2 * inv_optimal_mom2;
		data_max *= mom2 * inv_optimal_mom2;
	}
	
	// Regenerate histogram with new scaling
	free(histogram);
	histogram = create_histogram_flt(data, size, n_bins, data_min, data_max, cadence);
	
	// Fit Gaussian function using linear regression
	// (excluding first and last point in case of edge effects)
	float mean_x = 0.0;
	float mean_y = 0.0;
	size_t counter = 0;
	
	for(size_t i = n_bins - 1; i --> 1;)
	{
		if(histogram[i])
		{
			long int ii = i - origin;
			mean_x += (float)(ii * ii);
			mean_y += log((float)(histogram[i]));
			++counter;
		}
	}
	
	mean_x /= counter;
	mean_y /= counter;
	
	float upper_sum = 0.0;
	float lower_sum = 0.0;
	
	for(size_t i = n_bins - 1; i --> 1;)
	{
		if(histogram[i])
		{
			long int ii = i - origin;
			const float x = (float)(ii * ii);
			const float y = log((float)(histogram[i]));
			upper_sum += (x - mean_x) * (y - mean_y);
			lower_sum += (x - mean_x) * (x - mean_x);
		}
	}
	
	// Determine standard deviation from slope
	const float sigma = sqrt(-0.5 * lower_sum / upper_sum) * (data_max - data_min) / (n_bins - 1);
	
	// Clean up
	free(histogram);
	
	return sigma;
}



/// @brief Skewness
///
/// Calculates the skewness of the input data array values.
/// Skewness is defined as
///
///   mom3 / sqrt(mom2^3)
///
/// where mom2 and mom3 are the second and third moment of
/// the data as returned by the moment() and moments()
/// functions. The skewness of normally distributed data
/// with a mean of zero should be 0.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
///
/// @return Skewness of the array values. `NaN` will be returned for
///         `NaN`-only input arrays.

double skewness_flt(const float *data, const size_t size)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	return m3 / sqrt(m2 * m2 * m2);
}



/// @brief Kurtosis
///
/// Calculates the kurtosis of the input data array values.
/// Kurtosis is defined as
///
///   mom4 / mom2^2
///
/// where mom2 and mom4 are the second and fourth moment of
/// the data as returned by the moment() and moments()
/// functions. The kurtosis of normally distributed data
/// with a mean of zero should be 3.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
///
/// @return Kurtosis of the array values. `NaN` will be returned for
///         `NaN`-only input arrays.

double kurtosis_flt(const float *data, const size_t size)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	return m4 / (m2 * m2);
}



/// @brief Skewness and kurtosis
///
/// Calculates the skewness and kurtosis of the input data
/// array simultaneously. This is faster than calculating
/// them individually in situations where both are needed.
/// See the skewness() and kurtosis() functions for more
/// details on how these parameters are defined. Both skew
/// and kurt will be set to `NaN` if the array contains only
/// `NaN` values.
///
/// @param data  Pointer to the data array to be sorted.
/// @param size  Size of the input array.
/// @param skew  Pointer to variable for returning skewness.
/// @param kurt  Pointer to variable for returning kurtosis.

void skew_kurt_flt(const float *data, const size_t size, double *skew, double *kurt)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	*skew = m3 / sqrt(m2 * m2 * m2);
	*kurt = m4 / (m2 * m2);
	return;
}



/// @brief 1D boxcar filter
///
/// Applies a boxcar filter to the data array. `NaN` values
/// will be set to 0 prior to filtering only if `replace_nan`
/// is set to `true` (this is to avoid a check within this
/// function for reasons of speed). Values outside of the
/// boundaries of the array are assumed to be 0.
///
/// @param data           Pointer to data array to be filtered.
/// @param data_copy      Pointer to data array to be used for
///                       storing a copy of the data during
///                       filtering. Its size must be equal to
///                       `size` + 2 * `filter_radius`.
/// @param size           Size of input array.
/// @param filter_radius  Radius of boxcar filter.
///
/// @note This function will modify the original data array.

void filter_boxcar_1d_flt(float *data, float *data_copy, const size_t size, const size_t filter_radius)
{
	// Define filter size
	const size_t filter_size = 2 * filter_radius + 1;
	const float inv_filter_size = 1.0 / filter_size;
	size_t i;
	
	// Make copy of data, taking care of NaN
	for(i = size; i--;) data_copy[filter_radius + i] = FILTER_NAN(data[i]);
	
	// Fill overlap regions with 0
	for(i = filter_radius; i--;) data_copy[i] = data_copy[size + filter_radius + i] = 0.0;
	
	// Apply boxcar filter to last data point
	data[size - 1] = 0.0;
	for(i = filter_size; i--;) data[size - 1] += data_copy[size + i - 1];
	data[size - 1] *= inv_filter_size;
	
	// Recursively apply boxcar filter to all previous data points
	for(i = size - 1; i--;) data[i] = data[i + 1] + (data_copy[i] - data_copy[filter_size + i]) * inv_filter_size;
	
	return;
}



/// @brief 2D Gaussian filter
///
/// Applies a pseudo-Gaussian filter to the two-dimensional
/// data array by running a series of `n_iter` boxcar filters
/// of radius `filter_radius` across the data array in both
/// dimensions. The function optimal_filter_size() can
/// be used to automatically determine the required number
/// of iterations and boxcar filter radius for a given
/// standard deviation of the Gaussian.
///
/// For reasons of speed several data arrays must have been
/// pre-allocated and passed on to this function:
///
/// * `data_copy`: Used to store a single column of the input
///   data. Must be of size `size_y`.
/// * `data_row`:  Used to store a copy of the data passed on
///   to the boxcar filter. Must be of size `size_x` +
///   2 * `filter_radius`.
/// * `data_col`:  Used to store a copy of the data passed on
///   to the boxcar filter. Must be of size `size_y` +
///   2 * `filter_radius`.
///
/// The sole purpose of having these array created externally
/// and then passed on to the function is to improve
/// the speed of the algorithm in cases where it needs to
/// be invoked a large number of times.
///
/// Note that the function will only be able to approximate
/// a Gaussian, with the approximation becoming better for
/// a larger number of iterations of the boxcar filter. 3-4
/// iterations should already provide a reasonable approximation
/// for the purpose of image filtering. Also note that the value
/// of the standard deviation of the Gaussian will only be
/// approximated by the optimal number of iterations and filter
/// radius, usually to within several percent.
///
/// @param data           Pointer to data array to be
///                       filtered.
/// @param data_copy      Pointer to data array to be used
///                       for storing a single column of the
///                       data array. Its size must be equal
///                       to `size_y`.
/// @param data_row       Pointer to data array to be used by
///                       the boxcar filter to be employed.
///                       Its size must be equal to `size_x` +
///                       2 * `filter_radius`.
/// @param data_col       Pointer to data array to be used by
///                       the boxcar filter to be employed.
///                       Its size must be equal to `size_y` +
///                       2 * `filter_radius`.
/// @param size_x         Size of the first dimension of the
///                       input data array.
/// @param size_y         Size of the second dimension of the
///                       input data array.
/// @param n_iter         Number of iterations of boxcar
///                       filtering to be carried out.
/// @param filter_radius  Radius of the boxcar filter to be
///                       applied. The filter width will be
///                       defined as 2 * `filter_radius` + 1.

void filter_gauss_2d_flt(float *data, float *data_copy, float *data_row, float *data_col, const size_t size_x, const size_t size_y, const size_t n_iter, const size_t filter_radius)
{
	// Set up a few variables
	const size_t size_xy = size_x * size_y;
	float *ptr = data + size_xy;
	float *ptr2;
	
	// Run row filter (along x-axis)
	// This is straightforward, as the data are contiguous in x.
	while(ptr > data)
	{
		ptr -= size_x;
		for(size_t i = n_iter; i--;) filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
	}
	
	// Run column filter (along y-axis)
	// This is more complicated, as the data are non-contiguous in y.
	for(size_t x = size_x; x--;)
	{
		// Copy data into column array
		ptr = data + size_xy - size_x + x;
		ptr2 = data_copy + size_y;
		while(ptr2 --> data_copy)
		{
			*ptr2 = *ptr;
			ptr -= size_x;
		}
		
		// Apply filter
		for(size_t i = n_iter; i--;) filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);
		
		// Copy column array back into data array
		ptr = data + size_xy - size_x + x;
		ptr2 = data_copy + size_y;
		while(ptr2 --> data_copy)
		{
			*ptr = *ptr2;
			ptr -= size_x;
		}
	}
	
	return;
}



/// @brief Shift and subtract data from itself
///
/// The function will subtract a shifted copy of the data
/// array from itself. Processing will start from the end
/// such that `data[size - 1] -= data[size - 1 - shift]`,
/// etc., until position `data + shift` is reached. The first
/// `shift` elements will be left unchanged.
///
/// @param data   Data array to be processed.
/// @param size   Size of the data array.
/// @param shift  Number of positions by which to
///               shift before subtraction.
///
/// @note This function will modify the input data array.

void shift_and_subtract_flt(float *data, const size_t size, const size_t shift)
{
	for(float *ptr = data + size; ptr --> data + shift;) *ptr -= *(ptr - shift);
	return;
}



/// @brief Determine optimal boxcar filter size for Gaussian filter
///
/// Based on the specified standard deviation of the Gaussian
/// filter, this function will determine the optimal radius and
/// number of iterations for the boxcar filter that is to be used
/// to approximate the Gaussian kernel. Both `filter_radius` and
/// `n_iter` are pointers and will be set to the correct values by
/// this function. They can then be fed into the Gaussian filtering
/// function, filter_gauss_2d_flt().
///
/// @param sigma          Standard deviation of the Gaussian.
/// @param filter_radius  Radius of the boxcar filter that is
///                       to be used in approximation.
/// @param n_iter         Number of iterations of boxcar filtering
///                       to be used in approximation.

void optimal_filter_size_flt(const double sigma, size_t *filter_radius, size_t *n_iter)
{
	*n_iter = 0;
	*filter_radius = 0;
	double tmp = -1.0;
	size_t i;
	
	for(i = BOXCAR_MIN_ITER; i <= BOXCAR_MAX_ITER; ++i)
	{
		const double radius = sqrt((3.0 * sigma * sigma / i) + 0.25) - 0.5;
		const double diff = fabs(radius - floor(radius + 0.5));
		
		if(tmp < 0.0 || diff < tmp)
		{
			tmp = diff;
			*n_iter = i;
			*filter_radius = (size_t)(radius + 0.5);
		}
	}
	
	// Print some information
	/*const double sigma_approx = sqrt((double)(*n_iter) * ((2.0 * (double)(*filter_radius) + 1.0) * (2.0 * (double)(*filter_radius) + 1.0) - 1.0) / 12.0);
	message("Requested filter size:    sigma = %.2f\n", sigma);
	message("Approximated filter size: sigma = %.2f\n", sigma_approx);
	message("  using N = %zu and R = %zu\n", *n_iter, *filter_radius);*/
	
	return;
}



/// @brief Fit ellipse to moment-0 map
///
/// Function for fitting an ellipse to the specified moment-0 map.
/// The moment and count map arrays must be contiguous in x. The
/// source centroid must be relative to the boundaries of the
/// moment map. The ellipse fit is carried out using 2nd-order
/// moment analysis.
///
/// Two types of ellipse are fitted: the first one is a fit to the
/// flux-weighted positive pixels in the moment maps, while the
/// second fit is to all pixels greater than 3 times the specified
/// RMS noise level with equal weights. The results will be written
/// to the user-specified pointer variables.
///
/// @param moment_map  Pointer to moment-0 map
/// @param count_map   Pointer to map containing number of channels
///                    for each pixel in moment map.
/// @param size_x      Size of moment map in x.
/// @param size_y      Size of moment map in y.
/// @param centroid_x  x position of source centroid relative to
///                    moment map boundary.
/// @param centroid_y  y position of source centroid relative to
///                    moment map boundary.
/// @param rms         RMS noise level of the data.
/// @param ell_maj     Pointer to variable for ellipse major axis.
/// @param ell_min     Pointer to variable for ellipse minor axis.
/// @param ell_pa      Pointer to variable for ellipse position angle.
/// @param ell3s_maj   Same, but for 3-sigma ellipse.
/// @param ell3s_min   Same, but for 3-sigma ellipse.
/// @param ell3s_pa    Same, but for 3-sigma ellipse.

void moment_ellipse_fit_flt(const float *moment_map, const size_t *count_map, const size_t size_x, const size_t size_y, const float centroid_x, const float centroid_y, const float rms, float *ell_maj, float *ell_min, float *ell_pa, float *ell3s_maj, float *ell3s_min, float *ell3s_pa)
{
	float ell_momX    = 0.0;
	float ell_momY    = 0.0;
	float ell_momXY   = 0.0;
	float ell_sum     = 0.0;
	float ell3s_momX  = 0.0;
	float ell3s_momY  = 0.0;
	float ell3s_momXY = 0.0;
	float ell3s_sum   = 0.0;
	
	for(size_t y = 0; y < size_y; ++y)
	{
		for(size_t x = 0; x < size_x; ++x)
		{
			float value = moment_map[x + size_x * y];
			size_t count = count_map [x + size_x * y];
			
			if(value > 0.0)
			{
				// Calculate moments for flux-weighted ellipse fitting
				// NOTE: Only positive pixels considered here!
				ell_momX  += ((float)x - centroid_x) * ((float)x - centroid_x) * value;
				ell_momY  += ((float)y - centroid_y) * ((float)y - centroid_y) * value;
				ell_momXY += ((float)x - centroid_x) * ((float)y - centroid_y) * value;
				ell_sum += value;
				
				if(value > 3.0 * rms * sqrt((float)count))
				{
					// Calculate moments for 3-sigma ellipse fitting
					ell3s_momX  += ((float)x - centroid_x) * ((float)x - centroid_x);
					ell3s_momY  += ((float)y - centroid_y) * ((float)y - centroid_y);
					ell3s_momXY += ((float)x - centroid_x) * ((float)y - centroid_y);
					ell3s_sum += 1.0;
				}
			}
		}
	}
	
	if(ell_sum > 0.0)
	{
		ell_momX  /= ell_sum;
		ell_momY  /= ell_sum;
		ell_momXY /= ell_sum;
		*ell_pa  = 0.5 * atan2(2.0 * ell_momXY, ell_momX - ell_momY);
		*ell_maj = sqrt(2.0 * (ell_momX + ell_momY + sqrt((ell_momX - ell_momY) * (ell_momX - ell_momY) + 4.0 * ell_momXY * ell_momXY)));
		*ell_min = sqrt(2.0 * (ell_momX + ell_momY - sqrt((ell_momX - ell_momY) * (ell_momX - ell_momY) + 4.0 * ell_momXY * ell_momXY)));
		*ell_pa  = *ell_pa   * 180.0 / M_PI - 90.0;
		while(*ell_pa < -90.0) *ell_pa += 180.0;
	}
	
	if(ell3s_sum > 0.0)
	{
		ell3s_momX  /= ell3s_sum;
		ell3s_momY  /= ell3s_sum;
		ell3s_momXY /= ell3s_sum;
		*ell3s_pa  = 0.5 * atan2(2.0 * ell3s_momXY, ell3s_momX - ell3s_momY);
		*ell3s_maj = sqrt(2.0 * (ell3s_momX + ell3s_momY + sqrt((ell3s_momX - ell3s_momY) * (ell3s_momX - ell3s_momY) + 4.0 * ell3s_momXY * ell3s_momXY)));
		*ell3s_min = sqrt(2.0 * (ell3s_momX + ell3s_momY - sqrt((ell3s_momX - ell3s_momY) * (ell3s_momX - ell3s_momY) + 4.0 * ell3s_momXY * ell3s_momXY)));
		*ell3s_pa = *ell3s_pa * 180.0 / M_PI - 90.0;
		while(*ell3s_pa < -90.0) *ell3s_pa += 180.0;
	}
	
	// NOTE: Converting PA from radians to degrees.
	// NOTE: Subtracting 90 deg from PA, because astronomers like to have 0 deg pointing up.
	// NOTE: This means that PA will no longer have the mathematically correct orientation!
	// NOTE: PA should then be between -90 deg (right) and +90 deg (left), with 0 deg pointing up.
	// WARNING: PA will be relative to the pixel grid, not the coordinate system!
	
	return;
}



/// @brief Determine w20 and w50 line widths from spectrum
///
/// Function for measuring the w20 and w50 line width of the specified
/// spectrum. Both are measured by moving inwards from the edges of
/// the spectrum until the point where the flux density increases
/// above 20% (or 50%) of the peak flux density. Linear interpolation
/// is used to improved the accuracy of the measured line widths. The
/// results will be written to the user-specified `w20` and `w50`
/// pointers.
///
/// @param spectrum  Pointer to spectrum.
/// @param size      Size of spectrum.
/// @param w20       Pointer to variable for w20 line width.
/// @param w50       Pointer to variable for w50 line width.

void spectral_line_width_flt(const float *spectrum, const size_t size, float *w20, float *w50)
{
	// Auxiliary parameters
	size_t index;
	float maximum = -INFINITY;
	
	// Determine maximum
	for(size_t i = size; i--;) if(spectrum[i] > maximum) maximum = spectrum[i];
	
	// Determine w20 from spectrum (moving inwards)
	for(index = 0; index < size && spectrum[index] < 0.2 * maximum; ++index);
	if(index < size)
	{
		*w20 = (float)index;
		if(index > 0) *w20 -= (spectrum[index] - 0.2 * maximum) / (spectrum[index] - spectrum[index - 1]);
		for(index = size - 1; index < size && spectrum[index] < 0.2 * maximum; --index); // index is unsigned
		*w20 = (float)index - *w20;
		if(index < size - 1) *w20 += (spectrum[index] - 0.2 * maximum) / (spectrum[index] - spectrum[index + 1]);
	}
	else
	{
		*w20 = 0.0;
		warning("Failed to measure w20.");
	}
	
	// Determine w50 from spectrum (moving inwards)
	for(index = 0; index < size && spectrum[index] < 0.5 * maximum; ++index);
	if(index < size)
	{
		*w50 = (float)index;
		if(index > 0) *w50 -= (spectrum[index] - 0.5 * maximum) / (spectrum[index] - spectrum[index - 1]);
		for(index = size - 1; index < size && spectrum[index] < 0.5 * maximum; --index); // index is unsigned
		*w50 = (float)index - *w50;
		if(index < size - 1) *w50 += (spectrum[index] - 0.5 * maximum) / (spectrum[index] - spectrum[index + 1]);
	}
	else
	{
		*w50 = 0.0;
		warning("Failed to measure w50.");
	}
	
	return;
}



/// @brief Determine wm50 line width from spectrum
///
/// Function for measuring the wm50 line width from a spectrum
/// following the algorithm described in Courtois et al. (2009).
/// If the measurement fails for some reason, then a value of 0
/// will be returned and a warning message printed.
///
/// @param spectrum  Pointer to spectrum.
/// @param size      Size of spectrum.
///
/// @return Value of wm50.

double wm50_line_width_flt(const float *spectrum, const size_t size)
{
	// Determine integrated flux
	double flux = 0.0;
	for(size_t i = size; i--;) if(IS_NOT_NAN(spectrum[i])) flux += (double)(spectrum[i]);
	
	// Sanity check
	if(flux <= 0.0 || IS_NAN(flux))
	{
		warning("Failed to measure wm50.");
		return 0.0;
	}
	
	// Determine lower boundary
	size_t lower = 0;
	double flux_lower = spectrum[lower];
	while(flux_lower < 0.05 * flux) flux_lower += spectrum[++lower];
	
	// Determine upper boundary
	size_t upper = size - 1;
	double flux_upper = spectrum[upper];
	while(flux_upper < 0.05 * flux) flux_upper += spectrum[--upper];
	
	// Sanity check
	if(upper <= lower)
	{
		warning("Failed to measure wm50.");
		return 0.0;
	}
	
	// Linearly interpolate for higher precision
	double lower_intp = (double)lower - (flux_lower - 0.05 * flux) / spectrum[lower];
	double upper_intp = (double)upper + (flux_upper - 0.05 * flux) / spectrum[upper];
	
	// Determine mean flux per channel within window
	const double flux_mean = 0.9 * flux / (upper_intp - lower_intp);
	
	// Determine lower and upper boundaries for wm50
	lower = 0;
	upper = size - 1;
	while(spectrum[lower] < 0.5 * flux_mean) ++lower;
	while(spectrum[upper] < 0.5 * flux_mean) --upper;
	
	// Linearly interpolate for higher precision
	lower_intp = (double)lower;
	upper_intp = (double)upper;
	if(lower > 0)        lower_intp -= (spectrum[lower] - 0.5 * flux_mean) / (spectrum[lower] - spectrum[lower - 1]);
	if(upper < size - 1) upper_intp += (spectrum[upper] - 0.5 * flux_mean) / (spectrum[upper] - spectrum[upper + 1]);
	
	return upper_intp - lower_intp;
}



/// @brief Determine kinematic major axis of galaxy
///
/// Function for measuring the position angle of the kinematic
/// major axis of a galaxy from an array of flux-weighted centroid
/// measurements per spectral channel. The position angle is
/// determined by fitting a straight line to the centroid values using
/// orthogonal (Deming) regression. Missing or invalid centroids in
/// certain channels can be indicated by setting the corresponding
/// value in the sum array to zero. The resulting position angle
/// will point towards the side of the galaxy that sits at the upper
/// end of the channel range occupied by the galaxy.
///
/// @param centroidX  Array of x centroid positions per channel.
/// @param centroidY  Array of y centroid positions per channel.
/// @param sum        Array of summed flux densities per channel.
/// @param size       Size of `centroidX`, `centroidY` and `sum` arrays.
/// @param first      Index of first valid centroid position.
/// @param last       Index of last valid centroid position.
///
/// @return Position angle of kinematic major axis.

float kin_maj_axis_flt(const float *centroidX, const float *centroidY, const float *sum, const size_t size, const size_t first, const size_t last)
{
	// Fit a straight line to the set of centroids (using user-defined weights)
	float sumW   = 0.0;
	float sumX   = 0.0;
	float sumY   = 0.0;
	float sumXX  = 0.0;
	float sumYY  = 0.0;
	float sumXY  = 0.0;
	float weight = 1.0;
	
	// Using Deming (orthogonal) regression:
	for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			// Set the desired weights here (defaults to 1 otherwise):
			weight = sum[i] * sum[i];
			
			sumW += weight;
			sumX += weight * centroidX[i];
			sumY += weight * centroidY[i];
		}
	}
	
	sumX /= sumW;
	sumY /= sumW;
	
	for(size_t i = 0; i < size; ++i)
	{
		if(sum[i] > 0.0)
		{
			// Set the desired weights here (defaults to 1 otherwise):
			weight = sum[i] * sum[i];
			
			sumXX += weight * (centroidX[i] - sumX) * (centroidX[i] - sumX);
			sumYY += weight * (centroidY[i] - sumY) * (centroidY[i] - sumY);
			sumXY += weight * (centroidX[i] - sumX) * (centroidY[i] - sumY);
		}
	}
	
	const float slope = (sumYY - sumXX + sqrt((sumYY - sumXX) * (sumYY - sumXX) + 4.0 * sumXY * sumXY)) / (2.0 * sumXY);
	//const float inter = sumY - slope * sumX;
	
	// Calculate position angle of kinematic major axis:
	float pa = atan(slope);
	
	// Check orientation of approaching/receding side of disc:
	// (by averaging the lower and upper half of centroids)
	float x_first = 0.0;
	float x_last  = 0.0;
	float y_first = 0.0;
	float y_last  = 0.0;
	size_t counter_first = 0;
	size_t counter_last  = 0;
	float fullAngle     = 0.0;
	
	for(size_t i = 0; i < size / 2; ++i)
	{
		if(sum[i] > 0.0)
		{
			x_first += centroidX[i];
			y_first += centroidY[i];
			++counter_first;
		}
		if(sum[size - i - 1] > 0.0)
		{
			x_last += centroidX[size - i - 1];
			y_last += centroidY[size - i - 1];
			++counter_last;
		}
	}
	
	if(counter_first > 0 && counter_last > 0)
	{
		// User lower and upper half averages by default
		x_first /= counter_first;
		y_first /= counter_first;
		x_last  /= counter_last;
		y_last  /= counter_last;
		fullAngle = atan2(y_last - y_first, x_last - x_first);
	}
	else
	{
		// Otherwise use first and last point as a back-up option
		fullAngle = atan2(centroidY[last] - centroidY[first], centroidX[last] - centroidX[first]);
	}
	
	// Flip PA by 180° if necessary:
	const float difference = fabs(atan2(sin(fullAngle) * cos(pa) - cos(fullAngle) * sin(pa), cos(fullAngle) * cos(pa) + sin(fullAngle) * sin(pa)));
	if(difference > M_PI / 2.0) pa += M_PI;
	
	// Correct for astronomers' favourite definition of PA:
	pa = 180.0 * pa / M_PI - 90.0;
	while(pa <    0.0) pa += 360.0;
	while(pa >= 360.0) pa -= 360.0;
	// NOTE: PA should now be between 0° (pointing up) and 360° and refer to the side of the
	//       disc that occupies the upper end of the channel range covered by the source.
	//       PAs are relative to the coordinate system imposed by the centroid arrays.
	
	return pa;
}
