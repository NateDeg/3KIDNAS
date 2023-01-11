// ____________________________________________________________________ //
//                                                                      //
// SoFiA 2.5.1 (sofia.c) - Source Finding Application                   //
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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>

// WARNING: The following is compiler-specific and
//          not part of the C standard library.
#ifdef _OPENMP
	#include <omp.h>
#endif

// WARNING: The following will only work on POSIX-compliant
//          systems, but is needed for mkdir().
#include <errno.h>
#include <sys/stat.h>

#include "src/common.h"
#include "src/Table.h"
#include "src/Path.h"
#include "src/Array_dbl.h"
#include "src/Array_siz.h"
#include "src/Map.h"
#include "src/Matrix.h"
#include "src/Parameter.h"
#include "src/Catalog.h"
#include "src/WCS.h"
#include "src/DataCube.h"
#include "src/LinkerPar.h"



// ----------------------------------------------------------------- //
// This file contains the actual SoFiA pipeline that will read user  //
// parameters and data files, call the requested processing modules  //
// and write out catalogues and images.                              //
// ----------------------------------------------------------------- //

int main(int argc, char **argv)
{
	// ---------------------------- //
	// Record starting time         //
	// ---------------------------- //
	
	const time_t start_time = time(NULL);
	const clock_t start_clock = clock();
	
	
	
	// ---------------------------- //
	// A few global definitions     //
	// ---------------------------- //
	
	const char *noise_stat_name[] = {"standard deviation", "median absolute deviation", "Gaussian fit to flux histogram", "mean", "median"};
	const char *flux_range_name[] = {"negative", "full", "positive"};
	double global_rms = 1.0;
	#ifdef _OPENMP
		const int n_cpu_cores = omp_get_num_procs();
	#endif
	
	
	
	// ---------------------------- //
	// Print basic information      //
	// ---------------------------- //
	
	status("Pipeline started");
	message("Using:    Source Finding Application (SoFiA)");
	message("Version:  %s (%s)", SOFIA_VERSION, SOFIA_CREATION_DATE);
	#ifdef _OPENMP
		message("CPU:      %d %s available", n_cpu_cores, n_cpu_cores == 1 ? "thread" : "threads");
	#else
		message("CPU:      OpenMP disabled");
	#endif
	message("Time:     %s", ctime(&start_time));
	
	
	
	// ---------------------------- //
	// Check command line arguments //
	// ---------------------------- //
	
	ensure(argc >= 2, ERR_USER_INPUT, "Unexpected number of command line arguments.\nUsage: %s <parameter_file>", argv[0]);
	
	
	
	// ---------------------------- //
	// Set default parameters       //
	// ---------------------------- //
	
	status("Loading parameter settings");
	
	message("Activating SoFiA default parameter settings.");
	
	Parameter *par = Parameter_new(false);
	Parameter_default(par);
	
	
	
	// ---------------------------- //
	// Load user parameters         //
	// ---------------------------- //
	
	message("Loading user-specified parameters.");
	
	for(int i = 1; i < argc; ++i)
	{
		bool unknown_parameter = false;
		
		if(strchr(argv[i], '=') == NULL)
		{
			// Load parameter file
			message("- Loading user parameter file: %s", argv[i]);
			Parameter_load(par, argv[i], PARAMETER_UPDATE);
		}
		else
		{
			// Read command-line setting
			String *key = String_new("");
			String *val = String_new("");
			
			String_set_delim(key, argv[i], '=', true, true);
			String_set_delim(val, argv[i], '=', true, false);
			
			if(Parameter_exists(par, String_get(key), NULL))
			{
				message("- Setting parameter: %s = %s", String_get(key), String_get(val));
				Parameter_set(par, String_get(key), String_get(val));
			}
			else
			{
				warning("Unknown parameter: \'%s\'.", String_get(key));
				unknown_parameter = true;
			}
			
			String_delete(key);
			String_delete(val);
		}
		
		ensure(!unknown_parameter || !Parameter_get_bool(par, "pipeline.pedantic"), ERR_USER_INPUT, "Unknown parameter settings encountered. Please check\n       your input or change \'pipeline.pedantic\' to \'false\'.");
	}
	
	
	
	// ---------------------------- //
	// Set number of threads        //
	// ---------------------------- //
	
	#ifdef _OPENMP
		const int n_threads = Parameter_get_int(par, "pipeline.threads");
		if(n_threads > 0)
		{
			if(n_threads < n_cpu_cores)
			{
				omp_set_num_threads(n_threads);
				message("Using %d out of %d available CPU threads.\n", n_threads, n_cpu_cores);
			}
			else
			{
				omp_set_num_threads(n_cpu_cores);
				message("Using all %d available CPU threads.\n", n_cpu_cores);
			}
		}
		else
		{
			char *env_omp_num_threads = getenv("OMP_NUM_THREADS");
			message("Number of CPU threads controlled by OMP_NUM_THREADS = %s.\n", env_omp_num_threads == NULL ? "[undefined]" : env_omp_num_threads);
		}
	#else
		warning("Multi-threading is currently disabled. To enable it, please re-\n         install SoFiA with the '-fopenmp' option.");
	#endif
	
	
	
	// ---------------------------- //
	// Extract important settings   //
	// ---------------------------- //
	
	const bool verbosity         = Parameter_get_bool(par, "pipeline.verbose");
	const bool use_region        = strlen(Parameter_get_str(par, "input.region")) ? true : false;
	const bool use_gain          = strlen(Parameter_get_str(par, "input.gain"))   ? true : false;
	const bool use_noise         = strlen(Parameter_get_str(par, "input.noise"))  ? true : false;
	const bool use_weights       = strlen(Parameter_get_str(par, "input.weights"))? true : false;
	const bool use_mask          = strlen(Parameter_get_str(par, "input.mask"))   ? true : false;
	const bool use_invert        = Parameter_get_bool(par, "input.invert");
	      bool use_flagging      = strlen(Parameter_get_str(par, "flag.region"))  ? true : false;
	const bool use_flagging_cat  = strlen(Parameter_get_str(par, "flag.catalog")) ? true : false;
	const bool autoflag_log      = Parameter_get_bool(par, "flag.log");
	const bool use_cont_sub      = Parameter_get_bool(par, "contsub.enable");
	const bool use_noise_scaling = Parameter_get_bool(par, "scaleNoise.enable");
	const bool use_local_scaling = (strcmp(Parameter_get_str(par, "scaleNoise.mode"), "local") == 0);
	const bool use_sc_scaling    = Parameter_get_bool(par, "scaleNoise.scfind");
	const bool use_ripple_filter = Parameter_get_bool(par, "rippleFilter.enable");
	const bool use_scfind        = Parameter_get_bool(par, "scfind.enable");
	const bool use_threshold     = Parameter_get_bool(par, "threshold.enable");
	const bool use_linker        = Parameter_get_bool(par, "linker.enable");
	const bool linker_pos_pix    = Parameter_get_bool(par, "linker.positivity");
	const bool keep_negative     = Parameter_get_bool(par, "linker.keepNegative");
	const bool use_reliability   = Parameter_get_bool(par, "reliability.enable");
	const bool use_rel_plot      = Parameter_get_bool(par, "reliability.plot");
	const bool use_rel_cat       = strlen(Parameter_get_str(par, "reliability.catalog")) ? true : false;
	const bool use_rel_debug     = Parameter_get_bool(par, "reliability.debug");
	const bool use_mask_dilation = Parameter_get_bool(par, "dilation.enable");
	const bool use_parameteriser = Parameter_get_bool(par, "parameter.enable");
	const bool use_wcs           = Parameter_get_bool(par, "parameter.wcs");
	const bool use_physical      = Parameter_get_bool(par, "parameter.physical");
	const bool use_pos_offset    = Parameter_get_bool(par, "parameter.offset");
	
	const bool write_ascii       = Parameter_get_bool(par, "output.writeCatASCII");
	const bool write_xml         = Parameter_get_bool(par, "output.writeCatXML");
	const bool write_sql         = Parameter_get_bool(par, "output.writeCatSQL");
	const bool write_noise       = Parameter_get_bool(par, "output.writeNoise");
	const bool write_filtered    = Parameter_get_bool(par, "output.writeFiltered");
	const bool write_mask        = Parameter_get_bool(par, "output.writeMask");
	const bool write_mask2d      = Parameter_get_bool(par, "output.writeMask2d");
	const bool write_rawmask     = Parameter_get_bool(par, "output.writeRawMask");
	const bool write_moments     = Parameter_get_bool(par, "output.writeMoments");
	const bool write_cubelets    = Parameter_get_bool(par, "output.writeCubelets");
	const bool overwrite         = Parameter_get_bool(par, "output.overwrite");
	
	const double thresh_mom      = Parameter_get_flt(par, "output.thresholdMom12");
	const double rel_threshold   = Parameter_get_flt(par, "reliability.threshold");
	const double rel_snr_min     = Parameter_get_flt(par, "reliability.minSNR");
	const size_t rel_min_pix     = Parameter_get_int(par, "reliability.minPixels");
	
	unsigned int autoflag_mode = 0;
	if     (strcmp(Parameter_get_str(par, "flag.auto"), "channels") == 0) autoflag_mode = 1;
	else if(strcmp(Parameter_get_str(par, "flag.auto"), "pixels")   == 0) autoflag_mode = 2;
	else if(strcmp(Parameter_get_str(par, "flag.auto"), "true")     == 0) autoflag_mode = 3;
	
	// Statistic and range for noise measurement in
	// noise scaling, S+C finder and threshold finder
	noise_stat sn_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "scaleNoise.statistic"), "mad") == 0) sn_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "scaleNoise.statistic"), "gauss") == 0) sn_statistic = NOISE_STAT_GAUSS;
	
	int sn_range = 0;
	if(strcmp(Parameter_get_str(par, "scaleNoise.fluxRange"), "negative") == 0) sn_range = -1;
	else if(strcmp(Parameter_get_str(par, "scaleNoise.fluxRange"), "positive") == 0) sn_range = 1;
	
	noise_stat sc_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "scfind.statistic"), "mad") == 0) sc_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "scfind.statistic"), "gauss") == 0) sc_statistic = NOISE_STAT_GAUSS;
	
	int sc_range = 0;
	if(strcmp(Parameter_get_str(par, "scfind.fluxRange"), "negative") == 0) sc_range = -1;
	else if(strcmp(Parameter_get_str(par, "scfind.fluxRange"), "positive") == 0) sc_range = 1;
	
	noise_stat tf_statistic = NOISE_STAT_STD;
	if(strcmp(Parameter_get_str(par, "threshold.statistic"), "mad") == 0) tf_statistic = NOISE_STAT_MAD;
	else if(strcmp(Parameter_get_str(par, "threshold.statistic"), "gauss") == 0) tf_statistic = NOISE_STAT_GAUSS;
	
	int tf_range = 0;
	if(strcmp(Parameter_get_str(par, "threshold.fluxRange"), "negative") == 0) tf_range = -1;
	else if(strcmp(Parameter_get_str(par, "threshold.fluxRange"), "positive") == 0) tf_range = 1;
	
	// Ripple filter statistic
	int ripple_filter_statistic = NOISE_STAT_MEDIAN;
	if(strcmp(Parameter_get_str(par, "rippleFilter.statistic"), "mean") == 0) ripple_filter_statistic = NOISE_STAT_MEAN;
	
	// Noise and weights sanity check
	if(use_noise && use_weights) warning("Applying both a weights cube and a noise cube.");
	
	// Negative detections sanity check
	ensure(!(keep_negative && use_reliability), ERR_USER_INPUT, "With the reliability filter enabled, negative detections would always\n       be discarded irrespective of the value of linker.keepNegative! Please\n       set linker.keepNegative = false or disable reliability filtering.");
	ensure(!(linker_pos_pix && use_reliability), ERR_USER_INPUT, "With linker.positivity = true, there would be no negative\n       detections for the reliability filter to work on. Please either\n       disable the reliability filter or set linker.positivity = false.");
	
	// Linker sanity check
	//ensure(use_linker || write_noise || write_filtered || write_rawmask, ERR_USER_INPUT, "When disabling the linker, you will want to write either the\n       noise cube, the filtered cube or the raw mask, as otherwise\n       no output would be produced at all.");
	
	// Print parameterisation warning
	if(use_parameteriser && use_physical)
	{
		warning( "┌──────────────────────────────────────────────────────────┐\n"
		"         │ You have set parameter.physical = true. SoFiA will try   │\n"
		"         │ to convert some parameters to physical units under the   │\n"
		"         │ following fundamental assumptions:                       │\n"
		"         │                                                          │\n"
		"         │  * The beam information in the FITS header (BMAJ, BMIN)  │\n"
		"         │    is correct and accurate across the entire data cube.  │\n"
		"         │                                                          │\n"
		"         │  * The spectral channels of the data cube are uncorrela- │\n"
		"         │    ted, i.e. spectral resolution equals channel width.   │\n"
		"         │                                                          │\n"
		"         │ Should any of these assumptions be incorrect then the    │\n"
		"         │ measurement of parameters such as flux or line width may │\n"
		"         │ be wrong. Note that SoFiA will in principle not correct  │\n"
		"         │ line widths for any form of instrumental broadening.     │\n"
		"         └──────────────────────────────────────────────────────────┘\n");
	}
	
	
	
	// ---------------------------- //
	// Define file names and paths  //
	// ---------------------------- //
	
	const char *base_dir  = Parameter_get_str(par, "output.directory");
	const char *base_name = Parameter_get_str(par, "output.filename");
	
	// Create input paths
	Path *path_data_in    = Path_new();
	Path *path_gain_in    = Path_new();
	Path *path_noise_in   = Path_new();
	Path *path_weights_in = Path_new();
	Path *path_mask_in    = Path_new();
	Path *path_flag_cat   = Path_new();
	
	// Set up input paths
	                     Path_set(path_data_in,    Parameter_get_str(par, "input.data"));
	if(use_gain)         Path_set(path_gain_in,    Parameter_get_str(par, "input.gain"));
	if(use_noise)        Path_set(path_noise_in,   Parameter_get_str(par, "input.noise"));
	if(use_weights)      Path_set(path_weights_in, Parameter_get_str(par, "input.weights"));
	if(use_mask)         Path_set(path_mask_in,    Parameter_get_str(par, "input.mask"));
	if(use_flagging_cat) Path_set(path_flag_cat,   Parameter_get_str(par, "flag.catalog"));
	
	// Choose appropriate output file and directory names depending on user input
	String *output_file_name = String_new(strlen(base_name) ? base_name : Path_get_file(path_data_in));
	String *output_dir_name  = String_new(strlen(base_dir) ? base_dir : (strlen(Path_get_dir(path_data_in)) ? Path_get_dir(path_data_in) : "."));
	
	// Ensure that output file name ends with .fits or similar
	String *check_mime_type = String_new("");
	String_to_lower(String_set_delim(check_mime_type, String_get(output_file_name), '.', false, false));
	if(!String_compare(check_mime_type, "fits") && !String_compare(check_mime_type, "fit")) String_append(output_file_name, ".fits");
	String_delete(check_mime_type);
	
	// Create global output paths
	Path *path_cat_ascii = Path_new();
	Path *path_cat_xml   = Path_new();
	Path *path_cat_sql   = Path_new();
	Path *path_noise_out = Path_new();
	Path *path_filtered  = Path_new();
	Path *path_mask_out  = Path_new();
	Path *path_mask_2d   = Path_new();
	Path *path_mask_raw  = Path_new();
	Path *path_mom0      = Path_new();
	Path *path_mom1      = Path_new();
	Path *path_mom2      = Path_new();
	Path *path_chan      = Path_new();
	Path *path_rel_plot  = Path_new();
	Path *path_rel_cat_n = Path_new();
	Path *path_rel_cat_p = Path_new();
	Path *path_skel_plot = Path_new();
	Path *path_flag      = Path_new();
	Path *path_cubelets  = Path_new();
	
	// Set up global output directory names
	Path_set_dir(path_cat_ascii, String_get(output_dir_name));
	Path_set_dir(path_cat_xml,   String_get(output_dir_name));
	Path_set_dir(path_cat_sql,   String_get(output_dir_name));
	Path_set_dir(path_noise_out, String_get(output_dir_name));
	Path_set_dir(path_filtered,  String_get(output_dir_name));
	Path_set_dir(path_mask_out,  String_get(output_dir_name));
	Path_set_dir(path_mask_2d,   String_get(output_dir_name));
	Path_set_dir(path_mask_raw,  String_get(output_dir_name));
	Path_set_dir(path_mom0,      String_get(output_dir_name));
	Path_set_dir(path_mom1,      String_get(output_dir_name));
	Path_set_dir(path_mom2,      String_get(output_dir_name));
	Path_set_dir(path_chan,      String_get(output_dir_name));
	Path_set_dir(path_rel_plot,  String_get(output_dir_name));
	Path_set_dir(path_rel_cat_n, String_get(output_dir_name));
	Path_set_dir(path_rel_cat_p, String_get(output_dir_name));
	Path_set_dir(path_skel_plot, String_get(output_dir_name));
	Path_set_dir(path_flag,      String_get(output_dir_name));
	Path_set_dir(path_cubelets,  String_get(output_dir_name));
	
	// Set up global output file names
	Path_set_file_from_template(path_cat_ascii,  String_get(output_file_name), "_cat",         ".txt");
	Path_set_file_from_template(path_cat_xml,    String_get(output_file_name), "_cat",         ".xml");
	Path_set_file_from_template(path_cat_sql,    String_get(output_file_name), "_cat",         ".sql");
	Path_set_file_from_template(path_noise_out,  String_get(output_file_name), "_noise",       use_local_scaling ? ".fits" : ".txt");
	Path_set_file_from_template(path_filtered,   String_get(output_file_name), "_filtered",    ".fits");
	Path_set_file_from_template(path_mask_out,   String_get(output_file_name), "_mask",        ".fits");
	Path_set_file_from_template(path_mask_2d,    String_get(output_file_name), "_mask-2d",     ".fits");
	Path_set_file_from_template(path_mask_raw,   String_get(output_file_name), "_mask-raw",    ".fits");
	Path_set_file_from_template(path_mom0,       String_get(output_file_name), "_mom0",        ".fits");
	Path_set_file_from_template(path_mom1,       String_get(output_file_name), "_mom1",        ".fits");
	Path_set_file_from_template(path_mom2,       String_get(output_file_name), "_mom2",        ".fits");
	Path_set_file_from_template(path_chan,       String_get(output_file_name), "_chan",        ".fits");
	Path_set_file_from_template(path_rel_plot,   String_get(output_file_name), "_rel",         ".eps");
	Path_set_file_from_template(path_rel_cat_n,  String_get(output_file_name), "_rel_cat_neg", ".xml");
	Path_set_file_from_template(path_rel_cat_p,  String_get(output_file_name), "_rel_cat_pos", ".xml");
	Path_set_file_from_template(path_skel_plot,  String_get(output_file_name), "_skellam",     ".eps");
	Path_set_file_from_template(path_flag,       String_get(output_file_name), "_flags",       ".log");
	
	// Set up cubelet directory and file base name
	Path_append_dir_from_template(path_cubelets, String_get(output_file_name), "_cubelets");
	Path_set_file_from_template(path_cubelets,   String_get(output_file_name), "", "");
	
	// Delete temporary strings again
	String_delete(output_file_name);
	String_delete(output_dir_name);
	
	
	
	// ---------------------------- //
	// Check input files            //
	// ---------------------------- //
	
	ensure(Path_file_is_readable(path_data_in), ERR_FILE_ACCESS,
		"Failed to read data cube. Please ensure\n       that the file exists and is readable.");
	
	if(use_gain) {
		ensure(Path_file_is_readable(path_gain_in), ERR_FILE_ACCESS,
			"Failed to read gain cube. Please ensure\n       that the file exists and is readable.");
	}
	if(use_noise) {
		ensure(Path_file_is_readable(path_noise_in), ERR_FILE_ACCESS,
			"Failed to read noise cube. Please ensure\n       that the file exists and is readable.");
	}
	if(use_weights) {
		ensure(Path_file_is_readable(path_weights_in), ERR_FILE_ACCESS,
			"Failed to read weights cube. Please ensure\n       that the file exists and is readable.");
	}
	if(use_mask) {
		ensure(Path_file_is_readable(path_mask_in), ERR_FILE_ACCESS,
			"Failed to read mask cube. Please ensure\n       that the file exists and is readable.");
	}
	if(use_flagging_cat) {
		ensure(Path_file_is_readable(path_flag_cat), ERR_FILE_ACCESS,
			"Failed to read flagging catalogue. Please ensure\n       that the file exists and is readable.");
	}
	
	
	
	// ---------------------------- //
	// Check output settings        //
	// ---------------------------- //
	
	// Try to create cubelet directory
	if(write_cubelets)
	{
		errno = 0;
		mkdir(Path_get_dir(path_cubelets), 0755);
		ensure(errno == 0 || errno == EEXIST, ERR_FILE_ACCESS, "Failed to create cubelet directory; please check write permissions.");
	}
	
	// Check overwrite conditions
	if(!overwrite)
	{
		if(write_cubelets) {
			ensure(errno != EEXIST, ERR_FILE_ACCESS,
				"Cubelet directory already exists. Please delete the directory\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_ascii) {
			ensure(!Path_file_is_readable(path_cat_ascii), ERR_FILE_ACCESS,
				"ASCII catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_xml) {
			ensure(!Path_file_is_readable(path_cat_xml), ERR_FILE_ACCESS,
				"XML catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_sql) {
			ensure(!Path_file_is_readable(path_cat_sql), ERR_FILE_ACCESS,
				"SQL catalogue file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(use_noise_scaling && write_noise) {
			ensure(!Path_file_is_readable(path_noise_out), ERR_FILE_ACCESS,
				"Noise cube/spectrum already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_filtered) {
			ensure(!Path_file_is_readable(path_filtered), ERR_FILE_ACCESS,
				"Filtered cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_mask) {
			ensure(!Path_file_is_readable(path_mask_out), ERR_FILE_ACCESS,
				"Mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_mask2d) {
			ensure(!Path_file_is_readable(path_mask_2d), ERR_FILE_ACCESS,
				"2-D mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_rawmask) {
			ensure(!Path_file_is_readable(path_mask_raw), ERR_FILE_ACCESS,
				"Raw mask cube already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(write_moments) {
			ensure(!Path_file_is_readable(path_mom0) && !Path_file_is_readable(path_mom1) && !Path_file_is_readable(path_mom2), ERR_FILE_ACCESS,
				"Moment maps already exist. Please delete the files\n"
				"       or set \'output.overwrite = true\'.");
			ensure(!Path_file_is_readable(path_chan), ERR_FILE_ACCESS,
				"Channel map already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(use_reliability && use_rel_plot) {
			ensure(!Path_file_is_readable(path_rel_plot), ERR_FILE_ACCESS,
				"Reliability plot already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
			ensure(!Path_file_is_readable(path_skel_plot), ERR_FILE_ACCESS,
				"Skellam plot already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
		if(use_reliability && use_rel_debug) {
			ensure(!Path_file_is_readable(path_rel_cat_n) && !Path_file_is_readable(path_rel_cat_p), ERR_FILE_ACCESS,
				   "Reliability debugging catalogue already exists. Please delete\n"
				   "       the file or set \'output.overwrite = true\'.");
		}
		if(autoflag_log) {
			ensure(!Path_file_is_readable(path_flag), ERR_FILE_ACCESS,
				"Flagging log file already exists. Please delete the file\n"
				"       or set \'output.overwrite = true\'.");
		}
	}
	
	
	
	// ---------------------------- //
	// Load data cube               //
	// ---------------------------- //
	
	// Set up region if required
	Array_siz *region = use_region ? Array_siz_new_str(Parameter_get_str(par, "input.region")) : NULL;
	
	// Set up flagging region if required
	Array_siz *flag_regions = use_flagging ? Array_siz_new_str(Parameter_get_str(par, "flag.region")) : Array_siz_new(0);
	
	// Load data cube
	status("Loading data cube");
	DataCube *dataCube = DataCube_new(verbosity);
	DataCube_load(dataCube, Path_get(path_data_in), region);
	
	// Check for CELLSCAL = '1/F' setting
	if(DataCube_cmphd(dataCube, "CELLSCAL", "1/F", 3))
	{
		warning( "┌──────────────────────────────────────────────────────────┐\n"
		"         │ FITS header keyword of CELLSCAL = '1/F' detected.  SoFiA │\n"
		"         │ does not support cubes with variable spatial pixel size. │\n"
		"         │ While basic source finding will work, there could be un- │\n"
		"         │ foreseen side effects, including wrong coordinates in 2D │\n"
		"         │ output images from SoFiA (such as moment maps) and inac- │\n"
		"         │ curate celestial coordinate measurements.                │\n"
                "         │ It is possible to avoid this problem by regridding the   │\n"
                "         │ cube to a fixed pixel size before running SoFiA. This    │\n"
                "         │ can be done with, e.g., MIRIAD's task 'regrid' setting   │\n"
                "         │ 'options=noscale'. However, note that the regridded cube │\n"
                "         │ may have a varying beam size (in pixels), which may lead │\n"
                "         │ to errors in, e.g., the measurement of total fluxes.     │\n"
		"         └──────────────────────────────────────────────────────────┘\n");
	}
	
	// Search for values of infinity and append affected pixels to flagging region
	// (Yes, some data cubes do contain those!)
	if(DataCube_flag_infinity(dataCube, flag_regions)) use_flagging = true;
	
	// Apply flags if required
	if(use_flagging) DataCube_flag_regions(dataCube, flag_regions);
	
	// Invert cube if requested
	if(use_invert)
	{
		message("Inverting data cube");
		DataCube_multiply_const(dataCube, -1.0);
	}
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Apply flagging catalogue     //
	// ---------------------------- //
	
	if(use_flagging_cat)
	{
		status("Loading and applying flagging catalogue");
		message("Catalogue file:   %s", Parameter_get_str(par, "flag.catalog"));
		message("Flagging radius:  %ld", Parameter_get_int(par, "flag.radius"));
		DataCube_continuum_flagging(dataCube, Parameter_get_str(par, "flag.catalog"), 1, Parameter_get_int(par, "flag.radius"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load and apply noise cube    //
	// ---------------------------- //
	
	if(use_noise)
	{
		status("Loading and applying noise cube");
		DataCube *noiseCube = DataCube_new(verbosity);
		DataCube_load(noiseCube, Path_get(path_noise_in), region);
		
		// Divide data by noise cube
		DataCube_divide(dataCube, noiseCube);
		
		// Delete noise cube again
		DataCube_delete(noiseCube);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load and apply weights cube  //
	// ---------------------------- //
	
	if(use_weights)
	{
		status("Loading and applying weights cube");
		DataCube *weightsCube = DataCube_new(verbosity);
		DataCube_load(weightsCube, Path_get(path_weights_in), region);
		
		// Multiply data by square root of weights cube
		DataCube_apply_weights(dataCube, weightsCube);
		
		// Delete weights cube again
		DataCube_delete(weightsCube);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Continuum subtraction        //
	// ---------------------------- //
	
	if(use_cont_sub)
	{
		status("Continuum subtraction");
		message("Subtracting residual continuum emission.");
		message("- Polynomial order:  %ld",   Parameter_get_int(par, "contsub.order"));
		message("- Clip threshold:    %.1f",  Parameter_get_flt(par, "contsub.threshold"));
		message("- Shift:             %ld",   Parameter_get_int(par, "contsub.shift"));
		message("- Padding:           %ld\n", Parameter_get_int(par, "contsub.padding"));
		
		DataCube_contsub(dataCube, Parameter_get_int(par, "contsub.order"), Parameter_get_int(par, "contsub.shift"), Parameter_get_int(par, "contsub.padding"), Parameter_get_flt(par, "contsub.threshold"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Scale data by noise level    //
	// ---------------------------- //
	
	if(use_noise_scaling)
	{
		status("Scaling data by noise");
		
		if(use_local_scaling)
		{
			// Local noise scaling
			message("Correcting for local noise variations.");
			message("- Noise statistic:  %s", noise_stat_name[sn_statistic]);
			message("- Flux range:       %s\n", flux_range_name[sn_range + 1]);
			
			DataCube *noiseCube = DataCube_scale_noise_local(
				dataCube,
				sn_statistic,
				sn_range,
				Parameter_get_int(par, "scaleNoise.windowXY"),
				Parameter_get_int(par, "scaleNoise.windowZ"),
				Parameter_get_int(par, "scaleNoise.gridXY"),
				Parameter_get_int(par, "scaleNoise.gridZ"),
				Parameter_get_bool(par, "scaleNoise.interpolate")
			);
			
			if(write_noise)
			{
				// Apply flags to noise cube
				if(use_flagging) DataCube_flag_regions(noiseCube, flag_regions);
				DataCube_add_history(noiseCube, par);
				DataCube_save(noiseCube, Path_get(path_noise_out), overwrite, DESTROY);
			}
			DataCube_delete(noiseCube);
		}
		else
		{
			// Global noise scaling along spectral axis
			message("Correcting for noise variations along spectral axis.");
			message("- Noise statistic:  %s", noise_stat_name[sn_statistic]);
			message("- Flux range:       %s\n", flux_range_name[sn_range + 1]);
			Array_dbl *noise_spectrum = DataCube_scale_noise_spec(dataCube, sn_statistic, sn_range);
			
			if(write_noise)
			{
				FILE *fp;
				if(overwrite) fp = fopen(Path_get(path_noise_out), "wb");
				else fp = fopen(Path_get(path_noise_out), "wxb");
				ensure(fp != NULL, ERR_FILE_ACCESS, "Failed to open output file: %s", Path_get(path_noise_out));
				message("Writing noise spectrum: %s", Path_get_file(path_noise_out));
				
				fprintf(fp, "# Measured noise spectrum\n# Creator: %s\n#\n# Channel\tNoise\n", SOFIA_VERSION_FULL);
				for(size_t i = 0; i < Array_dbl_get_size(noise_spectrum); ++i) fprintf(fp, "%zu\t%.15f\n", i, Array_dbl_get(noise_spectrum, i));
				
				fclose(fp);
			}
			
			Array_dbl_delete(noise_spectrum);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Automatic data flagging      //
	// ---------------------------- //
	
	if(autoflag_mode)
	{
		status("Auto-flagging");
		
		// Set up auto-flagging if requested
		Array_siz *autoflag_regions = Array_siz_new(0);
		DataCube_autoflag(dataCube, Parameter_get_flt(par, "flag.threshold"), autoflag_mode, autoflag_regions);
		
		const size_t size = Array_siz_get_size(autoflag_regions);
		
		// Apply flags if necessary
		if(size)
		{
			DataCube_flag_regions(dataCube, autoflag_regions);  // Apply auto-flagging regions
			Array_siz_cat(flag_regions, autoflag_regions);      // Append auto-flagging regions to general flagging regions
			use_flagging = true;                                // Update flagging switch
		}
		else message("No flagging required.");
		
		// Write auto-flags to log file if requested
		if(size && autoflag_log)
		{
			// Try to open output file
			FILE *fp;
			if(overwrite) fp = fopen(Path_get(path_flag), "wb");
			else fp = fopen(Path_get(path_flag), "wxb");
			
			// If successful...
			if(fp != NULL)
			{	
				// ...write out flags...
				message("Writing log file:     %s", Path_get_file(path_flag));
				fprintf(fp, "# Auto-flagging log file\n");
				fprintf(fp, "# Creator: %s\n#\n", SOFIA_VERSION_FULL);
				fprintf(fp, "# Flagging codes:\n");
				fprintf(fp, "#   C z            =  spectral channel (z)\n");
				fprintf(fp, "#   P x y          =  spatial pixel (x, y)\n");
				fprintf(fp, "#   R x1 x2 y1 y2  =  spatial region (x1:x2, y1:y2)\n");
				fprintf(fp, "# Note that coordinates will be relative to subregion\n");
				fprintf(fp, "# unless parameter.offset was set to true.\n\n");
				
				for(size_t i = 0; i < size; i += 6)
				{
					const size_t x_min = Array_siz_get(autoflag_regions, i);
					const size_t x_max = Array_siz_get(autoflag_regions, i + 1);
					const size_t y_min = Array_siz_get(autoflag_regions, i + 2);
					const size_t y_max = Array_siz_get(autoflag_regions, i + 3);
					const size_t z_min = Array_siz_get(autoflag_regions, i + 4);
					const size_t z_max = Array_siz_get(autoflag_regions, i + 5);
					
					// NOTE: Subregion offset will be added if requested (use_pos_offset == true).
					if(z_min == z_max)
					{
						fprintf(fp, "C %zu\n", z_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 4) : 0));
					}
					else if(x_min == x_max && y_min == y_max)
					{
						fprintf(fp, "P %zu %zu\n", x_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), y_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0));
					}
					else if(z_min == 0 && z_max == DataCube_get_axis_size(dataCube, 2) - 1)
					{
						fprintf(fp, "R %zu %zu %zu %zu\n", x_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), x_max + ((use_region && use_pos_offset) ? Array_siz_get(region, 0) : 0), y_min + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0), y_max + ((use_region && use_pos_offset) ? Array_siz_get(region, 2) : 0));
					}
				}
				
				// ...and close output file again
				fclose(fp);
			}
			else warning("Failed to write flagging log file: %s", Path_get_file(path_flag));
		}
		
		// Clean up
		Array_siz_delete(autoflag_regions);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Ripple filter                //
	// ---------------------------- //
	
	if(use_ripple_filter)
	{
		status("Applying ripple filter");
		
		// Ripple filter
		message("Subtracting offset from data.");
		message("- Statistic:    %s", noise_stat_name[ripple_filter_statistic]);
		
		DataCube *rippleCube = DataCube_scale_noise_local(
			dataCube,
			ripple_filter_statistic,
			0,
			Parameter_get_int(par, "rippleFilter.windowXY"),
			Parameter_get_int(par, "rippleFilter.windowZ"),
			Parameter_get_int(par, "rippleFilter.gridXY"),
			Parameter_get_int(par, "rippleFilter.gridZ"),
			Parameter_get_bool(par, "rippleFilter.interpolate")
		);
		
		DataCube_delete(rippleCube);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Write filtered cube          //
	// ---------------------------- //
	
	if(write_filtered && (use_region || use_flagging || use_flagging_cat || use_cont_sub || use_noise || use_weights || use_noise_scaling || use_ripple_filter))  // ALERT: Add conditions here as needed.
	{
		status("Writing filtered cube");
		DataCube_add_history(dataCube, par);
		DataCube_save(dataCube, Path_get(path_filtered), overwrite, PRESERVE);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Measure global noise level   //
	// ---------------------------- //
	
	// NOTE: This is necessary so the linker and reliability module can
	//       divide all flux values by the RMS later on.
	// NOTE: This is currently being applied even when a noise cube has 
	//       been applied before or noise scaling is enabled.
	//       This is needed, as other algorithms, such as continuum sub-
	//       traction, might alter the noise level.
	
	status("Measuring global noise level");
	
	size_t cadence = DataCube_get_size(dataCube) / NOISE_SAMPLE_SIZE;          // Stride for noise calculation
	if(cadence < 2) cadence = 1;
	else if(cadence % DataCube_get_axis_size(dataCube, 0) == 0) cadence -= 1;  // Ensure stride is not equal to multiple of x-axis size
	
	global_rms = MAD_TO_STD * DataCube_stat_mad(dataCube, 0.0, cadence, -1);
	message("Global RMS:  %.3e  (using stride of %zu)", global_rms, cadence);
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Run source finder            //
	// ---------------------------- //
	
	// Terminate if no source finder is to be run, but no input mask is provided either
	ensure(use_scfind || use_threshold || use_mask, ERR_USER_INPUT, "No mask provided and no source finder selected. Cannot proceed.");
	
	// Create temporary 8-bit mask to hold source finding output
	DataCube *maskCubeTmp = DataCube_blank(DataCube_get_axis_size(dataCube, 0), DataCube_get_axis_size(dataCube, 1), DataCube_get_axis_size(dataCube, 2), 8, verbosity);
	DataCube_copy_wcs(dataCube, maskCubeTmp);
	DataCube_puthd_str(maskCubeTmp, "BUNIT", " ");
	
	// S+C finder
	if(use_scfind)
	{
		status("Running S+C finder");
		message("Using the following parameters:");
		message("- Kernels");
		message("  - spatial:        %s", Parameter_get_str(par, "scfind.kernelsXY"));
		message("  - spectral:       %s", Parameter_get_str(par, "scfind.kernelsZ"));
		message("- Flux threshold:   %s * rms", Parameter_get_str(par, "scfind.threshold"));
		message("- Noise statistic:  %s", noise_stat_name[sc_statistic]);
		message("- Flux range:       %s\n", flux_range_name[sc_range + 1]);
		
		// Extract and sort kernel sizes to ensure that smallest kernel comes first
		Array_dbl *kernels_spat = Array_dbl_new_str(Parameter_get_str(par, "scfind.kernelsXY"));
		Array_siz *kernels_spec = Array_siz_new_str(Parameter_get_str(par, "scfind.kernelsZ"));
		Array_dbl_sort(kernels_spat);
		Array_siz_sort(kernels_spec);
		
		// Sanity checks
		for(size_t i = 0; i < Array_dbl_get_size(kernels_spat); ++i)
		{
			const double ks = Array_dbl_get(kernels_spat, i);
			ensure(ks >= 0.0 && ks < DataCube_get_axis_size(dataCube, 0) && ks < DataCube_get_axis_size(dataCube, 1), ERR_USER_INPUT, "Illegal spatial kernel size encountered.");
			if(ks > 0.0 && ks < 3.0) warning("Spatial kernel sizes of < 3 cannot be accurately modelled.");
		}
		for(size_t i = 0; i < Array_siz_get_size(kernels_spec); ++i)
		{
			const size_t ks = Array_siz_get(kernels_spec, i);
			ensure(ks < DataCube_get_axis_size(dataCube, 2), ERR_USER_INPUT, "Illegal spectral kernel size encountered.");
			if(ks != 0 && ks % 2 == 0) warning("Spectral kernel size of %zu is even, will be treated as %zu!", ks, ks + 1);
			else if(ks == 1) warning("Spectral kernel size of 1 found, will be treated as 0!");
		}
		if(Array_dbl_get(kernels_spat, 0) > 0.0) warning("Including spatial kernel size of 0 is strongly advised.");
		if(Array_siz_get(kernels_spec, 0) > 0) warning("Including spectral kernel size of 0 is strongly advised.");
		
		// Run S+C finder to obtain mask
		DataCube_run_scfind(
			dataCube,
			maskCubeTmp,
			kernels_spat,
			kernels_spec,
			Parameter_get_flt(par, "scfind.threshold"),
			Parameter_get_flt(par, "scfind.replacement"),
			sc_statistic,
			sc_range,
			(use_noise_scaling && use_sc_scaling) ? (strcmp(Parameter_get_str(par, "scaleNoise.mode"), "local") == 0 ? 2 : 1) : 0,
			sn_statistic,
			sn_range,
			Parameter_get_int(par, "scaleNoise.windowXY"),
			Parameter_get_int(par, "scaleNoise.windowZ"),
			Parameter_get_int(par, "scaleNoise.gridXY"),
			Parameter_get_int(par, "scaleNoise.gridZ"),
			Parameter_get_bool(par, "scaleNoise.interpolate"),
			start_time,
			start_clock
		);
		
		// Clean up
		Array_dbl_delete(kernels_spat);
		Array_siz_delete(kernels_spec);
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCubeTmp, flag_regions);
	}
	
	// Threshold finder
	if(use_threshold)
	{
		// Determine mode
		const bool absolute = (strcmp(Parameter_get_str(par, "threshold.mode"), "absolute") == 0);
		
		status("Running threshold finder");
		message("Using the following parameters:");
		message("- Mode:             %s", absolute ? "absolute" : "relative");
		message("- Flux threshold:   %s%s", Parameter_get_str(par, "threshold.threshold"), absolute ? "" : " * rms");
		if(!absolute)
		{
			message("- Noise statistic:  %s", noise_stat_name[tf_statistic]);
			message("- Flux range:       %s", flux_range_name[tf_range + 1]);
		}
		
		// Run threshold finder
		DataCube_run_threshold(
			dataCube,
			maskCubeTmp,
			absolute,
			Parameter_get_flt(par, "threshold.threshold"),
			tf_statistic,
			tf_range
		);
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCubeTmp, flag_regions);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Load mask cube if specified  //
	// ---------------------------- //
	
	DataCube *maskCube = NULL;
	
	if(use_mask)
	{
		// Load mask cube
		status("Loading mask cube");
		DataCube *inputMaskCube = DataCube_new(verbosity);
		DataCube_load(inputMaskCube, Path_get(path_mask_in), region);
		
		// Ensure that mask has the right type and size
		ensure(
			DataCube_gethd_int(inputMaskCube, "BITPIX") == 8 ||
			DataCube_gethd_int(inputMaskCube, "BITPIX") == 16 ||
			DataCube_gethd_int(inputMaskCube, "BITPIX") == 32 ||
			DataCube_gethd_int(inputMaskCube, "BITPIX") == 64,
			ERR_USER_INPUT, "Mask cube must be of integer type."
		);
		ensure(
			DataCube_gethd_int(inputMaskCube, "NAXIS1") == DataCube_gethd_int(dataCube, "NAXIS1") &&
			DataCube_gethd_int(inputMaskCube, "NAXIS2") == DataCube_gethd_int(dataCube, "NAXIS2") &&
			DataCube_gethd_int(inputMaskCube, "NAXIS3") == DataCube_gethd_int(dataCube, "NAXIS3"),
			ERR_USER_INPUT, "Data cube and mask cube have different sizes."
		);
		
		if(DataCube_gethd_int(inputMaskCube, "BITPIX") == 32)
		{
			// If 32-bit, then make this the new source mask
			maskCube = inputMaskCube;
			
			// Set all masked pixels to -1
			DataCube_reset_mask_32(maskCube, -1);
			
			// NOTE: This will accept whatever WCS is defined in the input mask
			//       without any sanity checks! This should not be critical,
			//       though, as the WCS information from the mask is presumably
			//       not used anywhere.
		}
		else
		{
			// If 8, 16 or 64 bit, then create new mask and copy pixels across
			maskCube = DataCube_blank(DataCube_get_axis_size(dataCube, 0), DataCube_get_axis_size(dataCube, 1), DataCube_get_axis_size(dataCube, 2), 32, verbosity);
			
			// Copy WCS header elements from data cube to mask cube
			DataCube_copy_wcs(dataCube, maskCube);
			
			// Set BUNIT keyword of mask cube
			DataCube_puthd_str(maskCube, "BUNIT", " ");
			
			// Copy pixels across
			const size_t n_pix_det = DataCube_copy_mask_32(maskCube, inputMaskCube, -1);
			message("%zu pixels copied from input mask (%.3f%%).", n_pix_det, 100.0 * (double)(n_pix_det) / (double)(DataCube_get_size(inputMaskCube)));
			
			// Delete input mask again
			DataCube_delete(inputMaskCube);
		}
		
		// Apply flags to mask cube
		if(use_flagging) DataCube_flag_regions(maskCube, flag_regions);
	}
	else
	{
		// Else create an empty mask cube
		status("Creating source mask cube");
		maskCube = DataCube_blank(DataCube_get_axis_size(dataCube, 0), DataCube_get_axis_size(dataCube, 1), DataCube_get_axis_size(dataCube, 2), 32, verbosity);
		
		// Copy WCS header elements from data cube to mask cube
		DataCube_copy_wcs(dataCube, maskCube);
		
		// Set BUNIT keyword of mask cube
		DataCube_puthd_str(maskCube, "BUNIT", " ");
	}
	
	
	
	// ---------------------------- //
	// Copy SF mask across          //
	// ---------------------------- //
	
	// Copy SF mask prior to linking
	const size_t n_pix_det = DataCube_copy_mask_32(maskCube, maskCubeTmp, -1);
	message("%zu pixels detected by source finder (%.3f%%).", n_pix_det, 100.0 * (double)(n_pix_det) / (double)(DataCube_get_size(maskCube)));
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Write SF mask if requested   //
	// ---------------------------- //
	
	if(write_rawmask)
	{
		status("Writing raw binary mask");
		DataCube_add_history(maskCubeTmp, par);
		DataCube_save(maskCubeTmp, Path_get(path_mask_raw), overwrite, DESTROY);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	// Delete temporary SF mask again
	DataCube_delete(maskCubeTmp);
	
	
	
	// ---------------------------- //
	// Run linker                   //
	// ---------------------------- //
	
	ensure(use_linker, ERR_NO_SRC_FOUND, "Terminating pipeline, as linker is disabled.\n       No source catalogue has been created.");
	
	status("Running Linker");
	
	const bool remove_neg_src = !use_reliability && !keep_negative;  // ALERT: Add conditions here as needed.
	
	LinkerPar *lpar = DataCube_run_linker(
		dataCube,
		maskCube,
		Parameter_get_int(par, "linker.radiusXY"),
		Parameter_get_int(par, "linker.radiusXY"),
		Parameter_get_int(par, "linker.radiusZ"),
		Parameter_get_int(par, "linker.minSizeXY"),
		Parameter_get_int(par, "linker.minSizeXY"),
		Parameter_get_int(par, "linker.minSizeZ"),
		Parameter_get_int(par, "linker.minPixels"),
		Parameter_get_flt(par, "linker.minFill"),
		Parameter_get_int(par, "linker.maxSizeXY"),
		Parameter_get_int(par, "linker.maxSizeXY"),
		Parameter_get_int(par, "linker.maxSizeZ"),
		Parameter_get_int(par, "linker.maxPixels"),
		Parameter_get_flt(par, "linker.maxFill"),
		linker_pos_pix,
		remove_neg_src,
		global_rms
	);
	
	// Print time
	timestamp(start_time, start_clock);
	
	// Terminate pipeline if no sources left after linking
	ensure(LinkerPar_get_size(lpar), ERR_NO_SRC_FOUND, "No sources left after linking. Terminating pipeline.");
	
	
	
	// ---------------------------- //
	// Run reliability filter       //
	// ---------------------------- //
	
	Map *rel_filter = Map_new();  // Empty container for storing old and new labels of reliable sources
	
	if(use_reliability)
	{
		status("Measuring reliability");
		
		// Extract parameter space and dimensionality
		Array_siz *rel_par_space = Array_siz_new(0);
		char *rel_par_str = (char *)Parameter_get_str(par, "reliability.parameters");  // (strtok() doesn't like 'const')
		const char *token = strtok(rel_par_str, ", ");
		
		while(token != NULL)
		{
			if     (strcmp(token, "peak") == 0) Array_siz_push(rel_par_space, LINKERPAR_PEAK);
			else if(strcmp(token, "sum")  == 0) Array_siz_push(rel_par_space, LINKERPAR_SUM);
			else if(strcmp(token, "mean") == 0) Array_siz_push(rel_par_space, LINKERPAR_MEAN);
			else if(strcmp(token, "chan") == 0) Array_siz_push(rel_par_space, LINKERPAR_CHAN);
			else if(strcmp(token, "pix")  == 0) Array_siz_push(rel_par_space, LINKERPAR_PIX);
			else if(strcmp(token, "fill") == 0) Array_siz_push(rel_par_space, LINKERPAR_FILL);
			else if(strcmp(token, "std")  == 0) Array_siz_push(rel_par_space, LINKERPAR_STD);
			else if(strcmp(token, "skew") == 0) Array_siz_push(rel_par_space, LINKERPAR_SKEW);
			else if(strcmp(token, "kurt") == 0) Array_siz_push(rel_par_space, LINKERPAR_KURT);
			else ensure(false, ERR_USER_INPUT, "Unknown reliability parameter: '%s'.", token);
			token = strtok(NULL, ", ");
		}
		
		ensure(Array_siz_get_size(rel_par_space) > 1, ERR_USER_INPUT, "Reliability parameter space has < 2 dimensions.");
		
		// Check if catalogue supplied
		Table *rel_cat = NULL;
		if(use_rel_cat)
		{
			message("Reading in reliability catalogue.");
			
			// Read catalogue into table
			rel_cat = Table_from_file(Parameter_get_str(par, "reliability.catalog"), " \t,|");
			
			if(Table_rows(rel_cat) == 0 || Table_cols(rel_cat) != 2)
			{
				warning("Reliability catalogue non-compliant; must contain 2 data columns.\n         Catalogue file will be ignored.");
				Table_delete(rel_cat);
				rel_cat = NULL;
			}
			else
			{
				message("Extracting %zu position%s from catalogue.", Table_rows(rel_cat), Table_rows(rel_cat) > 1 ? "s" : "");
				
				// Extract WCS information
				WCS *wcs = DataCube_extract_wcs(dataCube);
				
				if(wcs != NULL)
				{
					// Loop over all rows and convert WCS to pixels
					for(size_t row = 0; row < Table_rows(rel_cat); ++row)
					{
						double lon = -1e+30;
						double lat = -1e+30;
						WCS_convertToPixel(wcs, Table_get(rel_cat, row, 0), Table_get(rel_cat, row, 1), 0.0, &lon, &lat, NULL);
						Table_set(rel_cat, row, 0, lon);
						Table_set(rel_cat, row, 1, lat);
					}
				}
				else
				{
					warning("WCS conversion failed; cannot apply reliability catalogue.");
					Table_delete(rel_cat);
					rel_cat = NULL;
				}
				
				// Delete WCS object again
				WCS_delete(wcs);
			}
		}
		
		// Extract beam information from header
		// WARNING: Implicitly assumes the same unit for BMAJ, BMIN and CDELT1!
		const double bmaj = DataCube_gethd_flt(dataCube, "BMAJ");
		const double bmin = DataCube_gethd_flt(dataCube, "BMIN");
		const double cdelt = fabs(DataCube_gethd_flt(dataCube, "CDELT1"));
		double sqrt_beam_area = 1.0;
		if(IS_NAN(bmaj) || IS_NAN(bmin) || IS_NAN(cdelt)) warning("Failed to determine beam area from header information.\n         Value of reliability.minSNR will not be meaningful.");
		else sqrt_beam_area = sqrt(M_PI * bmaj * bmin / (4.0 * log(2.0) * cdelt * cdelt));
		const double rel_fmin = rel_snr_min * sqrt_beam_area;
		
		// Calculate reliability values
		double scale_kernel = Parameter_get_flt(par, "reliability.scaleKernel");
		Array_dbl *skellam = NULL;
		Matrix *covar = LinkerPar_reliability(lpar, rel_par_space, &scale_kernel, rel_fmin, rel_min_pix, rel_cat, use_rel_plot || scale_kernel == 0 ? &skellam : NULL, Parameter_get_bool(par, "reliability.autoKernel"), Parameter_get_int(par, "reliability.iterations"), Parameter_get_flt(par, "reliability.tolerance"));
		
		// Create plots if requested
		if(use_rel_plot)
		{
			LinkerPar_rel_plots(lpar, rel_par_space, rel_threshold, rel_fmin, rel_snr_min, covar, Path_get(path_rel_plot), overwrite);
			LinkerPar_skellam_plot(skellam, Path_get(path_skel_plot), overwrite, scale_kernel);
		}
		
		// Clean up
		Array_dbl_delete(skellam);
		Array_siz_delete(rel_par_space);
		Matrix_delete(covar);
		Table_delete(rel_cat);
		
		// Set up relabelling filter by recording old and new label pairs of reliable sources
		size_t new_label = 1;
		
		for(size_t i = 0; i < LinkerPar_get_size(lpar); ++i)
		{
			const size_t old_label = LinkerPar_get_label(lpar, i);
			
			// Keep source if reliability > threshold and fmin = snr_min * sqrt(beam) parameter satisfied
			if(LinkerPar_get_rel(lpar, old_label) >= rel_threshold && LinkerPar_get_flux(lpar, old_label) / sqrt(LinkerPar_get_npix(lpar, old_label)) > rel_fmin) Map_push(rel_filter, old_label, new_label++);
		}
		
		// Check if any reliable sources left
		ensure(Map_get_size(rel_filter), ERR_NO_SRC_FOUND, "No reliable sources found. Terminating pipeline.");
		message("%zu reliable %s found.", Map_get_size(rel_filter), Map_get_size(rel_filter) > 1 ? "sources" : "source");
		
		// Apply filter to mask cube, so unreliable sources are removed
		// and reliable ones relabelled in consecutive order
		DataCube_filter_mask_32(maskCube, rel_filter);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create initial catalogue     //
	// ---------------------------- //
	
	status("Creating initial catalogue");
	
	// Extract flux unit from header
	String *unit_flux = String_trim(DataCube_gethd_string(dataCube, "BUNIT"));
	if(!String_size(unit_flux))
	{
		warning("No flux unit (\'BUNIT\') defined in header.");
		String_set(unit_flux, "???");
	}
	
	// Generate catalogue of reliable sources from linker output
	Catalog *catalog = LinkerPar_make_catalog(lpar, rel_filter, String_get(unit_flux));
	message("Initial source catalogue created.");
	
	// Create and save reliability parameter catalogues if requested
	if(use_reliability && use_rel_debug)
	{
		message("Writing reliability debugging catalogues:");
		message(" - %s", Path_get_file(path_rel_cat_n));
		message(" - %s", Path_get_file(path_rel_cat_p));
		
		Catalog *cat_rel_par_neg = Catalog_new();
		Catalog *cat_rel_par_pos = Catalog_new();
		
		LinkerPar_get_rel_cat(lpar, String_get(unit_flux), &cat_rel_par_neg, &cat_rel_par_pos);
		Catalog_save(cat_rel_par_neg, Path_get(path_rel_cat_n), CATALOG_FORMAT_XML, overwrite, par);
		Catalog_save(cat_rel_par_pos, Path_get(path_rel_cat_p), CATALOG_FORMAT_XML, overwrite, par);
		
		Catalog_delete(cat_rel_par_neg);
		Catalog_delete(cat_rel_par_pos);
	}
	
	// Delete linker parameters, reliability filter and flux unit string, as they are no longer needed
	LinkerPar_delete(lpar);
	Map_delete(rel_filter);
	String_delete(unit_flux);
	
	// Terminate if catalogue is empty
	ensure(Catalog_get_size(catalog), ERR_NO_SRC_FOUND, "No reliable sources found. Terminating pipeline.");
	
	// Print time
	timestamp(start_time, start_clock);
	
	
	
	// ---------------------------- //
	// Mask dilation if requested   //
	// ---------------------------- //
	// NOTE: It is not yet clear if mask dilation should happen in the noise-normalised data cube
	//       or the original data cube. Some more though will need to go into this...
	
	if(use_mask_dilation)
	{
		status("Mask dilation");
		
		message("Spectral dilation");
		DataCube_dilate_mask_z(dataCube, maskCube, catalog, Parameter_get_int(par, "dilation.iterationsZ"), Parameter_get_flt(par, "dilation.threshold"));
		
		message("Spatial dilation");
		DataCube_dilate_mask_xy(dataCube, maskCube, catalog, Parameter_get_int(par, "dilation.iterationsXY"), Parameter_get_flt(par, "dilation.threshold"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	// ---------------------------- //
	// Reload data cube if required //
	// ---------------------------- //
	
	if(use_noise || use_weights || use_noise_scaling)  // ALERT: Add conditions here as needed.
	{
		status("Reloading data cube for parameterisation");
		DataCube_load(dataCube, Path_get(path_data_in), region);
		
		// Apply flags if required
		if(use_flagging) DataCube_flag_regions(dataCube, flag_regions);
		
		// Apply flagging catalogue if required
		if(use_flagging_cat) DataCube_continuum_flagging(dataCube, Parameter_get_str(par, "flag.catalog"), 1, Parameter_get_int(par, "flag.radius"));
		
		// Invert cube if requested
		if(use_invert)
		{
			message("Inverting data cube");
			DataCube_multiply_const(dataCube, -1.0);
		}
		
		// Apply gain cube if provided
		if(use_gain)
		{
			status("Loading and applying gain cube");
			DataCube *gainCube = DataCube_new(verbosity);
			DataCube_load(gainCube, Path_get(path_gain_in), region);
			
			// Divide by gain cube
			DataCube_divide(dataCube, gainCube);
			
			// Delete gain cube again
			DataCube_delete(gainCube);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Parameterise sources         //
	// ---------------------------- //
	
	if(use_parameteriser)
	{
		status("Measuring source parameters");
		DataCube_parameterise(dataCube, maskCube, catalog, use_wcs, use_physical, Parameter_get_str(par, "parameter.prefix"));
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create and save cubelets     //
	// ---------------------------- //
	
	if(write_cubelets)
	{
		status("Creating cubelets");
		message("Flux threshold (moment 1 and 2): %.2e", thresh_mom);
		DataCube_create_cubelets(
			dataCube,
			maskCube,
			catalog,
			Path_get(path_cubelets),
			overwrite,
			use_wcs,
			use_physical,
			Parameter_get_int(par, "output.marginCubelets"),
			thresh_mom,
			(use_region && use_pos_offset) ? Array_siz_get(region, 4) : 0,
			par
		);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Create and save moment maps  //
	// ---------------------------- //
	
	if(write_moments)
	{
		status("Creating moment maps");
		
		// Generate moment maps
		DataCube *mom0 = NULL;
		DataCube *mom1 = NULL;
		DataCube *mom2 = NULL;
		DataCube *chan = NULL;
		DataCube *snr  = NULL;
		DataCube_create_moments(dataCube, maskCube, &mom0, &mom1, &mom2, &chan, &snr, NULL, use_wcs, 0.0, 0.0);
		// NOTE: snr will not actually be created or written, as this would not work in general
		//       unless the noise was guaranteed to be constant across the entire data cube.
		
		// Save moment maps to disk
		if(mom0 != NULL)
		{
			DataCube_add_history(mom0, par);
			DataCube_save(mom0, Path_get(path_mom0), overwrite, DESTROY);
		}
		if(mom1 != NULL)
		{
			DataCube_add_history(mom1, par);
			DataCube_save(mom1, Path_get(path_mom1), overwrite, DESTROY);
		}
		if(mom2 != NULL)
		{
			DataCube_add_history(mom2, par);
			DataCube_save(mom2, Path_get(path_mom2), overwrite, DESTROY);
		}
		if(chan != NULL)
		{
			DataCube_add_history(chan, par);
			DataCube_save(chan, Path_get(path_chan), overwrite, DESTROY);
		}
		
		// Delete moment maps again
		DataCube_delete(mom0);
		DataCube_delete(mom1);
		DataCube_delete(mom2);
		DataCube_delete(chan);
		DataCube_delete(snr);
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Save mask cube               //
	// ---------------------------- //
	
	if(write_mask || write_mask2d)
	{
		status("Writing mask cube");
		
		// Create and save projected 2-D mask image
		if(write_mask2d)
		{
			DataCube *maskImage = DataCube_2d_mask(maskCube);
			DataCube_add_history(maskImage, par);
			DataCube_save(maskImage, Path_get(path_mask_2d), overwrite, DESTROY);
			DataCube_delete(maskImage);
		}
		
		// Write 3-D mask cube
		if(write_mask)
		{
			DataCube_add_history(maskCube, par);
			DataCube_save(maskCube, Path_get(path_mask_out), overwrite, DESTROY);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Save catalogue(s)            //
	// ---------------------------- //
	
	if(write_ascii || write_xml || write_sql)
	{
		status("Writing source catalogue");
		
		// Correct x, y and z for subregion offset if requested
		// WARNING: This will alter the original x, y and z positions!
		if(use_region && use_pos_offset)
		{
			for(size_t i = Catalog_get_size(catalog); i--;)
			{
				Source *src = Catalog_get_source(catalog, i);
				Source_offset_xyz(src, Array_siz_get(region, 0), Array_siz_get(region, 2), Array_siz_get(region, 4));
			}
		}
		
		if(write_ascii)
		{
			message("Writing ASCII file:   %s", Path_get_file(path_cat_ascii));
			Catalog_save(catalog, Path_get(path_cat_ascii), CATALOG_FORMAT_ASCII, overwrite, NULL);
		}
		
		if(write_xml)
		{
			message("Writing VOTable file: %s", Path_get_file(path_cat_xml));
			Catalog_save(catalog, Path_get(path_cat_xml), CATALOG_FORMAT_XML, overwrite, par);
		}
		
		if(write_sql)
		{
			message("Writing SQL file:     %s", Path_get_file(path_cat_sql));
			Catalog_save(catalog, Path_get(path_cat_sql), CATALOG_FORMAT_SQL, overwrite, NULL);
		}
		
		// Print time
		timestamp(start_time, start_clock);
	}
	
	
	
	// ---------------------------- //
	// Clean up and exit            //
	// ---------------------------- //
	
	// Delete data cube and mask cube
	DataCube_delete(maskCube);
	DataCube_delete(dataCube);
	
	// Delete sub-cube region
	Array_siz_delete(region);
	
	// Delete flagging regions
	Array_siz_delete(flag_regions);
	
	// Delete input parameters
	Parameter_delete(par);
	
	// Delete input file paths
	Path_delete(path_data_in);
	Path_delete(path_gain_in);
	Path_delete(path_noise_in);
	Path_delete(path_weights_in);
	Path_delete(path_mask_in);
	Path_delete(path_flag_cat);
	
	// Delete output file paths
	Path_delete(path_cat_ascii);
	Path_delete(path_cat_xml);
	Path_delete(path_cat_sql);
	Path_delete(path_mask_out);
	Path_delete(path_mask_2d);
	Path_delete(path_mask_raw);
	Path_delete(path_noise_out);
	Path_delete(path_filtered);
	Path_delete(path_mom0);
	Path_delete(path_mom1);
	Path_delete(path_mom2);
	Path_delete(path_chan);
	Path_delete(path_rel_plot);
	Path_delete(path_rel_cat_n);
	Path_delete(path_rel_cat_p);
	Path_delete(path_skel_plot);
	Path_delete(path_flag);
	Path_delete(path_cubelets);
	
	// Delete source catalogue
	Catalog_delete(catalog);
	
	// Print status message
	status("Pipeline finished.");
	
	return ERR_SUCCESS;
}
