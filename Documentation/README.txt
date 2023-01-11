Simple Mock Galaxy Generator

*		PURPOSE
This suite of programs is intended to generate mock galaxies using a basic tilted ring
model.  There are three main use cases for the code:
1) Making a mock data cube of a single galaxy from tilted ring parameters
2) Making a mock data cube of a single galaxy based on the HI mass and scaling relations
3) Generating a suite of mock galaxy observations by sampling from a suite of masses.

*		STRUCTURE
The code is broken into three programs:
1) MockCubeGenerator — a Fortran program designed to make mock data cubes from
			tilted ring model parameters.
2) make_galaxy — a python program that generates tilted ring parameters 
			from scaling relations.
3) make_catalogue - a python program that generates a suite of tilted ring models.

*		DEPENDENCIES
The dependencies are, at the time of writing, fairly common and standard libraries.

MockCubeGenerator requires openmpi, the fftw3 library, and the cfitsio library.

make_galaxy requires python3, numpy, scipy, astropy, and matplotlib.  It also requires a successful install of MockCubeGenerator.

make_catalogue has the same requirements as make_galaxy.

*		INSTALLATION/COMPILATION
1) Enter /src/ and open the file makeflags.
2) Set the lines FitsLibLoc, FFTW_INCL, and FITS_LIBS to the locations of 
		cfitsio and fftw3 respectively.
3) Type make

This should build MockCubeGenerator and install it into the folder labelled Programs.

make_galaxy and make_catalogue are python programs and do not need compilation.


*		BASIC USAGE
***			MockCubeGenerator

The MockCubeGenerator program is controlled by a set of input text files.  Example
input files are located in Inputs/

Each run uses three files; a control file, a datacube header file, and a tilted-ring
parameter file.  The program is run from the command line via:
./MockCubeGenerator ControlFile

The example control file is ‘MockCubeGeneratorInput.in’.  It contains the name
of a folder that will hold all the outputs, a switch for the output volume, the 
name of the datacube header file, the name of the file containing the 
tilted ring parameters, and the noise.

The datacube header contains parameters that specify the datacube.  Similarly, the 
tilted ring parameter file contains the specific galaxy model that will be built.
The example files are ’MockDataCubeHeader.in’ and ‘TiltedRingModel.in’.  A more
detailed description of these files and associated parameters is located in a 
later section of this document.


***			make_galaxy

make_galaxy has a decidedly simpler interface than MockCubeGenerator.  It is 
simply run from the command line by:
./make_galaxy

There are two files that control the behavior of this program.  The most important
one is galaxyconfig_MCG.py.  In this file, the user specifies the HI mass of the 
mock galaxy, the number of beams across the major axis, the inclination and 
position angle of the mock galaxy, and whether the code will produce a set of 
plots and mock data cubes. 

The secondary configuration file is:
galaxy_generation_from_scaling_relations/secondary_config.py
That file contains the datacube parameters, the beam parameters, the specific
noise, the volume of outputs, and the more complicated tilted ring parameters.
For the most part the user will not need to adjust these values

***			make_catalogue
N/A

*			DETAILED DESCRIPTIONS
