Installation

Installation Overview
-------------

3KIDNAS is a combination of two distinct Fortran programs linked my a set of python scripts.  The Fortran programs require some 3rd party software and there may be some compatibility issues that have yet to be solved.

The Fortran programs are:
	SingleGalaxyFitter
	BootStrapSampler


Requirements
-------------

The Fortran code requires:
1) gfortran
2) The cfitsio library
3) The fftw3 library

For the users convenience, versions of the cfitsio and fftw3 library have been included.

The python scripts require:
1) python3.9 (the code is currently locked to this version)
2) numpy
3) scipy
4) matplotlib 
5) astropy
6) pandas 
7) spectral_cube 
8) CosmosCanavas
9) SoFiA-2 
10) multiprocessing

Installation
--------------

To install the Fortran programs, begin by installing the third party packages.  This can be done in your own machine and the makeflags file in the src/ directory can be adjusted.  However, the simplest solution is to build them using the specific packages.  In detail:
1) cd third_party/cfitsio
2) make clean
3) ./configure
4) make 
5) cd third_party/fftw-3.3.8
6) make clean
7) ./configure
8) make
9) cd third_party/SoFiA-2-master_2_5_1
10) make clean
11) ./compile.sh
12) make

Once all third part software is installed locally, the Fortran programs can be installed via
1) cd src/
2) make clean
3) make

Once these are run, the executables should be located in the Programs/ folder.  

Issues
--------------

There is a known issue in the Fortran programs related to the interface between the cfitsio library and whether the os is a Mac M1+ or a Linux environment.  We are working to solve this issue.  While the programs will install, they will not make cubes of the correct size.  If you are using Linux, please ensue that:

1. In src/Inputs/DataCubeInput.f that line 161 is
.. code-block:: bash
	integer nInts(3)
In a Mac, that line should be
.. code-block:: bash
	integer*8 nInts(3)
2. Similarly, in src/Outputs/DataCubeOutputs, line 118 should be
.. code-block:: bash
	integer naxesT(3),naxisT
for Linux and 
.. code-block:: bash
	integer*8 naxesT(3),naxisT
for a Mac.  

