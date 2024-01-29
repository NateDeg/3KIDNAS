Installation

3KIDNAS is a combination of two distinct Fortran programs linked my a set of python scripts.  The Fortran programs require some 3rd party software and there may be some compatibility issues that have yet to be solved.

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
