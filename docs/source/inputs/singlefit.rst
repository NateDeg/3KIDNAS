Single Galaxy Fitter Input Files
=================================

The SingleGalaxyFitter program is perhaps the most complex of the 3KIDNAS components and the various input files reflect this complexity.  However, when running 3KIDNAS in pipeline mode, most of these options do not need to be adjusted.

As with the BootstrapSampler program, the format of the SingleGalaxyFitter files is fixed.  It ultimately takes 2 distinct input files.

Main Input
--------------

An example of the main input file for the SingleGalaxyFitter code is /Inputs/SingleFitInput.txt.  The key inputs in the main file are

1. Line 2 gives the name on the data cubelet.
2. Line 4 gives the name of the second input file.
3. Line 6 gives the switch for using a mask or whether the code will attempt to construct a threshold mask.  Note that this mask is only used for initial estimates.
4. Line 8 gives the name of a mask when using one.
5. Line 9/10 gives the switch for reading in catalogue values
6. Line 11/12 gives the name of the catalogue file when reading it in.
7. Line 13/14 gives the random seed for the code.
8. Line 15/16 gives the name of the output folder.
9. Line 17/18 gives the name of the galaxy itself.
10. Line 19/20 is a legacy line that gives the noise.  This is overwritten in the code itself.

