BootstrapSampler Input File
=================================


The BoostrapSampler code requires the user to set up a single input file.  An example of this input file is /Inputs/BootStrapSampleInput.txt.  Unlike the pipeline and single galaxy driver inputs, this file is a text file and must retain the exact format of the example file.  

The contents of the input file are as follows:
1. Line 2 of the file must contain the name of the data cube
2. Line 4 of the file must contain the name of the model cube.
3. Line 6 of the file must contain the base name to be used for naming the bootstrap cube
4. Line 8 of the file is legacy code and does not impact the results.
5. Line 10 of the file gives the cube centre in pixel/cells as well as the inclination and position angle in radians.


