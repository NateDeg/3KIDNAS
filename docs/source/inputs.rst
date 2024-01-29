Inputs

Input files for 3KIDNAS
-------------


Pipeline Input
-------------

An example of the input file needed for the full pipeline is 
/Inputs/Sample_KIDNASPipelineInput.py

The input file is a python file and requires that 7 variables be set.  They are

1) CatName == The name of the SoFiA-2 catalogue containing the source names and parameters.  This file should be a fits file.

2) SourceFolder == The name of the folder containing all data cubelets and mask files.

3) TargFolder == The name of the folder that the pipeline will place all outputs into.

4) KinTR == This is a keyword used to sort out different runs of data.  It is important when considering how to share the results of a pipeline run to a data archive.

5) nBootstraps == The number of bootstrap resamples generated and fit for each galaxy.  This number must be large enough to get an accurate measurement of the model uncertainties.

6) nProcessors == The number of processors used to fit individual galaxies at the same time.  This variable is currently ill-named as it is the really the number of individual galaxies that will be fit at one time.

7) nProcessors_Bootstrap == The number of processors used to fit independent bootstrap resamples for a particular galaxy.

The total number of processors used in a run of the full 3KIDNAS pipeline is nProcessors * nProcessors_Bootstrap



