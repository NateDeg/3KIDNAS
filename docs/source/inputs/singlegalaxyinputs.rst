Single Galaxy Input File
=================================



An example of an input file needed to run 3KIDNAS for a single galaxy is 
/Inputs/Sample_KIDNAS_SingleGalaxyInput.py

As with the full pipeline input, there are a number of parameters that must be set in this file.

1) CubeName == The name of the cube that will be fit.

2) MaskName == The name of the mask of the data cube.  These must have the same dimensions.

3) ObjName == The name of the galaxy that will be fit.  This will be used for naming output files and folders as well as encoded into the model cubelets.

4) TargFolder == The name of the parent folder that will hold the outputs from the fit.

5) PA_Estimate == The initial estimate for the position angle in degrees

6) Inc_Estimate == The initial estimate for the inclination in degrees

7) nBootstraps == The number of bootstraps fits to be done to calculate the model uncertainties.  This number must be large enough to produce accurate uncertainties.

8) nProcessors_Bootstraps == The number of processors to be used to fit the bootstraps in parallel.

