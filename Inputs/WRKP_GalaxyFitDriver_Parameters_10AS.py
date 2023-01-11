#           This file contains definitions needed for the WRKP_GalaxyFitDriver
#
#   First give the name of the cube that will be fit
CubeName="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT2/BootstrapDriverTesting/10arc_cubes_Nathan/HIPASSJ1255+04b_Vel.fits"
#   Next give the name of the base cube mask
MaskName="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT2/BootstrapDriverTesting/10arc_cubes_Nathan/sofia_HIPASSJ1255+04b/HIPASSJ1255+04b_cubelets/HIPASSJ1255+04b_2_mask.fits"
#   Also give a name for the object for naming purposes
ObjName="HIPASSJ1255+04b"
#   And give a name for the folder that will contain the fitting results
TargFolder="BootstrapDriver_Trial_v1/"
#   And get the initial estimates for the inclination and position angle in degrees
PA_Estimate= 303.441
Inc_Estimate=60.0
#   Also set the number of bootstraps to make and run
nBootstraps= 50

