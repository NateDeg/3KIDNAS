#           This file contains definitions needed for the WRKP_GalaxyFitDriver
#
#   First give the name of the cube that will be fit
CubeName="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Hydra_DR2_VelocityCubes_Jan2022/WALLABY_PS_Hya_DR2_J102019-285220_cube_Vel.fits"
#   Next give the name of the base cube mask
MaskName="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/WALLABY Hydra DR2/WALLABY_PS_Hya_DR2_source_products/WALLABY_J102019-285220_mask.fits.gz"
#   Also give a name for the object for naming purposes
ObjName="WALLABY_PS_Hya_DR2_J102019-285220"
#   And give a name for the folder that will contain the fitting results
TargFolder="BootstrapDriver_Trial_v1/"
#   And get the initial estimates for the inclination and position angle in degrees
PA_Estimate= 67.34
Inc_Estimate=61.2
#   Also set the number of bootstraps to make and run
nBootstraps= 50
#   And set the number of processors to use for the bootstraps
nProcessors_Bootstraps=5
