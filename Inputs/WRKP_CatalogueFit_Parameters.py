#           This file contains definitions needed for the WRKP_CatalogueFitDriver
#
#   First give the name of the Catalogue that will be fit
CatName="/Users/nate/Dropbox/WALLABY/NGC4808_TR1/NGC_4808_DR1_catalog.fits"
#   Next give the name of the source folder
SourceFolder="/Users/nate/Dropbox/WALLABY/NGC4808_TR1/NGC_4808_DR1_products/"
#   And give a name for the folder that will contain the fitting results
#TargFolder="WRKP_Trial/"
TargFolder="NGC4808TR1_WRKPTrial_Sept9_2022/"
#   Also set the number of bootstraps to make and run
nBootstraps= 50
#   Set the number of processors to use for different galaxies
nProcessors=1
#   And the number of processors to use for bootstraps for each individual galaxy
nProcessors_Bootstrap=1

#CatName="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/WALLABY Hydra DR2/WALLABY_PS_Hya_DR2_catalog.xls"
