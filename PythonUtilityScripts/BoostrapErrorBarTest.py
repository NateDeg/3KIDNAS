import numpy as np
import multiprocessing as mp

import LoadCatalogue as LC
import CatalogueClass as CC
import CubeAnalysis as CA

import FitComparisons as FC


#       Initialize location of data files, folder and catalogue
DataFolder="WALLABY_Hydra_DR1/"
BaseName="WALLABY_PS_Hya_DR1"
CatalogueFile=DataFolder+BaseName+"_catalog.xls"
WallabyCat=LC.LoadWallabyCat(CatalogueFile)

#Define the pixel and beam size if needed
pixelSize=1.66666666667E-03
beamSize=8.33333333333E-03
PixelBeamRatio=beamSize/pixelSize
BeamArea=(np.pi*(beamSize**2.)/4.)/(pixelSize**2.)
#print("Beam Area in Pixels^2",BeamArea)

RestFreq=1.42040575179E+09

#   Load in the full catalogue
WallabyCat.ObjectSize_Beams=WallabyCat.ell_maj/PixelBeamRatio

ModellingSwitch=np.array([113])
ModellingSwitch=ModellingSwitch-1


FitPaths=[
          ,["/Users/nate/Dropbox/KinematicPipeline_Fortran_Development/BBaroloBatchFits_MCG_J103523-244506/",0,1, "BS Trial",3,100]
          ]



VelPath=["VelocityCubes"]

StrDict={'DataFolder':DataFolder,'BaseFileName':BaseName,'FitPaths':FitPaths,'VelCubePath':VelPath}



#for i in range(17,18,1):
for i in range(len(ModellingSwitch)):
    j=ModellingSwitch[i]
    Name=WallabyCat.name[j]
    FC.BEP.BootstrapErrors(j,StrDict,WallabyCat)
