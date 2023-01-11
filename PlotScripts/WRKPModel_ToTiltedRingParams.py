import numpy as np
from decimal import Decimal
import argparse

import LoadModel as LM
import CubeAnalysis as CA
import DiagnosticPlotFncs as DPF

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator


def FileNames():
    ParentFolder="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT2/WRKP/WALLABY_HydraTests_New/"
    
    nModels=1
    ModelFolders=[None]*nModels
    Labels=[None]*nModels
    FinalFiles=[None]*nModels
    
    for i in range(nModels):
        if i ==0:
            ModelFolders[i]="/Users/nate/Dropbox/WALLABY/DataReleases/Wallaby_NormaDR1V1_HydraDR2V2_KinematicModels/Wallaby_Hydra_DR2_KinematicModels_v2/WALLABY_J100342-270137/"
            FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_AvgModel_v2.txt"
            #ModelFolders[i]="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/WALLABY_J100342-270137/"
            Labels[i]="Model"
            #FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_AvgModel_v2.txt"

        
    FolderDict=locals()
    return FolderDict


def WriteTiltedRingFile(ParamDict):
    print("Writing TR file")
    MCGFileName="TestTRFile.txt"
    
    nRings=np.shape(ParamDict['R'])[0]
    Rwidth=ParamDict['R'][1]-ParamDict['R'][0]
    file=open(MCGFileName,"w")
    
    
    file.write("#   Number of Rings in the model\n")
    file.write(str(nRings)+"\n")
    
    file.write("#    The cloud mode you are using\n")
    file.write("0 \n")

    file.write("#    The base cloud surface density\n")
    file.write("100. \n")
    
    file.write("#    Central Position Switch (0=degrees, 1=arcsec,2=pixels), Inc/PA Unit switch (0=degrees, 1=arcsec), Velocity Units (0=m/s, 1=km/s), Brightness Units (0=Jy km/s arcsec^-2)\n")
    file.write("2 \t 0 \t 1 \t  0 \n")
    
    file.write("#    The parameters in each radial bin\n")
    file.write("#    Rmid    Rwidth    Xcent    Ycent    Inc    PA    VSys    VRot    VRad    Vvert    VDisp    dvdz    Sigma        z0    zGradStart\n")
    
    ParamDict['VRAD']=[0]*nRings
    ParamDict['VDISPERSION']=[0]*nRings
    ParamDict['Z0']=[0]*nRings
    
    ConvFactor=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    print("SD Conv Factor", ConvFactor)
    
    for i in range(nRings):
        ParamDict['VDISPERSION'][i]=10.
        file.write(str(ParamDict['R'][i])+"\t"+str(Rwidth) +"\t" + str(ParamDict['XCENTER'][0]) \
                   + "\t"+ str(ParamDict['YCENTER'][0]) +"\t"\
                   + str(ParamDict['INCLINATION'][i])\
                   + "\t" +str(ParamDict['POSITIONANGLE'][i]) +"\t"\
                   + str(ParamDict['VSYS'][i])+ "\t" + str(ParamDict['VROT'][i]) \
                   +"\t" + str(ParamDict['VRAD'][i]) + "\t" + "0." +"\t" \
                   + str(ParamDict['VDISPERSION'][i]) + "\t"+ "0." + "\t" \
                   + str(ParamDict['SURFDENS'][i])+"\t" + str(ParamDict['Z0'][i])\
                   + "\t" + str(5.*ParamDict['Z0'][i]) +   "\n")
    
    
    file.close()
    

def Main():
    #   Get the set of arguments needed to produce the model cube
    
    FileDict=FileNames()
    #print(FileDict)
    

    FinalModel=[None]*FileDict['nModels']
    
    for i in range(FileDict['nModels']):
        FinalModel[i]=LM.LoadBestFitModelFile(FileDict['FinalFiles'][i])
    
    
    print(FinalModel[0]['SURFDENS'])
    FinalModel[0]['SD_Obs']=FinalModel[0]['SURFDENS']
    FinalModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    print(FinalModel[0]['SURFDENS'])
    
    WriteTiltedRingFile(FinalModel[0])
    

    
    
    
    
Main()

