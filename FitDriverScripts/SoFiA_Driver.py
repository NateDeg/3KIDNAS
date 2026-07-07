import numpy as np
from decimal import Decimal
import argparse
import os as os
import multiprocessing as mp
import pandas as pd
import copy as copy

import astropy
from astropy.io import fits

from . import GeometryEstimates as GE
 
def LoadSoFiATemplate(DefaultTemplate):
    DefFile=open(DefaultTemplate,'r')
    SoFiATemplate=DefFile.readlines()
    
    DefFile.close()
    
    #for i in range(len(SoFiATemplate)):
    #    print(i,SoFiATemplate[i])
    return SoFiATemplate
    
def WriteTempSoFiAFile(SoFiATemplate,GalaxyDict):
    #print("writing SoFiA file",ObjDict['SourceCubeName'],os.path.isfile(ObjDict['SourceCubeName']))
    CubeUse=GalaxyDict['CubeNameU']
    SoFiATemplate[17]="input.data                 = "+CubeUse+"\n"
    
    SoFiATemplate[142]="output.directory           = "+ GalaxyDict['SoFiAFolder']+"\n"
    SoFiATemplate[143]="output.filename            = "+ GalaxyDict['ObjNameU']+"\n"
    SoFiATemplate[149]="output.writeMask           =  true \n"
    
    TempFileName="SofiAParams_"+str(GalaxyDict['ObjNameU'])+".txt"
    TempFile=open(TempFileName,'w')
    for i in range(len(SoFiATemplate)):
        TempFile.write(SoFiATemplate[i])
    TempFile.close()
    return TempFileName

def RunSoFiA(GeneralDict,GalaxyDict):
    #   Copy the SoFiA template for modification
    CurrTemplate=copy.deepcopy(GeneralDict['SoFiATemplate'])
    #   Make the SoFiA folder needed
    os.makedirs(GalaxyDict['SoFiAFolder'], exist_ok=True)
    #   Now write the SoFiA parameter file
    SoFiARunFile=WriteTempSoFiAFile(CurrTemplate,GalaxyDict)
    RunCMD=GeneralDict['SoFiAExecPath']+" "+SoFiARunFile
    os.system(RunCMD)
    ClnCmd="rm "+SoFiARunFile
    os.system(ClnCmd)
    
    #   After cleaning everything, store the name of the final catalogue file and mask file
    GalaxyDict['SoFiA_CatFileName']=GalaxyDict['SoFiAFolder']+GalaxyDict['ObjNameU']+"_cat.txt"
    GalaxyDict['MaskNameU']=GalaxyDict['SoFiAFolder']+GalaxyDict['ObjNameU']+"_mask.fits"
    return GalaxyDict
    
def LoadSoFiAOutput(ObjDict):

    #   For now, we need to set the method for getting the inclination here.  In the future this should be a runtime option.
    IncMethod='WALLABY_Like'
    #   Start by setting the file name of the SoFiA catalogue
    File=ObjDict['SoFiA_CatFileName']
    #   Check that the SoFiA file exists
    if os.path.exists(File):
        #   If so, load in the relevant inputs
        FSum,Ell_Maj,Ell_Min,PA,X,Y,Z=np.loadtxt(File,unpack=True,skiprows=13,usecols=(15,27,28,33,3,4,5))
        #   Note that the SoFiA run was successful
        SoFiASucces=True
    else:
        #   If the file doesn't exist, send back an empty geometry dictionary.
        Geo={'Name':None ,'Ell_Maj':None,'Ell_Min':None,'PA':None,'IncApprox':None,'X':None,'Y':None,'Z':None,'SoFiASucces':False,'MaskVal':0}
        return Geo
    #   It is possible that the SoFiA run will find multiple entries rather than a single run.
    #       As such, we'll start by checking if Ell_Maj is a float or not
    if isinstance(Ell_Maj,float)  :
        #   If it is a float, we don't have to sort anything and can set the useable params
        WorkingDict={'Ell_MajU':Ell_Maj,'Ell_MinU':Ell_Min,'XU':X,'YU':Y,'ZU':Z,'MaxFluxIndx':0,'PAU':PA}
    else:
        #   If it is an array, we need to get the entry with the maximum flux
        MaxFluxIndx=np.where(FSum==np.max(FSum))[0][0]
        #   Now we'll set the working dictionary based on this entry
        WorkingDict={'Ell_MajU':Ell_Maj[MaxFluxIndx],'Ell_MinU':Ell_Min[MaxFluxIndx],'XU':X[MaxFluxIndx],'YU':Y[MaxFluxIndx],'ZU':Z[MaxFluxIndx],'MaxFluxIndx':MaxFluxIndx,'PAU':PA[MaxFluxIndx]}
    #   Now that we have these set up, we do need the beam size in pixels
    BeamPix=np.abs(ObjDict['CubeHeader']['BMAJ']/ObjDict['CubeHeader']['CDELT1'])
    #   And set up a temporary dictionary needed for the GeometryEstimates function
    TempSoFiADict={'ell_maj':WorkingDict['Ell_MajU'],'ell_min':WorkingDict['Ell_MinU'],'kin_pa':WorkingDict['PAU'],'ell_pa':WorkingDict['PAU']}
    #   Get the position angle and inclination estimate
    PA,IncApprox=GE.GetGeometryEstimates(TempSoFiADict,BeamPix,IncMethod)
    #   Store everything needed in a dictionary
    Geo={'Ell_Maj':WorkingDict['Ell_MajU'],'Ell_Min':WorkingDict['Ell_MajU'],'PA':PA,'IncApprox':IncApprox,'X':WorkingDict['XU'],'Y':WorkingDict['YU'],'Z':WorkingDict['ZU'],'SoFiASucces':SoFiASucces,'MaskVal':WorkingDict['MaxFluxIndx']+1}
    #   And return the dictionary
    return Geo


def AdjustMaskFile(MaskFile,TargVal):
    print(MaskFile,os.path.isfile(MaskFile),TargVal)
    Mask=fits.open(MaskFile)
    
    MData=Mask[0].data
    #print(np.nansum(MData))
    MDataNew=(MData == TargVal).astype(int)
    #print(np.nansum(MDataNew))
    Mask[0].data=MDataNew
    Mask.writeto(MaskFile,overwrite=True)
    Mask.close()

def WriteSoFiACatFileForWRKP(GalaxyDict):
    print("Writing WRKP SoFiA catalogue file")
    #   Start by creating a string variable that holds all the estimated values and write in the estimated SoFiA inclination
    CatString="#\t Inclination \n"
    CatString+=str(GalaxyDict['Inc_EstimateU'])+"\n"
    #   Next write in the position angle esimtations
    CatString+="#\t Position Angle \n"
    CatString+=str(GalaxyDict['PA_EstimateU'])+"\n"
    
    #   Now we need to name the SoFiA Catalogue file
    CatFile=GalaxyDict['ObjNameU']+"_SoFiAShapeEstimate.txt"
    file=open(CatFile,'w')
    file.write(CatString)
    file.close()
    return CatFile
