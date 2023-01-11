import numpy as np
from decimal import Decimal
import argparse
import os as os
import multiprocessing as mp
import pandas as pd
import copy as copy

import astropy
from astropy.io import fits
 
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
    File=ObjDict['SoFiA_CatFileName']
    if os.path.exists(File):
        FSum,Ell_Maj,Ell_Min,PA,X,Y,Z=np.loadtxt(File,unpack=True,skiprows=13,usecols=(15,27,28,33,3,4,5))
    else:
        Geo={'Name':None ,'Ell_Maj':None,'Ell_Min':None,'PA':None,'IncApprox':None,'X':None,'Y':None,'Z':None,'SoFiASucces':False,'MaskVal':0}
        return Geo
    
    
    IncApprox=np.arccos(Ell_Min/Ell_Maj)*180./np.pi
    SoFiASucces=True
#IncApprox=np.array([5.1,3.2])
    #print(isinstance(IncApprox,float))
    if isinstance(IncApprox,float)  :
        Geo={'Ell_Maj':Ell_Maj,'Ell_Min':Ell_Min,'PA':PA,'IncApprox':IncApprox,'X':X,'Y':Y,'Z':Z,'SoFiASucces':SoFiASucces,'MaskVal':1}
    else:
        MaxFluxIndx=np.where(FSum==np.max(FSum))[0][0]
        Geo={'Ell_Maj':Ell_Maj[MaxFluxIndx],'Ell_Min':Ell_Min[MaxFluxIndx],'PA':PA[MaxFluxIndx],'IncApprox':IncApprox[MaxFluxIndx],'X':X[MaxFluxIndx],'Y':Y[MaxFluxIndx],'Z':Z[MaxFluxIndx],'SoFiASucces':SoFiASucces,'MaskVal':MaxFluxIndx+1}

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
