import sys as sys
import os as os
import copy as copy
import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

from . import MakeBootstrapSample as MBS
from . import SoFiA_Driver as SD
from . import RunWRKP as RW
from . import ReadWRKPFit as RWF

def GetCubeHeader(CFile):
    CubeHDU=fits.open(CFile)
    CubeHeader=CubeHDU[0].header
    CubeHDU.close()
    return CubeHeader
    
def SetErrorAnalysisVariables(GalaxyDict):
    print(GalaxyDict['ObjName'])
    
    #   Start by setting the name of the WRKP results folder
    GalaxyDict['WRKP_ResultsFolder']=GalaxyDict['TargFolder']+GalaxyDict['ObjName']+"/"
   
    #   Name a folder for the bootstrap files
    GalaxyDict['BootstrapFolder']=GalaxyDict['WRKP_ResultsFolder']+"BootstrapCubes/"
    #   Also name a folder for the SoFiA results
    GalaxyDict['SoFiAFolder']=GalaxyDict['WRKP_ResultsFolder']+"SoFiARuns/"
    
    GalaxyDict['TargFolderU']=GalaxyDict['BootstrapFolder']+"Fits/"
    os.makedirs(GalaxyDict['TargFolderU'],exist_ok=True)
    
    #   The cube header is necessary for calculating the velocity channel of VSys
    GalaxyDict['CubeHeader']=GetCubeHeader(GalaxyDict['CubeName'])
    
    return GalaxyDict
    
def SetBootstrapVariableValues(GalaxyDict,Step):
    GalaxyDict['ObjNameU']=GalaxyDict['ObjName']+"_Bootstrap_"+str(Step)
    
    return GalaxyDict
    
 
def GetBootstrapModel(GeneralDict,GalaxyDict,Step):
    #   Start by naming the bootstrap files
    GalaxyDict=MBS.SetBootstrapVariableValues(GalaxyDict,Step)
    #   Now Make a bootstrap sample
    GalaxyDict=MBS.MakeBootstrapSample(GeneralDict,GalaxyDict,Step)
    #   Now that we have the bootstrap cube, run SoFiA on it
    GalaxyDict=SD.RunSoFiA(GeneralDict,GalaxyDict)
    #   Now we need to load in the SoFiA Catalogue
    SoFiAGeo=SD.LoadSoFiAOutput(GalaxyDict)
    #   Store the SoFiA geometry estimates into the Galaxy Dictionary
    if SoFiAGeo['SoFiASucces']:
        GalaxyDict['Inc_EstimateU']=SoFiAGeo['IncApprox']
        GalaxyDict['PS_EstimateU']=SoFiAGeo['PA']
        #   Also adjust the mask file as is needed
        SD.AdjustMaskFile(GalaxyDict['MaskNameU'],SoFiAGeo['MaskVal'])
        #   Now we can go through the process of running WRKP again
        GalaxyDict['SoFiAShapeFile']=SD.WriteSoFiACatFileForWRKP(GalaxyDict)
        #       The first step is create a 'catalogue' file with the initial estimates for the inclination and position angle
        GalaxyDict['SoFiAShapeFile']=SD.WriteSoFiACatFileForWRKP(GalaxyDict)
        #   Now run WRKP
        GalaxyDict=RW.RunWRKP(GeneralDict,GalaxyDict,1)
        #   With WRKP now completed, the fit must be ingested
        CurrModel=RWF.ReadWRKPOutputFile(GeneralDict,GalaxyDict)
    else:
        CurrModel={'FITAchieved':False}

    return CurrModel
    
    
def EstimateUncertaintiesFromBootstraps(GeneralDict,GalaxyDict,BootstrapModels):
    print("Estimating Uncertainties")
    
    print(GalaxyDict['BestFitModel'].keys())
    #   Start by listing the various model parameters that must be averaged together
    AvgKeys=['XCENTER','YCENTER','INCLINATION','POSITIONANGLE','VSYS','VDISP','VROT','SURFDENS','RA','DEC','SURFDENS_FACEON']
    
    #   Loop through each key
    #for i in range(1):
    #    key=AvgKeys[i]
    for key in AvgKeys:
        #   Set the error key
        ErrKey=key+"_ERR"
        GalaxyDict['BestFitModel']=BootsrapErrors_ForParam(key,ErrKey,GalaxyDict['BestFitModel'],BootstrapModels)
        #print(key,GalaxyDict['BestFitModel'][ErrKey])
    return GalaxyDict['BestFitModel']

def BootsrapErrors_ForParam(key,ErrKey,Model,BootstrapModels):
    
    #   First get the number of radial bins
    if key=='SURFDENS' or key == 'SURFDENS_FACEON':
        nR=len(Model['R_SD'])
    else:
        nR=len(Model['R'])
    #   Next set the number of bootstraps for the local calculations
    nBootstraps=len(BootstrapModels)
    #   Set the number of successful bootstraps needed for a ring
    TargSuccessFrac=0.
    nSuccessNeeded=TargSuccessFrac*nBootstraps
    #   Set a keyword for the mean value
    MeanKey=key+"_BS_MEAN"

    #   Initialize the mean and error arrays
    MeanArr=np.zeros(nR)
    ErrArr=np.zeros(nR)
    CurrIndx=-1
    #for i in range(1):
    for i in range(nR):
        CurrIndx+=1
        BootstrapArr=[]
        for j in range(len(BootstrapModels)):
            CurrModel=BootstrapModels[j]
            try:
                ModelVal=BootstrapModels[j][key][i]
                BootstrapArr.append(ModelVal)
            except:
                continue
        BootstrapArr=np.array(BootstrapArr)
        nFits=len(BootstrapArr)
        
        if key == 'POSITIONANGLE':
            for j in range(1,nFits):
                #print(i,j,BootstrapArr[j],BootstrapArr[0])
                if BootstrapArr[j]-BootstrapArr[0] > 180.:
                    BootstrapArr[j]=BootstrapArr[j]-360.
                elif BootstrapArr[j]-BootstrapArr[0] < -180.:
                    BootstrapArr[j]=BootstrapArr[j]+360.
                #print("After corr",i,j,BootstrapArr[j],BootstrapArr[0])
        MeanArr[CurrIndx]=np.mean(BootstrapArr)
        ErrArr[CurrIndx]=np.std(BootstrapArr)
        print(key,i,nFits,len(BootstrapModels),nSuccessNeeded)
        if nFits <= nSuccessNeeded:
            MeanArr=np.delete(MeanArr,CurrIndx)
            ErrArr=np.delete(ErrArr,CurrIndx)
            Model[key]=np.delete(Model[key],CurrIndx)
            CurrIndx-=1
            

    
    Model[ErrKey]=ErrArr
    Model[MeanKey]=MeanArr
    Model[ErrKey]=np.sqrt(ErrArr**2.+(MeanArr-Model[key])**2.)
    
    #   Round the averages before sending them back
    RoundMeasures(key,Model)
    Model['nFits']=nFits
    return Model


def RoundMeasures(key,Model):
    ErrKey=key+"_ERR"
    MeanKey=key+"_BS_MEAN"
    nRound=1
    if key=='XCENTER' or key=='YCENTER':
        nRound=4
    elif key=='RA' or key=='DEC':
        nRound=6
    elif key=='SURFDENS' or key=='SURFDENS_FACEON':
        nRound=2
    nRound_Err=nRound
    for i in range(len(Model[key])):
        Model[key][i]=round(Model[key][i],nRound)
        Model[ErrKey][i]=round(Model[ErrKey][i],nRound)
        Model[MeanKey][i]=round(Model[MeanKey][i],nRound)
    return Model
