#!/usr/bin/env python3
import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

from . import MCGGenerator as MCG
from . import CubeComparison as CC
from . import CubeInformation as CI

import os

def GetBBaroloFit(CatID,WallCat,Path,FitType,PathType,ObsDict):
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    if PathType == 0 or PathType==2 or PathType ==3:
        OutputDir=Path+"/"+NameSuffix+"/"
    elif PathType ==1 :
        OutputDir=Path+"/output/"+NameSuffix+"/out3/"
        if os.path.isdir(OutputDir) == False:
            OutputDir=Path+"/output/"+NameSuffix+"/out2/"
            if os.path.isdir(OutputDir) == False:
                OutputDir=Path+"/output/"+NameSuffix+"/out1/"

    
    ParamFile=[OutputDir+"ringlog2.txt",OutputDir+"ringlog1.txt"]
    DensFile=OutputDir+"densprof.txt"
    print(ParamFile[0])

    if PathType==2:
        NameSuffix="NONE"
    if PathType==3:
        NameSuffix="TestCube"

    if FitType == 0:
        CubeFile=OutputDir+NameSuffix+"mod_azim.fits"
    elif FitType == 1:
        CubeFile=OutputDir+NameSuffix+"mod_local.fits"



    MainParams,FITAchieved=LoadBBaroloParamFile(ParamFile)
    try:
        RR,SD=LoadBBaroloSurfDens(DensFile)
    except:
        FITAchieved='False'

    if FITAchieved=='True':
        BBaroloDict=BBaroloParamsToDict(MainParams,SD,ObsDict)
        BBaroloDict['R_SD']=RR

        MCG_BBaroloModelCube=MCG.MakeMCGModel(BBaroloDict,ObsDict['CubeHeader'])
        ResidTot,chi2=CC.CubeCompare(ObsDict,MCG_BBaroloModelCube)
        BBaroloDict['RESID']=ResidTot
        BBaroloDict['CHI2']=chi2

        ModelCube=CI.MakeModelPV(CubeFile,ObsDict,WallCat,CatID)
        BBaroloDict['CUBE']=ModelCube

    else:
        BBaroloDict=BadBBaroloFit()



    return BBaroloDict

def LoadBBaroloParamFile(FileName):
    FITAchieved='True'
    try:
        Params=np.loadtxt(FileName[0],skiprows=1)
    except:
        try:
            Params=np.loadtxt(FileName[1],skiprows=1)
        except:
            FITAchieved='False'
            Params=[]
    return Params,FITAchieved

def LoadBBaroloSurfDens(FileName):
    RR,SD=np.loadtxt(FileName,skiprows=15,usecols=(0,9),unpack='True')
    return RR,SD

def BBaroloParamsToDict(Params,SD,ObsDict):
    TRParamsDict={'R':Params[:,1],'XCENTER':Params[:,9] \
            ,'YCENTER':Params[:,10], 'INCLINATION':Params[:,4] \
            ,'POSITIONANGLE':Params[:,5], 'VSYS':Params[:,11], 'VROT':Params[:,2]\
            ,'VRAD':Params[:,12], 'VDISPERSION':Params[:,3],'Z0':Params[:,7]\
            ,'SURFDENS':SD[:],'FITAchieved':'True'}
  
    RA,DEC=BBaroloCentConversion(TRParamsDict['XCENTER'],TRParamsDict['YCENTER'],ObsDict)
    TRParamsDict['RA']=RA
    TRParamsDict['DEC']=DEC
    
    return TRParamsDict

def BBaroloCentConversion(X,Y,ObsDict):
    w=ObsDict['CubeWCS']
    RA=[]
    DEC=[]
    for i in range(np.shape(X)[0]):
        pixcrd=[[X[i], Y[i],0]]
        CentCoord=w.wcs_pix2world(pixcrd,0)
        RA.append(CentCoord[0,0])
        DEC.append(CentCoord[0,1])

    return RA,DEC



def BadBBaroloFit():
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':'False','RESID':[],'CHI2':[]}
    return Disk



