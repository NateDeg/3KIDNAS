#
#    This module contains routines to calculate an extend Surface density profile from the moment 0 map
#

import numpy as np
from decimal import Decimal
import argparse
import os as os
import multiprocessing as mp
import pandas as pd
import copy as copy

import astropy
from astropy.io import fits


def CalcScalingParams(GalaxyDict):

    #   Get RHI and limits
    ScalingDict=CalcRHI(GalaxyDict)
    #   Get VHI and limits
    ScalingDict=CalcVHI(GalaxyDict,ScalingDict)
    #   Convert RHI to kpc
    H0=70.
    Dist=DistEst(H0,GalaxyDict['BestFitModel']['VSYS'][0])
    print("Estimated distance", GalaxyDict['BestFitModel']['VSYS'][0],Dist)
    
    RHI_kpc=ScalingDict['RHI_CorrArr']/206265*Dist*1000.
    ScalingDict['RHI_kpc']=RHI_kpc
    
    GalaxyDict['ScalingDict']=ScalingDict
    return GalaxyDict
    


def CalcRHI(GalaxyDict):
    print("Calc RHI")
    
    SDLim=1.
    #   First check on the 3D modelled SD profile
    Model=GalaxyDict['BestFitModel']
    SD3D=Model['SURFDENS_FACEON']
    SD2D=GalaxyDict['ExtendedSDProfile']['SURFDENS_FACEON']
    SDMax=np.nanmax(SD3D)
    SDMin=np.nanmin(SD3D)
    print("3d Min Max check", SDMin,SDMax)
    
    if SDMax >= SDLim and SDMin <= SDLim:
        SDCalcMethod=0
        DictUse=Model
    #       If the 3D method doesn't bracket the SD values, check the extended profile
    else:
        SDCalcMethod=1
        DictUse=GalaxyDict['ExtendedSDProfile']
    R=DictUse['R_SD']
    SD=DictUse['SURFDENS_FACEON']
    SDUpper=DictUse['SURFDENS_FACEON']+DictUse['SURFDENS_FACEON_ERR']
    SDLower=DictUse['SURFDENS_FACEON']-DictUse['SURFDENS_FACEON_ERR']

    print("R check", R)
    print("SD CHeck", SD)

    if len(R) <=1:
        RHI=R[0]
        RHI_Found=False
        R_indx=0
        SDCalcMethod=5
        RHI_Upper=np.nan
        RHI_Lower=np.nan
    else:
        RHI,RHI_Found,R_indx=GetSD_Intecept(R,SD,SDLim)
    
        if np.isnan(DictUse['SURFDENS_FACEON_ERR']).any():
            RHI_Lower=np.nan
            RHI_Upper=np.nan
            SDCalcMethod=4
        else:
            RHI_Upper,RUFound,R_UpperIndx=GetSD_Intecept(R,SDUpper,SDLim)
            RHI_Lower,RLFound,R_LowerIndx=GetSD_Intecept(R,SDLower,SDLim)
        
        if SDCalcMethod==1 and RHI_Found==False:
            SDCalcMethod=2

     
    RHIArr=np.array([RHI_Lower,RHI,RHI_Upper])
    RCorrArr=np.array([RHI_Lower,RHI,RHI_Upper])
    
    
    
    if SDCalcMethod==1:
        Beam=GalaxyDict['ExtendedSDProfile']['BMAJ']
        for RCorr in RCorrArr:
            BU=Beam
            i=1
            while BU > RCorr:
                BU=Beam/float(i)
                i+=1
            RCorr=np.sqrt(RCorr**2.-BU**2.)

    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod}
    return RHDict

        
    
def GetSD_Intecept(R,SD_FO,SDLim):
    #   Check the FO limits
    if np.nanmin(SD_FO) > SDLim:
        RHI=R[-1]
        RHI_Found=False
        indx=-1
    elif np.nanmax(SD_FO) < SDLim:
        SDMax=np.nanmax(SD_FO)
        indx=np.argmax(SD_FO)
        RHI=R[indx]
        RHI_Found=False
    else:
    #       If the limits are good, get RHI
        RHI,indx=FindProfileIntersection(R,SD_FO,SDLim)
        RHI_Found=True
    return RHI,RHI_Found,indx
    
    
    
def FindProfileIntersection(X,Y,Lim):
    #   Loop through all radii and check for the point where we go below 1
    #       Go from out to in
    for i in range(len(X)-1,0,-1):
        #   Once the profile goes above one, find the point where the SD goes below the limit
        if Y[i-1] > Lim:
            #   Select the points
            x1=X[i-1]
            x2=X[i]
            y1=Y[i-1]
            y2=Y[i]
            #   Get the slope
            m=(y2-y1)/(x2-x1)
            #   Find RHI
            dY=Lim-y1
            dX=dY/m
            X_Int=x1+dX
            indx=i-1
            break
    return X_Int,indx

    
def GetProfilePoint(X,Y,XTarg):
    ProfLimitsFlag=True
    if XTarg<=X[0]:
        i=0
        ProfLimitsFlag=False
    elif XTarg>X[-1]:
        i=len(X)-2
        ProfLimitsFlag=False
    else:
        for j in range(len(X)-1):
            if XTarg > X[j] and XTarg <= X[j+1]:
                i=j
                break

    x1=X[i]
    x2=X[i+1]
    y1=Y[i]
    y2=Y[i+1]
    #   Get the slope
    m=(y2-y1)/(x2-x1)
    dX=XTarg-x1
    YTarg=y1+m*dX
    return YTarg,i,ProfLimitsFlag


def CalcVHI(GalaxyDict,ScalingDict):
    print("Calc VHI")

    Model=GalaxyDict['BestFitModel']
    
    R=Model['R']
    VProf=Model['VROT']
    VProfErr=Model['VROT_ERR']

    RHI=ScalingDict['RHIArr'][1]
    
    if len(VProf) <=1:
        VHI=VProf[0]
        RIndx=0
        VHIFlag=2
    else:
        VHI,RIndx,VHIFlag=GetProfilePoint(R,VProf,RHI)
        
        
    if ScalingDict['SDMethod']==4 or len(VProf)<=1:
        VMin=np.nan
        VMax=np.nan
    else:
        VUpperArr=np.zeros(3)
        VLowerArr=np.zeros(3)
        for i in range(3):
            RU=ScalingDict['RHIArr'][i]
            VUpperArr[i],RIndxSpec,VHIFlagSpec=GetProfilePoint(R,VProf+VProfErr,RU)
            VLowerArr[i],RIndxSpec,VHIFlagSpec=GetProfilePoint(R,VProf-VProfErr,RU)
        
        VMax=np.nanmax(VUpperArr)
        VMin=np.nanmin(VLowerArr)
    
    
    ScalingDict['VHIArr']=np.array([VMin,VHI,VMax])
    if VHIFlag:
        ScalingDict['VHIFlag']=0
    else:
        ScalingDict['VHIFlag']=1
    
    return ScalingDict

def DistEst(H0,Vel):
    Distance=Vel/H0
    return Distance
