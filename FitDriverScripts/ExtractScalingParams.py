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
    if GalaxyDict['ExtendedSDProfile']['ProfileAcceptFlag']==False:
        RHDict=NoWorkableSD()
        return RHDict
    
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
        SDCalcMethod=-1
        #DictUse=GalaxyDict['ExtendedSDProfile']
        RHDict=NoWorkableSD()
        return RHDict
    
    R=DictUse['R_SD']
    SD=DictUse['SURFDENS_FACEON']
    SDUpper=DictUse['SURFDENS_FACEON']+DictUse['SURFDENS_FACEON_ERR']
    SDLower=DictUse['SURFDENS_FACEON']-DictUse['SURFDENS_FACEON_ERR']

    #print("R check", R)
    #print("SD CHeck", SD)
    #print("errs",DictUse['SURFDENS_FACEON_ERR'])
    print("SD ConsistencyCheck", SD)
    print(SDUpper)
    print(SDLower)

   
    

    if len(R) <=1:
        RHI=R[0]
        RHI_Found=False
        R_indx=0
        RHI_Upper=np.nan
        RHI_Lower=np.nan
        RHIFlag=-1


    else:
        RHI,RHI_Found,R_indx=GetSD_Intecept(R,SD,SDLim)
    
        if np.isnan(DictUse['SURFDENS_FACEON_ERR']).any():
            RHI_Lower=np.nan
            RHI_Upper=np.nan
            RHIFlag=-1
        else:
            RHI_Upper,RUFound,R_UpperIndx=GetSD_Intecept(R,SDUpper,SDLim)
            RHI_Lower,RLFound,R_LowerIndx=GetSD_Intecept(R,SDLower,SDLim)
            RHIFlag=0
            if RUFound==False:
                RHI_Upper=np.nan
                RHIFlag=1
            if RLFound==False:
                RHI_Lower=np.nan
                RHIFlag=2
            if RUFound==False and RLFound==False:
                RHIFlag=-1
                RHI_Lower=np.nan
                RHI_Upper=np.nan
        
        if SDCalcMethod==1 and RHI_Found==False:
            RHIFlag=-1
            RHI_Lower=np.nan
            RHI_Upper=np.nan
            RHI=np.nan
     
    RHIArr=np.array([RHI_Lower,RHI,RHI_Upper])
    RCorrArr=np.array([RHI_Lower,RHI,RHI_Upper])
    print("RHI ARRs", RHIArr)
    print(RCorrArr)
    print("RHI Flag",RHIFlag)
    
    
    
    if SDCalcMethod==1:
        Beam=GalaxyDict['ExtendedSDProfile']['BMAJ']
        print(RCorrArr)
        for RCorr in RCorrArr:
            BU=Beam
            i=1
            while BU > RCorr:
                BU=Beam/float(i)
                i+=1
            RCorr=np.sqrt(RCorr**2.-BU**2.)
        print(RCorrArr)

    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
    return RHDict

        
    
def GetSD_Intecept(R,SD_FO,SDLim):
    print("SD Intercept SD", SD_FO)
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
        if indx==-1:
            RHI_Found=False
            RHI=R[-1]
        else:
            RHI_Found=True
    return RHI,RHI_Found,indx
    
    
    
def FindProfileIntersection(X,Y,Lim):
    #   Loop through all radii and check for the point where we go below 1
    #       Go from out to in
    for i in range(len(X)-1,0,-1):
        #print(i,X[i],Y[i],Lim)
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
            return X_Int,indx
    indx=-1
    X_Int=X[0]
    return X_Int,indx
    

    
def GetProfilePoint(X,Y,XTarg):
    print("len profile", X,XTarg)
    #if np.isnan(XTarg):
    #    ProfLimitsFlag=False
    #    YTarg=np.nan
    #    i=len(X)-1
    #    return YTarg,i,ProfLimitsFlag

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

    if ScalingDict['SDMethod']==-1:
        ScalingDict=BadVHResults(ScalingDict)
        return ScalingDict
    if ScalingDict['RHIFlag']==-1:
        ScalingDict=BadVHResults(ScalingDict)
        return ScalingDict
    Model=GalaxyDict['BestFitModel']
    
    R=Model['R']
    VProf=Model['VROT']
    VProfErr=Model['VROT_ERR']

    RHI=ScalingDict['RHIArr'][1]
    
    if len(VProf) <=1:
        VHI=VProf[0]
        RIndx=0
        VHIFlag=-1
    else:
        VHI,RIndx,VHIFlag_Ini=GetProfilePoint(R,VProf,RHI)
    
    RHILimsFlags=ScalingDict['RHIFlag']
    if ScalingDict['RHIFlag']==0:
        PtArray=np.array(ScalingDict['RHIArr'])
    elif ScalingDict['RHIFlag']==1:
        PtArray=np.array([ScalingDict['RHIArr'][0],ScalingDict['RHIArr'][1]])
    elif ScalingDict['RHIFlag']==2:
        PtArray=np.array([ScalingDict['RHIArr'][1],ScalingDict['RHIArr'][2]])
        
    nPts=len(PtArray)
    Mid=np.zeros(nPts)

    Diffs=np.zeros(nPts)
    RHI=ScalingDict['RHIArr'][1]
    Low,RIndx,VFlag=GetProfilePoint(R,VProf-VProfErr,RHI)
    High,RIndx,VFlag=GetProfilePoint(R,VProf+VProfErr,RHI)
    
    VErr_Local=np.abs((High-Low)/2.)

    for i in range(nPts):

        RUse=PtArray[i]
        print("pt",i,RUse,PtArray,ScalingDict['RHIFlag'])
        Mid[i],RIndx,VFlag=GetProfilePoint(R,VProf,RUse)
        
        Diffs[i]=np.abs(Mid[i]-VHI)
       
   
    Avg=np.mean(Diffs)
    
    TotErr=np.sqrt(VErr_Local**2.+Avg**2.)
    
    if VHIFlag_Ini==True and RHILimsFlags==0:
        ScalingDict['VHIFlag']=0
    elif VHIFlag_Ini==False and RHILimsFlags==0:
        ScalingDict['VHIFlag']=-1
        dR=R[1]-R[0]
        #print(RHI,R[-1],dR,(RHI-R[-1])/dR)
        if (RHI-R[-1])/dR <= 0.5:
            ScalingDict['VHIFlag']=-1
        
    elif VHIFlag_Ini==True and RHILimsFlags==1:
        ScalingDict['VHIFlag']=1
    elif VHIFlag_Ini==True and RHILimsFlags==2:
        ScalingDict['VHIFlag']=1
    elif VHIFlag_Ini==False:
        ScalingDict['VHIFlag']=-1

      
    ScalingDict['VHIArr']=np.array([VHI,TotErr])

    
    return ScalingDict

def DistEst(H0,Vel):
    Distance=Vel/H0
    return Distance


def NoWorkableSD():
    RCorrArr=np.array([np.nan,np.nan,np.nan])
    RHIArr=np.array([np.nan,np.nan,np.nan])
    
    SDCalcMethod=-1
    RHIFlag=-1
    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
    return RHDict

def BadVHResults(ScalingDict):

    ScalingDict['VHIArr']=np.array([np.nan,np.nan,np.nan])
    ScalingDict['VHIFlag']=-1
    
    return ScalingDict
