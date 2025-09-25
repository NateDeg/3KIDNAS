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


def CalcScalingParams(GalaxyDict,BootstrapModels):

    #   Get RHI and limits
    ScalingDict=CalcRHI(GalaxyDict,BootstrapModels)
    

    #   Get VHI and limits
    ScalingDict=CalcVHI(GalaxyDict,ScalingDict,BootstrapModels)

    #   Convert RHI to kpc
    H0=70.
    Dist=DistEst(H0,GalaxyDict['BestFitModel']['VSYS'][0])
    #print("Estimated distance", GalaxyDict['BestFitModel']['VSYS'][0],Dist)
    
    RHI_kpc=ScalingDict['RHI_CorrArr']/206265*Dist*1000.
    ScalingDict['RHI_kpc']=RHI_kpc
    
    GalaxyDict['ScalingDict']=ScalingDict

    return GalaxyDict
    


def CalcRHI(GalaxyDict,BootstrapModels):
    if GalaxyDict['ExtendedSDProfile']['ProfileAcceptFlag']==False:
        RHDict=NoWorkableSD(BootstrapModels)
        return RHDict
    
    SDLim=1.
    SDCalcMethod=0
    count=0
    #   Go through every model and get RHI
    for i in range(-1,len(BootstrapModels)):
        #   Start with the best fitting model
        if i==-1:
            CurrModel=GalaxyDict['BestFitModel']
        else:
            CurrModel=BootstrapModels[i]
        if CurrModel['FITAchieved']==False:
            CurrModel['RHI']=np.nan
            CurrModel['RHIFlag']=-1
            CurrModel['R_Indx']=np.nan
            continue
        #   Set the radius and SD arrays
        RU=CurrModel['R_SD']
        SDU=CurrModel['SURFDENS_FACEON']
        #   Attempt to get RHI fro the model
        CurrModel['RHI'],CurrModel['RHIFlag'],CurrModel['R_indx']=ExtractRHI_NoErr(RU,SDU,SDLim)
        #   Count up all successful bootstrap models
        if i>=0 and CurrModel['RHIFlag']==True:
            count+=1
    #   Check on the best fit model.  If it failed, then no need to get uncertainties, so end here
    if GalaxyDict['BestFitModel']['RHIFlag']==False:
        RHIArr=np.array([np.nan,np.nan,np.nan])
        RCorrArr=RHIArr
        RHIFlag=-1
        RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
        return RHDict
    
    #   Check if there are enough bootstraps with RHI values to get an uncertainty
    RequiredFrac=0.6
    nRequired=int(RequiredFrac*len(BootstrapModels))
    if count < nRequired:
        print("Too few fits with RHI - bad errors", count,nRequired)
        RHIArr=np.array([np.nan,CurrModel['RHI'],np.nan])
        RCorrArr=RHIArr
        RHIFlag=-1
        RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
        return RHDict
    
    #   If we make it to this point there are enough fits to get the uncertainty
    #       Start by initializing an array of RHI values
    RHI_BS_Arr=np.zeros(count)
    j=0
    #   Add all the RHI's to the array
    for i in range(len(BootstrapModels)):
        if BootstrapModels[i]['RHIFlag']==True:
            RHI_BS_Arr[j]=BootstrapModels[i]['RHI']
            j+=1
    #   Calculate the mean and the standard deviation of the bootstrap array
    RHI_BS_Mean=np.mean(RHI_BS_Arr)
    RHI_BS_Err=np.std(RHI_BS_Arr)
    #   Get the difference between the best fit and the mean
    RHI_Diff=GalaxyDict['BestFitModel']['RHI']-RHI_BS_Mean
    #   Get the uncertainty in RHI
    RHI_Err=np.sqrt(RHI_BS_Err**2.+RHI_Diff**2.)
    
    #   Store the appropriate values into an array
    RHIArr=np.array([CurrModel['RHI']-RHI_Err,CurrModel['RHI'],CurrModel['RHI']+RHI_Err])
    RCorrArr=RHIArr
    RHIFlag=0
    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
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
        if indx==-1:
            RHI_Found=False
            RHI=R[-1]
        else:
            RHI_Found=True
    return RHI,RHI_Found,indx
    
    
    
def FindProfileIntersection(X,Y,Lim):
    epsilon=1.e-7
    #   Loop through all radii and check for the point where we go below 1
    #       Go from out to in
    for i in range(len(X)-1,1,-1):
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
            if m < epsilon:
                X_Int=x1
                indx=i-1
            return X_Int,indx
    indx=-1
    X_Int=X[0]
    return X_Int,indx
    

    
def GetProfilePoint(X,Y,XTarg):
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


def CalcVHI(GalaxyDict,ScalingDict,BootstrapModels):

    #   First check if it's a bad model
    if ScalingDict['SDMethod']==-1:
        ScalingDict=BadVHResults(ScalingDict,BootstrapModels)
        return ScalingDict
    if ScalingDict['RHIFlag']==-1:
        ScalingDict=BadVHResults(ScalingDict,BootstrapModels)
        return ScalingDict
    Model=GalaxyDict['BestFitModel']
    
    
    # Now get VHI for the best fit model
    VHI,RIndx,VHIFlag=GetVHIFromProf(Model['R'],Model['VROT'],Model['RHI'])
    Model['VHI']=VHI
    Model['VHIFlag']=VHIFlag
    #   Now loop through the bootstraps and do the same
    count=0
    for i in range(len(BootstrapModels)):
        Model=BootstrapModels[i]
        if Model['RHIFlag']==True:
            VHI,RIndx,VHIFlag=GetVHIFromProf(Model['R'],Model['VROT'],Model['RHI'])
            Model['VHI']=VHI
            Model['VHIFlag']=VHIFlag
            count+=1
        else:
            Model['VHI']=np.nan
            Model['VHIFlag']=False
 
    #   Collect all VHI measures into an array
    VHI_Arr_BS=np.zeros(count)
    j=0
    for i in range(len(BootstrapModels)):
        Model=BootstrapModels[i]
        if Model['VHIFlag']==True:
            VHI_Arr_BS[j]=Model['VHI']
            j+=1
            
    #   Calculate the mean and the standard deviation of the bootstrap array
    VHI_BS_Mean=np.mean(VHI_Arr_BS)
    VHI_BS_Err=np.std(VHI_Arr_BS)
    #   Get the difference between the best fit and the mean
    VHI_Diff=GalaxyDict['BestFitModel']['VHI']-VHI_BS_Mean
    #   Get the uncertainty in RHI
    VHI_Err=np.sqrt(VHI_BS_Err**2.+VHI_Diff**2.)
    ScalingDict['VHIFlag']=0
    ScalingDict['VHIArr']=np.array([GalaxyDict['BestFitModel']['VHI'],VHI_Err])
    return ScalingDict


def DistEst(H0,Vel):
    Distance=Vel/H0
    return Distance


def NoWorkableSD(BootstrapModels):
    RCorrArr=np.array([np.nan,np.nan,np.nan])
    RHIArr=np.array([np.nan,np.nan,np.nan])
    
    SDCalcMethod=-1
    RHIFlag=-1
    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
    for i in range(len(BootstrapModels)):
        BootstrapModels[i]['RHI']=np.nan
        BootstrapModels[i]['RHIFlag']=-1
    
    return RHDict

def BadVHResults(ScalingDict,BootstrapModels):

    ScalingDict['VHIArr']=np.array([np.nan,np.nan,np.nan])
    ScalingDict['VHIFlag']=-1
    
    for i in range(len(BootstrapModels)):
        BootstrapModels[i]['VHI']=np.nan
    
    return ScalingDict


def ExtractRHI_NoErr(R,SD,SDLim):


    if len(R) <=1:
        RHI=R[0]
        RHI_Found=False
        R_indx=0
        
    else:
        RHI,RHI_Found,R_indx=GetSD_Intecept(R,SD,SDLim)
        
    return RHI,RHI_Found,R_indx


def CalcRHI_Old(GalaxyDict,BootstrapModels):

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
    
    
    
    if SDCalcMethod==1:
        Beam=GalaxyDict['ExtendedSDProfile']['BMAJ']
        for RCorr in RCorrArr:
            BU=Beam
            i=1
            while BU > RCorr:
                BU=Beam/float(i)
                i+=1
            RCorr=np.sqrt(RCorr**2.-BU**2.)

    RHDict={'RHI_CorrArr':RCorrArr,'RHIArr':RHIArr,'SDMethod':SDCalcMethod,'RHIFlag':RHIFlag}
    return RHDict



def CalcVHI_Old(GalaxyDict,ScalingDict):

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
        Mid[i],RIndx,VFlag=GetProfilePoint(R,VProf,RUse)
        
        Diffs[i]=np.abs(Mid[i]-VHI)
       
   
    Avg=np.mean(Diffs)
    
    TotErr=np.sqrt(VErr_Local**2.+Avg**2.)
    
    if VHIFlag_Ini==True and RHILimsFlags==0:
        ScalingDict['VHIFlag']=0
    elif VHIFlag_Ini==False and RHILimsFlags==0:
        ScalingDict['VHIFlag']=-1
        dR=R[1]-R[0]
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

def GetVHIFromProf(R,VProf,RHI):

    if len(VProf) <=1:
        VHI=VProf[0]
        RIndx=0
        VHIFlag=-1
    else:
        VHI,RIndx,VHIFlag=GetProfilePoint(R,VProf,RHI)

    return VHI,RIndx,VHIFlag
