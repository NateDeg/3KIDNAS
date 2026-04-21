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


from scipy.interpolate import interp1d
from scipy.optimize import brentq

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
    RHIArr=np.array([GalaxyDict['BestFitModel']['RHI']-RHI_Err,GalaxyDict['BestFitModel']['RHI'],GalaxyDict['BestFitModel']['RHI']+RHI_Err])
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
    """
    This function finds the outermost X-value where the Y-values of a profile crosses some target limit value.  Note that this function is/should only called when the min and max Y values cross the Lim value.
    
    Parameters
        ----------
        X : real array
            The X points for the array
        Y : real array
            The Y points making up the array
        Lim : real
            The target value for the intersection
    Returns
        ------- 
        X_Int : Real
            The X-value of the profile where Y=Lim
        indx : integer
            The lower index of the two points bracketing Lim.
    """
    #   Start by subtracting the limit from Y
    YNorm=Y-Lim
    #   Find the indices of all spots where YNorm crosses 0
    crossings = np.where(np.diff(np.sign(YNorm)))[0]
    #   Set the index to the outermost crossing point
    if len(crossings)==0:
        indx=-1
        X_Int=np.nan
        return X_Int,indx
    indx=crossings[-1]
    #   Now make a 1D linearly interpolated function for the normalized profile
    f2=interp1d( X, YNorm )
    #   Find the spot where YNorm=0 in the region between the two indices
    X_Int=brentq(f2,X[indx],X[indx+1])
    #   Return the X point and index
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


def GetVHIFromProf(R,VProf,RHI):

    if len(VProf) <=1:
        VHI=VProf[0]
        RIndx=0
        VHIFlag=-1
    else:
        VHI,RIndx,VHIFlag=GetProfilePoint(R,VProf,RHI)

    return VHI,RIndx,VHIFlag
