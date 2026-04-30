"""
ModelAcceptance.py
=======

**Author:** Nathan Deg

**Description:**
This module contains the functions that determine whether a fit is reliable enough to be placed into the accepted catalogue
"""


import numpy as np
import pandas as pd
import os.path
from os import path
import multiprocessing as mp

import subprocess

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs
    
def SetCutLimits():
    """
    This function sets the hard limits for the automated acceptance
    """
    #
    lim_nR=4
    #   The inclination limits are for flagging purposes
    lim_Inc=25.
    lim_IncU=75.
    
    lim_VSysErr=15.
    lim_PAErr=15.
    lim_IncErr=15.
    lim_deltaSinI=0.15
    
    lim_FluxDiff=1.25
    limVRot=75.
    
    limDict=locals()
    return limDict
    
def DetermineSuccess(Model,CutLimits,ModelNames,BeamSize_Pix):
    """
    This function determines whether a particular model will be accepted 
    Parameters:
    -------
    Model : Dictionary
        The dictionary containing all the information about the model
    CutLimits : Dictionary
        The dictionary containing the limits for the accept/reject catalogue
    ModelNames : Dictionary
        The dictionary containing various cube names
    BeamSize_Pix : Float
        The size of the beam in pixels
        
    Returns:
    -------
    Model : Dictionary
        The dictionary containing all the information about the model
    """

    #   Build a model check key that can be used to add to the flags file
    ModelCheckDict={}
    #CheckKeys=['FitAchieved','nRings','size','Inc','VSys_Err','PA_Err','deltaSinI','NaNErrs','VelLims']
    CheckKeys=['FitAchieved','nRings','Inc','Inc_Err','PA_Err','deltaSinI','NaNErrs','VelLims','VRotErrs','FluxCheck','VSys_Err']
    for key in CheckKeys:
        ModelCheckDict[key]=0

    #   Start by assuming success
    AutoSuccess=1
    #   Set the profile keys to check
    ProfKeys=['POSITIONANGLE_ERR','VSYS_ERR','INCLINATION_ERR']
    LimProfKeys=['lim_PAErr','lim_VSysErr','lim_IncErr']
   
    #  Check whether the model actually was fit
    if Model['ModelFitAchieved']==False:
        AutoSuccess=0
        for key in CheckKeys:
            ModelCheckDict[key]=0
        ModelCheckDict['FitAchieved']=1
    else:
        #   Check the length of the model
        nR=len(Model['Model']['R'])
        if nR <CutLimits['lim_nR']:
            AutoSuccess=0
            ModelCheckDict['nRings']=1
        #   If the limits are equal to the number of rings, we do want to flag it, but not remove it
        if nR == CutLimits['lim_nR']:
            ModelCheckDict['nRings']=2
        #   Check on the PA, Inc, and VSys errors
        i=0
        for key in ProfKeys:
            if Model['Model'][key][0] >= CutLimits[LimProfKeys[i]]:
                AutoSuccess=0
                if key =='VSYS_ERR':
                    ModelCheckDict['VSys_Err']=1
                elif key =='POSITIONANGLE_ERR':
                    ModelCheckDict['PA_Err']=1
                elif key == 'INCLINATION_ERR':
                    ModelCheckDict['Inc_Err']=1
            i+=1
        #   Check on the inclination -- note this is not acceptablity cut, but just a check
        if Model['Model']['INCLINATION'][0] <= CutLimits['lim_Inc']:
            ModelCheckDict['Inc']=1
        elif Model['Model']['INCLINATION'][0] >= CutLimits['lim_IncU']:
            ModelCheckDict['Inc']=2
        
    #   It is possible for all the bootstraps to fail and have NaN's in the errors.  Reject any model where the geometric errors are NaNs
        SuccessCheck=CheckGeoErrorForNaNs(Model)
        if SuccessCheck==0:
            AutoSuccess=0
            ModelCheckDict['NaNErrs']=1
    #   Make sure projected velocities are inside the cube
        SuccessCheck=CheckProjectedVel(Model,ModelNames)
        if SuccessCheck==0:
            AutoSuccess=0
            ModelCheckDict['VelLims']=1
        #   Check the median velocity errors
        MedVRotErr=np.median(Model['Model']['VROT_ERR'])
        if MedVRotErr > CutLimits['limVRot']:
           AutoSuccess=0
           ModelCheckDict['VRotErrs']=1
            
        #   Check on the flux
        FluxSuccess=CheckCubeFluxDiff(Model,CutLimits)
        if FluxSuccess==0:
            AutoSuccess=0
            ModelCheckDict['FluxCheck']=1
        
    #   Add the tag to the model
    Model['ModelSuccess']=AutoSuccess
    #   For clarity, go to the flags file and add all the checks to that file --- note that this file will only exist if the fit was achieved
    if Model['ModelFitAchieved']:
        AddChecksToFlagFile(Model,ModelNames,ModelCheckDict,CutLimits)
    
    return Model
    
def CheckProjectedVel(Model,ModelNames):
    #print("Check Velocity Limits",ModelNames['CubeFile'])
    #   Open up the velocity cube
    HDU=fits.open(ModelNames['CubeFile'])
    #   Get the Header
    Header=HDU[0].header
    #   Close the cube
    HDU.close()
    #   Get the Min and Max velocities of the cube
    VMin=(Header['NAXIS3']-Header['CRPIX3'])*Header['CDELT3']+Header['CRVAL3']
    VMax=(-Header['CRPIX3'])*Header['CDELT3']+Header['CRVAL3']
    #   Convert to km/s
    VMin=VMin/1000.
    VMax=VMax/1000.
    if VMax < VMin:
        Temp=VMax
        VMax=VMin
        VMin=Temp
    #   Now get the largest rotation velocity
    VModel=np.nanmax(Model['Model']['VROT'])
    #   Get the inclination in radians
    Inc_Radian=Model['Model']['INCLINATION'][0]*np.pi/180.
    #   Get the projected velocity
    VProj=VModel*np.sin(Inc_Radian)
    #   Get the velocity limits using the systemic velocity and projected velocity
    VMLow=Model['Model']['VSYS'][0]-VProj
    VMHigh=Model['Model']['VSYS'][0]+VProj
    #   Do the low and high velocity checks
    AutoSuccess=1
    
    
    if VMLow < VMin:
        AutoSuccess=0
    if VMHigh > VMax:
        AutoSuccess=0
    #print("VCheck",AutoSuccess,VMin,VMax,VMLow,VMHigh)

    return AutoSuccess
    
def CheckGeoErrorForNaNs(Model):
    #   List the various error keys to check
    ErrorKeys=['INCLINATION_ERR','VSYS_ERR','POSITIONANGLE_ERR','XCENTER_ERR','YCENTER_ERR','VDISP_ERR']
    #   Loop through all the keys
    for key in ErrorKeys:
        ErrTest=Model['Model'][key][0]
        if np.isnan(ErrTest):
            Success=0
            return Success

def CheckCubeFluxDiff(Model,CutLimits):

    AutoSuccess=1
    try:
        #   Open the difference cube and get the data
        DCube=fits.open(Model['DiffCube'])
        DData=DCube[0].data
        DCube.close()
        #   Do the same with the model cube
        Cube=fits.open(Model['ProcCube'])
        CData=Cube[0].data
        Cube.close()
        #   And the mask cube
        Mask=fits.open(Model['OriMask'])
        MData=Mask[0].data
        Mask.close()
        #   Now mask the difference and cube data
        DData=DData*MData
        CData=CData*MData
        #   Normalize the masked difference cube
        NormDiff=1000.*DData/Model['Model']['RMS']
        #   Get the total flux
        CTot=np.nansum(CData)
        #   Get the total number of cells
        nCells=np.nansum(MData)
        #   Get the rms of the normalized difference cube
        DiffTot=np.sqrt(np.nansum(NormDiff**2.)/(nCells))
    except:
        DiffTot=0.
        CTot=1.
        AutoSuccess=0
    #   If the the normalized difference is above the flux limit, cut it
    if np.abs(DiffTot) >=  CutLimits['lim_FluxDiff']:
        AutoSuccess=0
    return AutoSuccess

def AddChecksToFlagFile(Model,ModelNames,ModelCheckDict,CutLimits):
    #   The flags file should contain all the information about which flags the model passed to be accepted or rejected
    #   The Fotran code writes a text file, but we will want to convert this to a csv
    #   Get the name of the original flag file
    OriFlagFile=ModelNames['FlagFile'].rsplit('.',1)[0]+".txt"
    #   Now read in that file
    with open(OriFlagFile, 'r') as file:
        content = file.readlines()
    #   Store the key lines into a dictionary
    FlagDict={}
    FlagDict['BestFitGoodness']=float(content[1])
    FlagDict['RMS']=float(content[3])
    FlagDict['nCells']=int(content[5])
    FlagDict['NormGoodness']=float(content[7])
    #   Now set the additional flags based on the fits
    TargFlags=['FitFlag','nRingsFlag','PAErrFlag','MedianVErrFlag','VSysErrFlag','NormFluxDiffFlag','IncErrFlag','IncFlag','NaNErrFlag']
    MatchedFlags={'FitFlag':'FitAchieved','nRingsFlag':'nRings','PAErrFlag':'PA_Err','MedianVErrFlag':'VRotErrs','VSysErrFlag':'VSys_Err','NormFluxDiffFlag':'FluxCheck','IncErrFlag':'Inc_Err','IncFlag':'Inc','NaNErrFlag':'NaNErrs'}

 
    for key in TargFlags:
        mKey=MatchedFlags[key]
        FlagDict[key]=[ModelCheckDict[mKey]]
        
    DF=pd.DataFrame.from_dict(FlagDict)
    DF.to_csv(ModelNames['FlagFile'],index=False)
    
    #subprocess.run(["rm", OriFlagFile], capture_output=True, text=True)
