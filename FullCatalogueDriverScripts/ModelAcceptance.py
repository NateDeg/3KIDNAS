import numpy as np
import pandas as pd
import os.path
from os import path
import multiprocessing as mp

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs
    
def SetCutLimits():
    #   Set the limits for the automated acceptance
    lim_nR=2
    lim_Inc=25.
    lim_Size=1.5
    #lim_Size*=5.    #Convert the size limit into pixels
    lim_VSysErr=40.
    lim_PAErr=20.
    lim_deltaSinI=0.15
    
    limDict=locals()
    return limDict
    
def DetermineSuccess(Model,CutLimits,ModelNames,BeamSize_Pix):

    #   Build a model check key that can be used to add to the flags file
    ModelCheckDict={}
    CheckKeys=['FitAchieved','nRings','size','Inc','VSys_Err','PA_Err','deltaSinI','NaNErrs','VelLims']
    for key in CheckKeys:
        ModelCheckDict[key]=1

    #   Start by assuming success
    AutoSuccess=1
    #   Set the profile keys to check
    ProfKeys=['VSYS_ERR','POSITIONANGLE_ERR']
    LimProfKeys=['lim_VSysErr','lim_PAErr']
   
    #  Check whether the model actually was fit
    if Model['ModelFitAchieved']==False:
        AutoSuccess=0
        for key in CheckKeys:
            ModelCheckDict[key]=-1
        ModelCheckDict['FitAchieved']=0
    else:
        #   Check the length of the model
        nR=len(Model['Model']['R'])
        if nR <=CutLimits['lim_nR']:
            AutoSuccess=0
            ModelCheckDict['nRings']=0
        #   Check on the size
        if Model['ell_maj_SoFiA'] <= CutLimits['lim_Size']*BeamSize_Pix:
            AutoSuccess=0
            ModelCheckDict['size']=0
        #   Check on the VSys Error, and PA Error
        i=0
        for key in ProfKeys:
            if Model['Model'][key][0] >= CutLimits[LimProfKeys[i]]:
                AutoSuccess=0
                if key =='VSYS_ERR':
                    ModelCheckDict['VSys_Err']=0
                elif key =='POSITIONANGLE_ERR':
                    ModelCheckDict['PA_Err']=0
            #print(key,Model['Model'][key][0],AutoSuccess)
            i+=1
        #   Check on the inclination
        if Model['Model']['INCLINATION'][0] <= CutLimits['lim_Inc']:
            AutoSuccess=0
            ModelCheckDict['Inc']=0
        #   Check on the delta Sin I
        ILow=(Model['Model']['INCLINATION'][0]-Model['Model']['INCLINATION_ERR'][0])*np.pi/180.
        if ILow < 0.:
            ILow=0.
        IHigh=(Model['Model']['INCLINATION'][0]+Model['Model']['INCLINATION_ERR'][0])*np.pi/180.
        if IHigh > np.pi*0.5:
            IHigh=np.pi*0.5
        deltaSinI=np.sin(IHigh)-np.sin(ILow)
        if deltaSinI >=CutLimits['lim_deltaSinI']:
            AutoSuccess=0
            ModelCheckDict['deltaSinI']=0
    #   It is possible for all the bootstraps to fail and have NaN's in the errors.  Reject any model where the geometric errors are NaNs
        SuccessCheck=CheckGeoErrorForNaNs(Model)
        if SuccessCheck==0:
            AutoSuccess=0
            ModelCheckDict['NaNErrs']=0
    #   Make sure projected velocities are inside the cube
        SuccessCheck=CheckProjectedVel(Model,ModelNames)
        if SuccessCheck==0:
            AutoSuccess=0
            ModelCheckDict['VelLims']=0
        
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

def AddChecksToFlagFile(Model,ModelNames,ModelCheckDict,CutLimits):
    #   The flags file should contain all the information about which flags the model passed to be accepted or rejected
    #   Open the flags file
    Flags = open(ModelNames['FlagFile'], "a")
    #   State whether a model was produced at all during the initial fit
    Str="Model Acceptence/Rejection Checks (0=fail, 1=Success)\n"
    Str+="Initial Pipeline Fit Achieved \n"
    Str+=str(ModelCheckDict['FitAchieved'])+"\n"
    #   State whether there are enough rings
    Str+="The number of rings in the fits is >="+str(CutLimits['lim_nR'])+"\n"
    Str+=str(ModelCheckDict['nRings'])+"\n"
    #   State whether the estimate size is above the limit
    Str+="The SoFiA estimated size (ell_maj) in beams is >="+str(CutLimits['lim_Size'])+"\n"
    Str+=str(ModelCheckDict['size'])+"\n"
    #   State whether the inclination is acceptable
    Str+="The model inclination is >="+str(CutLimits['lim_Inc'])+"\n"
    Str+=str(ModelCheckDict['Inc'])+"\n"
    #   State whether the uncertainty in the inclination is acceptable
    Str+="sin(i_max)-sin(i_min) <="+str(CutLimits['lim_deltaSinI'])+"\n"
    Str+=str(ModelCheckDict['deltaSinI'])+"\n"
    #   State whether the error on the systemic velocity is acceptable
    Str+="The error on the systemic velocity is <="+str(CutLimits['lim_VSysErr'])+"\n"
    Str+=str(ModelCheckDict['VSys_Err'])+"\n"
    #   State whether the error on the position angle is acceptable
    Str+="The error on the position angle is <="+str(CutLimits['lim_PAErr'])+"\n"
    Str+=str(ModelCheckDict['PA_Err'])+"\n"
    #   Make sure all the errors are real
    Str+="No geometric error terms are NaN's\n"
    Str+=str(ModelCheckDict['NaNErrs'])+"\n"
    #   And that the projected RC falls within the curve
    Str+="The projected model rotation curve falls within the datacube\n"
    Str+=str(ModelCheckDict['VelLims'])
    
    Flags.write(Str)
    Flags.close()
    

