import numpy as np
from decimal import Decimal

import os.path
from os import path
import re

def LoadBaroloModel(FileName,DensFile):
    print("Loading Barolo Model")
    BaroloParams,FitAchieved=LoadBaroloParamFile(FileName)
    if FitAchieved:
        SDParams,FitAchieved=LoadBaroloSurfDens(DensFile)
    #   If the fit failed, set the bad fit dictionary
    if FitAchieved ==False:
        BaroloDict=BadBaroloFit()
    #   If the fit was successful, then there is a number of steps that need to happen
    else:
        BaroloDict=BaroloParamstoDict(BaroloParams,SDParams)
        BaroloDict['FitForAveraging']=True
    return BaroloDict
    
    
    
    
def LoadBaroloParamFile(ParamFile):
    FitAchieved=True
    try:
        Params=np.loadtxt(ParamFile,skiprows=1)
    except:
        FitAchieved=False
        Params=[]
    #   Check that the parameters were loaded correctly by making sure there is some shape to them
    if np.shape(Params)[0] == 0:
        FitAchieved=False
    return Params,FitAchieved
    
    
def LoadBaroloSurfDens(DensFile):
    try:
        SDParams=np.loadtxt(DensFile,skiprows=13,usecols=(0,9,10,11),unpack='True')
        FitAchieved=True
    except:
        FitAchieved=False
        SDParams=[]
    return SDParams,FitAchieved
    
    
def BaroloParamstoDict(Params,SDParams):
    #   First set the dictionary according to the parameters
    TRParamsDict={'R':Params[:,1],'R_SD':SDParams[0,:],'XCENTER':Params[:,9] \
            ,'YCENTER':Params[:,10], 'INCLINATION':Params[:,4] \
            ,'POSITIONANGLE':Params[:,5], 'VSYS':Params[:,11], 'VROT':Params[:,2]\
            ,'VRAD':Params[:,12], 'VDISPERSION':Params[:,3],'Z0':Params[:,7]\
                ,'SURFDENS':SDParams[1,:],'SURFDENS_ERR':SDParams[2,:],'FITAchieved':True}
    #   The face-on SD is needed for generating MCG models later on
    TRParamsDict['SURFDENS_FACEON']=SDParams[3,:]
    TRParamsDict['SURFDENS_FACEON_ERR']=SDParams[2,:]
    #   Try to error terms if they exist
    
        #       Get a rotation error if it exists
    try:
        VerrLow=Params[:,13]
        VerrHigh=Params[:,14]
    except:
        VerrLow=Params[:,1]*0.
        VerrHigh=Params[:,1]*0.
        #       Get the inclination error if it exists
    try:
        IncerrLow=Params[:,15]
        IncerrHigh=Params[:,16]
    except:
        IncerrLow=Params[:,1]*0.
        IncerrHigh=Params[:,1]*0.
        #       Get the position angle error if it exists
    try:
        PAerrLow=Params[:,17]
        PAerrHigh=Params[:,18]
    except:
        PAerrLow=Params[:,1]*0.
        PAerrHigh=Params[:,1]*0.

    #       Place the errors into the TR dictionary
    VErr=[VerrLow,VerrHigh]
    TRParamsDict['VROT_ERR']=VErr
    IncErr=[IncerrLow,IncerrHigh]
    TRParamsDict['INCLINATION_ERR']=IncErr
    PAErr=[PAerrLow,PAerrHigh]
    TRParamsDict['POSITIONANGLE_ERR']=PAErr
            
    return TRParamsDict
    
    
def BadBaroloFit():
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':False,'RESID':[]
    ,'FitForAveraging':False,'CHI2':[]}
    return Disk
