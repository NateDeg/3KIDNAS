import sys as sys
import os as os
import copy as copy
import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs


    
def GetGlobalPositionAngle(GalaxyDict):
    """
       The kinematic modelling codes get the position angle relative to the pixel axes, but these may be tilted with respect to RA and DEC.  This routine gets the position angle in the RA and DEC coordinates using the astropy wcslib functions.
    """
    #Cube=fits.open(GalaxyDict[])
    #print(GalaxyDict['CubeHeader'])
    CubeHeader=GalaxyDict['CubeHeader']
    CubeWCS = wcs.WCS(CubeHeader)
    
    #   Set the location of the center in pixels as well as the average model RA and DEC
    Model=GalaxyDict['BestFitModel']
    #print("Model Check", Model)
    XPix=Model['XCENTER'][0]
    YPix=Model['YCENTER'][0]
    RA_Cent=Model['RA'][0]
    DEC_Cent=Model['DEC'][0]
    
    

    #   Set a new DEC value one beam higher
    DEC_Delt=DEC_Cent+CubeHeader['BMAJ']
    #   Set the RA_Cent and DEC_Delt into a coordinate
    CoordDelt=[[RA_Cent,DEC_Delt,0]]
    #   Convert this point to pixel coordinates.
    PosDelt=CubeWCS.wcs_world2pix(CoordDelt,0)
    #   Now figure out how far away this new point is from the center point
    XDelt=PosDelt[0][0]-XPix
    YDelt=PosDelt[0][1]-YPix
    #   Use these points to get the angle the declination makes with the X-Y axis
    DEC_Angle=np.arctan2(YDelt,XDelt)*180./np.pi
    #   Get the global PA by adjusting the XY version of PA by the angle DEC makes with the XY axis
    PA_Global=Model['POSITIONANGLE'][0]-(DEC_Angle-90.)
    #   Ensure that the PA is between 0 and 360
    if PA_Global < 0.:
        PA_Global+=360.
    elif PA_Global > 360.:
        PA_Global-=360.
    #   Store the global PA and associated error into the tilted ring dictionary
    Model['PA_GLOBAL']=PA_Global
    Model['PA_GLOBAL_ERR']=Model['POSITIONANGLE_ERR'][0]
    #   Round the results
    Model['PA_GLOBAL']=round(Model['PA_GLOBAL'],1)
    Model['PA_GLOBAL_ERR']=round(Model['PA_GLOBAL_ERR'],1)
    return Model
 
