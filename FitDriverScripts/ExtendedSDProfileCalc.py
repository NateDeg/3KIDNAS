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

from . import CubeAnalysis as CA
from . import MomentMapPlotFncs as MMP

def CalcExtendedSDProfile(GalaxyDict):
    print("Calculating Extended profile")
    #   The first thing to do is to get the moment 0 map from the masked cubelet
    #       Start by loading in the data cube and mask
    DataCube=CA.BasicCubeAnalysis(GalaxyDict['CubeName'])
    MaskCube=CA.BasicCubeAnalysis(GalaxyDict['MaskName'])
    #       Now mask the cubelet
    DataCube['MaskedData']=DataCube['Data']*MaskCube['Data']
    #       Use the MakeMomMap function from MomentMapPlotFncs to make a Mom0 map
    Mom0=MMP.MakeMomMap(DataCube,0,1)
    
    ChanWidth=np.abs(DataCube['CubeHeader']['CDELT3']/1000.)
    BeamSize=DataCube['CubeHeader']['BMAJ']/abs(DataCube['CubeHeader']['CDELT1'])
    BeamArea3=2.*np.pi*(BeamSize/2.355)**2.
    #   Convert to Jy km/s pixel^-1
    Mom0=Mom0*ChanWidth/BeamArea3
    
    
    
    ExtendedSDProfile={}
    ExtendedSDProfile['Mom0']=ExtendedSDProfile

    #   Now let's grab the inclination and position angle
    Model=GalaxyDict['BestFitModel']
    ExtendedSDProfile['Inc']=Model['INCLINATION'][0]*np.pi/180.
    ExtendedSDProfile['PA']=(Model['POSITIONANGLE'][0]+90)*np.pi/180.
    ExtendedSDProfile['Ellipticity']=np.cos(ExtendedSDProfile['Inc'])

    #   The next step is to figure out the outermost radius
    maxR,dR=DetermineRMax(Model,DataCube,ExtendedSDProfile['PA'],Mom0)
    ExtendedSDProfile['maxR']=maxR
    ExtendedSDProfile['dR']=dR
    #   It's necessary to have the radial grid size in pixels rather than arcseconds for some calculations.
    dRPix=dR/(abs(DataCube['CubeHeader']['CDELT1'])*3600.)
    ExtendedSDProfile['dRPix']=dRPix
    #   Make a dictionary for the extend profile
    ExtendedSDProfile['FluxTot']=np.zeros(maxR)
    ExtendedSDProfile['nPixTot']=np.zeros(maxR)
    ExtendedSDProfile['FluxDiff']=np.zeros(maxR)

    #   Set up the radial grid
    ExtendedSDProfile=SetupRadialGrid(ExtendedSDProfile,Model,DataCube)
    #   Get the profile
    ExtendedSDProfile=GetEllipseFittingSDProfile(ExtendedSDProfile,Mom0,Model,DataCube)
    #   Now get the uncertainty on the profile
    ExtendedSDProfile=GetEllipseFittingUncertainties(ExtendedSDProfile,Mom0,Model,DataCube)
    
    
    
    ExtendedSDProfile['SURFDENS']=ProfUnitConvert(ExtendedSDProfile['SURFDENS'])
    ExtendedSDProfile['SURFDENS_ERR']=ProfUnitConvert(ExtendedSDProfile['SURFDENS_ERR'])
    
    print("Extend SD Calc", ExtendedSDProfile['SURFDENS'])
    print("Extend SD Err", ExtendedSDProfile['SURFDENS_ERR'])

    ExtendedSDProfile['SURFDENS_FACEON']=ExtendedSDProfile['SURFDENS']*np.cos(ExtendedSDProfile['Inc'])
    ExtendedSDProfile['SURFDENS_FACEON_ERR']=ExtendedSDProfile['SURFDENS_ERR']*np.cos(ExtendedSDProfile['Inc'])
    
    ExtendedSDProfile['BMAJ']=DataCube['CubeHeader']['BMAJ']*3600.

    GalaxyDict['ExtendedSDProfile']=ExtendedSDProfile
    return GalaxyDict
        
    
        
def ProfUnitConvert(SDs):
    SDConv=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    SDs=SDs/SDConv
    return SDs

def GetEllipseFittingSDProfile(ExtendedSDProfile,Mom0,Model,CubeInfo):

    for i in range(np.shape(Mom0)[0]):
        for j in range(np.shape(Mom0)[1]):
            #       Get the elliptical radial index for the current pixel
            RadIndx=DetermineREllipIndx(i,j,Model,ExtendedSDProfile)
            #   Check if the radius is inside the maximum size
            if RadIndx < ExtendedSDProfile['maxR']:
                #   Check if the Mom0 is a NaN
                if np.isnan(Mom0[i,j]):
                    ExtendedSDProfile['FluxTot'][RadIndx]+=0.
                    #   If it's not, add the moment 0 flux to the calculation
                else:
                    ExtendedSDProfile['FluxTot'][RadIndx]+=Mom0[i,j]
                ExtendedSDProfile['nPixTot'][RadIndx]+=1
                
                
    PixArea=(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)**2.
    ExtendedSDProfile['SURFDENS']=copy.copy(ExtendedSDProfile['FluxTot']/(ExtendedSDProfile['nPixTot']*PixArea))
    
    return ExtendedSDProfile
    
def GetEllipseFittingUncertainties(ExtendedSDProfile,Mom0,Model,CubeInfo):

    FluxAvg=copy.copy(ExtendedSDProfile['FluxTot']/(ExtendedSDProfile['nPixTot']))
    for i in range(np.shape(Mom0)[0]):
        for j in range(np.shape(Mom0)[1]):
            #       Get the elliptical radial index for the current pixel
            RadIndx=DetermineREllipIndx(i,j,Model,ExtendedSDProfile)
            #   Check if the radius is inside the maximum size
            if RadIndx < ExtendedSDProfile['maxR']:
                #   Check if the Mom0 is a NaN
                if np.isnan(Mom0[i,j]):
                    ExtendedSDProfile['FluxDiff'][RadIndx]+=(0.-FluxAvg[RadIndx])**2.
                    #   If it's not, add the square of the difference between the pixel flux and average radial flux
                else:
                    ExtendedSDProfile['FluxDiff'][RadIndx]+=(Mom0[i,j]-FluxAvg[RadIndx])**2.
    #       First get the number of beams in the ring -- since dRPix is half a beam, get the beam area via 2*pi*(FWHM/2.355)**2.
    PixArea=(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)**2.
    BeamArea_Pix=(2.*np.pi*(ExtendedSDProfile['dRPix']*2/2.355)**2.)
    #       Now get the number of beams in the ring
    nBeams=ExtendedSDProfile['nPixTot']/BeamArea_Pix
    #       Make sure that the number of beams per ring is at least 1
    for i in  range(len(nBeams)):
        if nBeams[i] <=1:
            nBeams[i]=1
    #       Get the st.dev. as the variance/sqrt(number of beams)
    ExtendedSDProfile['SURFDENS_ERR']=copy.copy(np.sqrt(ExtendedSDProfile['FluxDiff']/(ExtendedSDProfile['nPixTot']))/PixArea)/np.sqrt(nBeams)
    return ExtendedSDProfile
            
            

def SetupRadialGrid(SDProfile,Model,CubeInfo):
    #   Initialize the radial grid to zeros
    SDProfile['R_SD']=np.zeros(SDProfile['maxR'])
    #   Now set up the base grid
    for i in range(SDProfile['maxR']):
        SDProfile['R_SD'][i]=Model['R'][0]+i*SDProfile['dR']
        
    #   Set the first midpoint between the inner circle and the first ring in pixels as this is necessary for indexing
    if len(Model['R']) <=1 :
        SDProfile['RMid1']=SDProfile['dR']
    else:
        SDProfile['RMid1']=Model['R'][1]+SDProfile['dR']/2.
    #   And, as usual, this requires a value in pixels
    SDProfile['RMid1Pix']=SDProfile['RMid1']/(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
    return SDProfile


def DetermineREllipIndx(i,j,AvgModel,ExtendedSDProfile):
    """
        This function determines the index of a point in elliptical coordinates on the surface density radial grid.
    """
    #   Turn the pixel coordinates into X and Y
    Y=copy.copy(i)
    X=copy.copy(j)
    PARad=ExtendedSDProfile['PA']
    Ellipticity=ExtendedSDProfile['Ellipticity']
    #   Get the elliptical radius for this point
    REllip=CalculateEllipticalRadius(X,Y,AvgModel,PARad,Ellipticity)
    #   Get the index for this bin using the size of the bins in pixels
    RealRadBin=(REllip-ExtendedSDProfile['RMid1Pix'])/ExtendedSDProfile['dRPix']
    #   Check whether we managed to get an index less than zero
    if RealRadBin <= 0.:
        RadIndx=0
    else:
        #   Convert the real into an integer
        RadIndx=int((REllip-ExtendedSDProfile['RMid1Pix'])/ExtendedSDProfile['dRPix'])+1
    #   Return the radial index
    return RadIndx
    
def CalculateEllipticalRadius(X,Y,AvgModel,PARad,Ellipticity):
    """
        This function gets the elliptical radius of a point given a specific position angle and ellipticity
    """
    #   First get the pixel position relative to the model center
    XCent=X-AvgModel['XCENTER'][0]
    YCent=Y-AvgModel['YCENTER'][0]
    #   Next rotate the points using the PA in radians
    XRot=XCent*np.cos(-PARad)-YCent*np.sin(-PARad)
    YRot=XCent*np.sin(-PARad)+YCent*np.cos(-PARad)
    #   Get the Yprime coorddinate by divinding by the ellipticity
    YEllip=YRot/Ellipticity
    #   Calculate the elliptical radius
    REllip=np.sqrt(XRot**2.+YEllip**2.)
    return REllip

def DetermineRMax(AvgModel,CubeInfo,PARad,Mom0):
    """
        This function gets the maximum radius for the surface density profile
    """
        #   Get the radial size of the bins in pixels
    #dR=AvgModel['R'][1]-AvgModel['R'][0]
    dR=AvgModel['R'][0]*2.
    print("DR check", dR)
    
    #       Set the maximum radius by setting it to the length of the diaganol across the image
    maxX=(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)*CubeInfo['CubeHeader']['NAXIS1']
    maxY=(abs(CubeInfo['CubeHeader']['CDELT2'])*3600.)*CubeInfo['CubeHeader']['NAXIS2']
    maxR=int(np.sqrt(maxX**2.+maxY**2)/dR)
    
    #   Loop through the maximum possible length of the grid
    for i in range(maxR):
        #   Get the radius of the current point
        TestR=(AvgModel['R'][0]+(i+0.5)*dR)
        #   Convert it to a point in pixels
        TestR=TestR/(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
        #   Turn it to a point in X and Y for the pixels
        X,Y=GetPixelPoint(TestR,AvgModel,PARad)
        #   Only start checking the limits once we reach the end of the existing rotation curve grid
        if i >= len(AvgModel['R'])-1:
            #   Check whether X and Y are inside the image
            if X <0 or X >CubeInfo['CubeHeader']['NAXIS1']-1:
                print("Beyond X Lims SD break")
                break
            if Y <0 or Y >CubeInfo['CubeHeader']['NAXIS2']-1:
                print("Beyond Y Lims SD break")
                break
            #   If it's inside, get the indices of the X and Y values
            iTest=int(Y)
            jTest=int(X)
            #   Check whether the map has been masked out
            if Mom0[iTest,jTest] <=0.:
                print("Flux Lims SD break", X,Y,Mom0[iTest,jTest])
                break
    #   Once the loop ends (either from hitting the masked out points or the edge of the profile, i+1 will be the maximum number of radial points to use in the SD profile calculation.
    maxR=i
    #   Return both the maximum number of grid points and the size of the grid.
    return maxR,dR
    
    
def GetPixelPoint(R,AvgModel,PARad):
    """
        This function turn a radius + position angle into a pixel location
    """
    #   In elliptical coordinates, the point is at (R,0)
    XEllip=R
    YEllip=0.
    #   Rotate this to physical coordinates
    XRot=XEllip*np.cos(PARad)-YEllip*np.sin(PARad)
    YRot=XEllip*np.sin(PARad)+YEllip*np.cos(PARad)
    #   Go from coordinates relative to the center to absolut pixel coordinates.
    XPos=XRot+AvgModel['XCENTER'][0]
    YPos=YRot+AvgModel['YCENTER'][0]
    #   Return the position.
    return XPos,YPos
