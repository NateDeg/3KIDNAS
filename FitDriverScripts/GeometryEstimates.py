"""
GeometryEstimates.py
=======

**Author:** Nathan Deg

**Description:**
This module defines the functions used to calculate the inclination and position estimates used in 3KIDNAS.  It assumes that the systems have been analyzed with SoFiA for mask making.
"""
import numpy as np
import pandas as pd
import os.path
from os import path
import os as os

def GetGeometryEstimates(SoFiACat,BeamPix,IncMethod):
    """
    This function calculates the PA and Inc estimates used for 3KIDNAS from SoFiA parameters
    
    Parameters
        ----------

        SoFiACat : Dictionary
            The underlying SoFiA catalogue containing all estimates
        BeamPix : float
            The size of the beam in pixels
        IncMethod : string
            The method used for obtaining the inclination estimate

    Returns
        -------
        PAEst : float
            The position angle estimate
        IncEst : float
            The inclination estimate
        
    """
    #   Start by getting the PA estimate
    PAEst=PAEstimate(SoFiACat['kin_pa'],SoFiACat['ell_pa'])
    #   Then get the inclination estimate
    IncEst=IncEstimate(SoFiACat,BeamPix,IncMethod)
    return PAEst,IncEst
    
def IncEstimate(SoFiACat,BeamPix,IncMethod):
    #   There are three built in options for estimating the inclination from the SoFiA catalogue.  They are all based on the SoFiA estimated size and have differing beam corrections
    if IncMethod=='UnCorrected':
        #   The uncorrected version uses the base minor and major axis size estimate.
        MinUse2=SoFiACat['ell_min']**2.
        MajUse2=SoFiACat['ell_maj']**2.
    elif IncMethod=='BeamCorrected':
        #   The second option is to correct for the beam size.  This assumes the ell_min and ell_maj are approximately the diameter of the galaxy.
        MinUse2=BeamCorrSize2(SoFiACat['ell_min'],BeamPix)
        MajUse2=BeamCorrSize2(SoFiACat['ell_maj'],BeamPix)
    elif IncMethod=='WALLABY_Like':
        #   The third option is to correct for the beam size assuming that ell_min and ell_maj are approximately the radii of the galaxy.  This is true for WALLABY-like observations
        MinUse2=BeamCorrSize2(2.*SoFiACat['ell_min'],BeamPix)
        MajUse2=BeamCorrSize2(2.*SoFiACat['ell_maj'],BeamPix)
    else:
        print("No valid inclination estimate method provided.  Ending here")
        exit()

    #   With the corrected sizes, there are a few adjustments that may be needed.
    #   First, make sure the major axis used for the size calculation isn't negative
    if MajUse2 < 0.:
        MajUse2=BeamPix**2.
    #   If the corrected minor axis is unresolved, set the inclination to be edge on.
    if MinUse2 <=0.:
        IncEst=89.
    else:
        IncEst=np.arccos(np.sqrt(MinUse2)/np.sqrt(MajUse2))*180./np.pi

    return IncEst
    
def PAEstimate(kinPA,ellPA):
    """
    This function selects which position angle estimate to use for 3KIDNAS. Generally we always want the kinematic PA, but sometimes SoFiA fails to measure this. In such a case, we use the PA from the moment 0 map
    
    Parameters
        ----------
        kinPA : real
            The kinematic PA estimate from SoFiA for the galaxy
        ellPA : real
            The PA estimate from the galaxy's moment 0 map made by SoFiA

        Returns
        -------
        PA_Estimate: real
            The PA estimate that will be used to initialize the fitting run
    """
    #   Check whether the kinematic PA is a nan.  If it is, use ellPA.  Otherwise use kinPA
    if np.isnan(kinPA):
        PA_Estimate=ellPA
    else:
        PA_Estimate=kinPA
        
    return PA_Estimate

def BeamCorrSize2(Size,Beam):
    """
    This function calculates the beam corrected size
    Parameters
        ----------
        Size : real
            The size to be correct
        ellPA : real
            The beam size in the same unit as Size

        Returns
        -------
        BeamCorrSize2: real
            The Beam corrected size squared
    """
    
    BeamCorrSize2=Size**2.-Beam**2.
    return BeamCorrSize2
