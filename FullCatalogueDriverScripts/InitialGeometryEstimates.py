"""
InitialGeometryEstimates.py
=======

**Author:** Nathan Deg

**Description:**
This module defines the functions used to calculate the initial inclination and position estimates used in 3KIDNAS
It assumes that the input catalogue is from SoFiA and has all associated parameters.
"""
import numpy as np
import pandas as pd
import os.path
from os import path
import os as os

def GetGeometryEstimates(step,Cat,GalaxyDict,RTDict):
    """
    This function is actually gets the geometry estimates for 3KIDNAS.
    
    Parameters
        ----------
        step : integer
            The current step in the loop going through the catalogue 
        Cat : Dictionary
            The underlying SoFiA catalogue containing all estimates
        GalaxyDict: Dictionary
            The dictionary containing all parameters related to the specific galaxy about to be analyzed
        RTDict: Dictionary
            The dictionary containing runtime variables

    Returns
        -------
        GalaxyDict : Dictionary
            The updated dictionary now containing the initial geometry estimates.
    """
    #   Start by getting the PA estimate
    GalaxyDict['PA_Estimate']=PAEstimate(Cat['kin_pa'][step],Cat['ell_pa'][step])
    #   There are three built in options for estimating the inclination from the SoFiA catalogue.  They are all based on the SoFiA estimated size and have differing beam corrections
    if RTDict['IncIniMethod']=='UnCorrected':
        #   The uncorrected version uses the base minor and major axis size estimate.
        MinUse2=Cat['ell_min'][step]**2.
        MajUse2=Cat['ell_maj'][step]**2.
    elif RTDict['IncIniMethod']=='BeamCorrected':
        #   The second option is to correct for the beam size.  This assumes the ell_min and ell_maj are approximately the diameter of the galaxy.
        MinUse2=BeamCorrSize2(Cat['ell_min'][step],GalaxyDict['Beam_Pix'])
        MajUse2=BeamCorrSize2(Cat['ell_maj'][step],GalaxyDict['Beam_Pix'])
    elif RTDict['IncIniMethod']=='WALLABY_Like':
        #   The third option is to correct for the beam size assuming that ell_min and ell_maj are approximately the radii of the galaxy.  This is true for WALLABY-like observations
        MinUse2=BeamCorrSize2(2.*Cat['ell_min'][step],GalaxyDict['Beam_Pix'])
        MajUse2=BeamCorrSize2(2.*Cat['ell_maj'][step],GalaxyDict['Beam_Pix'])
    else:
        print("No valid inclination estimate method provided.  Ending here")
        exit()

    #   With the corrected sizes, there are a few adjustments that may be needed.
    #   First, make sure the major axis used for the size calculation isn't negative
    if MajUse2 < 0.:
        MajUse2=GalaxyDict['Beam_Pix']**2.
    #   If the corrected minor axis is unresolved, set the inclination to be edge on.
    if MinUse2 <=0.:
        IncEst=89.
    else:
        IncEst=np.arccos(np.sqrt(MinUse2)/np.sqrt(MajUse2))*180./np.pi
    #   Store the inclination estimate in the GalaxyDict and return it.
    GalaxyDict['IncEst']=IncEst
    return GalaxyDict
    
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
