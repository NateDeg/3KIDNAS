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
import sys as sys

from . import GenerateAcceptedModelCatalogue as GAMC


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
    #   The key functions are actually located in the FitDriverScripts location.  So the first step is getting that location
    ModDict=GAMC.GetDriverFolder()
    sys.path.append(ModDict['GalFitterDir'])
    #   Now that the path has been appended, we can import the geometry estimates.
    import GeometryEstimates as GE
    #   The geometry estimates function in GE uses a somewhat simpler dictionary, so we'll set that here.
    TempDict={'kin_pa':Cat['kin_pa'][step],'ell_pa':Cat['ell_pa'][step],'ell_min':Cat['ell_min'][step],'ell_maj':Cat['ell_maj'][step]}
    #   Now call the geometry estimate function.
    GalaxyDict['PA_Estimate'],GalaxyDict['IncEst']=GE.GetGeometryEstimates(TempDict,GalaxyDict['Beam_Pix'],RTDict['IncIniMethod'])
    #   And return the dictionary.
    return GalaxyDict
  
