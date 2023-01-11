#!/usr/bin/env python3
import numpy as np
import os

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

def CubeCompare(ObsData,CubeFile):
    ModelCube=fits.open(CubeFile)
    ModelData=ModelCube[0].data
    #print(np.shape(ModelData),np.shape(ObsData['Mask']))
    MaskedModelData=ModelData*ObsData['Mask']
    # print("First check", np.nansum(ObsData['MaskedCube']),np.nansum(MaskedModelData))
    ResidData=(ObsData['MaskedCube']-MaskedModelData)**2./ObsData['Noise']**2.
    ResidTot=np.nansum(ResidData)
    nCells=np.nansum(ObsData['Mask'])
    #    print(ResidTot,nCells,ResidTot/nCells)
    reducedChi2=ResidTot/nCells
    print("reduced chi2",reducedChi2)
    return ResidTot,reducedChi2
