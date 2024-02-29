import numpy as np
import os as os
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs


def CubePreprocessing(GalaxyDict):
    print("Generating velocity cube")
    CubeHeaderConvert(GalaxyDict)

def CubeHeaderConvert(ObjDict):
    """
        This function adjusts the header of the SoFiA cubelets so that the output cubes are in velocity space with a reference point in the center of the cubelet.
    """
    #   Get the target name for the converted velocity file.
    OutName=ObjDict['VelCubeName']
    #   Open up the frequency file
    Cubehdu=fits.open(ObjDict['FreqCubeName'])
    #   Read the frequency cube header
    CubeHeader=Cubehdu[0].header
    #   Get the wcs coordinate information for the frequency cube
    w = wcs.WCS(Cubehdu[0].header)
    #   Now to the frequency-velocity conversion
    Cubehdu=Frequency_VelocityConversion(Cubehdu,CubeHeader)
    #   Check that the header has the units for the cube header
    Cubehdu=UnitsCheck(Cubehdu,CubeHeader)
    #   Adjust the cube object name to match the dictionary object name
    Cubehdu[0].header.set('OBJECT',ObjDict['name_underscore'])
    
    #   Write out the cube with the new header to the given file name
    try:
        Cubehdu[0].writeto(OutName)
    except:
        os.system("rm " + OutName)
        Cubehdu[0].writeto(OutName)
    #   Close the file
    Cubehdu.close()



def Frequency_VelocityConversion(CubeHDU,CubeHeader):
    """
        This function converts the spectral axis from frequency to velocity.  It assumes that the units of the frequency axis are Hertz.
        
        Because the source cubes are so large, the size of the velocity channels must be calculated at the central reference point of the cube and not at the overall reference point as dV is not constant across the larger cube.
    """
    #      First get the frequecy cube reference channel, frequency, and bin size
    RefChannel=CubeHDU[0].header['CRPIX3']
    RefFrequency=CubeHDU[0].header['CRVAL3']
    DeltaFrequency=CubeHDU[0].header['CDELT3']
    
    #   Now set the new reference channel to be the middle of the cube
    NewRefChannel=int(CubeHeader['NAXIS3']/2)

    #   Get the number of channels being changed
    DeltaChannel=NewRefChannel-RefChannel

    #   Get the frequency at the new reference channel
    F0=RefFrequency+(DeltaChannel)*DeltaFrequency
    #   Get the frequency at the channel right after the new reference channel
    F1=RefFrequency+(DeltaChannel+1)*DeltaFrequency
    
    #   Set the rest frequency to that of the HI 21 cm line in Hz
    RestFreq=1.42040575179E+09
    #   Convert the reference frequency to a velocity in m/s
    V0=RedShiftConv(F0,RestFreq)
    #   Convert the frequency in the neighbouring channgel to a velocity in m/s
    V1=RedShiftConv(F1,RestFreq)
    #   Use these two velocities to get the velocity width for channels in the cubelet
    dV=V1-V0

    #   Save all the new velocity information to the CubeHDU
    CubeHDU[0].header.set('CRPIX3',NewRefChannel)
    CubeHDU[0].header.set('CRVAL3',V0)
    CubeHDU[0].header.set('CDELT3',dV)
    CubeHDU[0].header.set('CTYPE3','VOPT')
    CubeHDU[0].header.set('CUNIT3','m/s')
    #   return the Cube HDU
    return CubeHDU


def RedShiftConv(Freq,RestFreq):
    """
        Calculate the velocity of some object by it's redshifted spectrum
    """
    #   First get the redshift
    z=(RestFreq-Freq)/Freq
    #   Next get the velocity by multiplying the redshift by the speed of light in m/s.
    Vel=z*2.9979245e8
    return Vel


def UnitsCheck(CubeHDU,CubeHeader):

    keys=['CUNIT1','CUNIT2']
    for key in keys:
        if key not in CubeHeader:
            print("Missing angular unit definition", key)
            print("Assuming that it is degrees")
            CubeHDU[0].header.set(key,'deg')
    return CubeHDU
    
    
