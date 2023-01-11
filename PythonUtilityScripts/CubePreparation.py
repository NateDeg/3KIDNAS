import numpy as np
import argparse
import os

import astropy
from astropy.io import fits
from astropy import wcs


def RedShiftConv(Freq,RestFreq):
    z=(RestFreq-Freq)/Freq
    Vel=z*2.9979245e8
    return Vel

def WriteFitsFile(FITShdu,FileName):
    try:
        FITShdu.writeto(FileName)
    except:
        os.system("rm " + FileName)
        FITShdu.writeto(FileName)

def MomMapHeaderModification(CubeHeader,MomentMap,USwitch):
    MomentMap[0].header['CRPIX1']=CubeHeader['CRPIX1']
    MomentMap[0].header['CRVAL1']=CubeHeader['CRVAL1']
    MomentMap[0].header['CRPIX2']=CubeHeader['CRPIX2']
    MomentMap[0].header['CRVAL2']=CubeHeader['CRVAL2']
    if USwitch == 1:
        MomentMap[0].header.set('BUNIT','m/s')
    elif USwitch == 0:
        MomentMap[0].header.set('BUNIT','Jy/pixel')
    return MomentMap

def MomMapFreqConversion(Mom1,Mom2,RestFreq):
    Mom1Data=Mom1[0].data
    Mom1Vels=RedShiftConv(Mom1Data,RestFreq)
    Mom2DeltaFreq=Mom1Data-Mom2[0].data  #subtraction due to inverse frequency-velocity relation
    Mom2TotVels=RedShiftConv(Mom2DeltaFreq,RestFreq)
    Mom2Vels=Mom2TotVels-Mom1Vels
    Mom1[0].data=Mom1Vels
    Mom2[0].data=Mom2Vels
    return Mom1,Mom2



def CubeConversions(BasePath,RA,DEC):

                    #       Set the names of the files that need conversion
    Cube=BasePath+"cube.fits"
    Mom0=BasePath+"moment0.fits"
    Mom1=BasePath+"moment1.fits"
    Mom2=BasePath+"moment2.fits"
    Mask=BasePath+"mask.fits"

    #   Open up the main data cube

    Cube=fits.open(Cube)
    #   Get the header
    CubeHeader=Cube[0].header

    #       Get the Beam Area
    BeamMajPix=Cube[0].header['BMAJ']/np.abs(Cube[0].header['CDELT1'])
    BeamMinPix=Cube[0].header['BMAJ']/np.abs(Cube[0].header['CDELT1'])
    BeamArea_Pixels=((2.*np.pi))*BeamMajPix*BeamMinPix
    #print("Beam Area Check", BeamMajPix,BeamMinPix,BeamArea_Pixels)


    #       Get the coordinates for the first pixel
    w = wcs.WCS(CubeHeader)
    TestPix=[[1,1,1]]
    TestCoords=w.wcs_pix2world(TestPix,1)
    
    #Get the pixel coordinates of the SoFiA center
    #   (This is passed by the parent Fortran pipeline
    TestCent=[[RA,DEC,3000]]
    print(TestCent)
    TestCentPix=w.wcs_world2pix(TestCent,1)
#    print("Test Center Pixel", TestCentPix)
#    print(CubeHeader['CRPIX1'],CubeHeader['CRVAL1'])
#    print(TestCentPix[0][0]-CubeHeader['CRPIX1'])

    #       Set the central RA and DEC as the reference pixels
    #       Convert the spatial reference values
    CubeHeader.set('CRPIX1',TestCentPix[0][0])
    CubeHeader.set('CRVAL1',RA)
    CubeHeader.set('CRPIX2',TestCentPix[0][1])
    CubeHeader.set('CRVAL2',DEC)

    #       Now do a velocity conversion
    RefChannel=CubeHeader['CRPIX3']
    RefFrequency=CubeHeader['CRVAL3']
    DeltaFrequency=CubeHeader['CDELT3']
    F0=RefFrequency-((RefChannel-1))*DeltaFrequency
    F1=RefFrequency-((RefChannel-1)-1)*DeltaFrequency

    RestFreq=1.42040575179E+09
    V0=RedShiftConv(F0,RestFreq)
    V1=RedShiftConv(F1,RestFreq)
    dV=V1-V0

    CubeHeader.set('CRPIX3',0)
    CubeHeader.set('CRVAL3',V0)
    CubeHeader.set('CDELT3',dV)
    CubeHeader.set('CTYPE3','VELOHEL')
    CubeHeader.set('CUNIT3','m/s')
    CubeHeader.set('CUNIT1','DEGREE')
    CubeHeader.set('CUNIT2','DEGREE')
#    CubeHeader.set('BUNIT','Jy/pixel)
    


#   Save the cube to a temporary file
    ConvertedCubeName="TempDataCube.fits"
    WriteFitsFile(Cube,ConvertedCubeName)

#   Deal with the other useful quantities
    Mom0=fits.open(Mom0)
    Mom1=fits.open(Mom1)
    Mom2=fits.open(Mom2)
    Mask=fits.open(Mask)


#   Convert the mask header to the cubes and make a temporary mask
    Mask[0].header=CubeHeader
    ConvertedMaskName="TempMaskCube.fits"
    WriteFitsFile(Mask,ConvertedMaskName)
    
    print("Moment 0 flux units ini", Mom0[0].header['BUNIT'])
    Mom0[0].data /= DeltaFrequency

    
    Mom0[0].data /= BeamArea_Pixels
    print("Beam area in pixels",BeamMajPix, BeamArea_Pixels)

#   Modify the moment map headers
    Mom0=MomMapHeaderModification(CubeHeader,Mom0,0)
    Mom1=MomMapHeaderModification(CubeHeader,Mom1,1)
    Mom2=MomMapHeaderModification(CubeHeader,Mom2,1)
#   Do the moment 1 and 2 conversions
    Mom1,Mom2=MomMapFreqConversion(Mom1,Mom2,RestFreq)

    ConvertedMomNames=["TempMom0.fits","TempMom1.fits","TempMom2.fits"]
    WriteFitsFile(Mom0,ConvertedMomNames[0])
    WriteFitsFile(Mom1,ConvertedMomNames[1])
    WriteFitsFile(Mom2,ConvertedMomNames[2])


#       The main conversion driver
#           Get the basepath name
parser = argparse.ArgumentParser(description='Convert SoFiA cubes')
parser.add_argument('BasePath', metavar='Name',help='some string')
parser.add_argument('RA', metavar='Name',help='some string')
parser.add_argument('DEC', metavar='Name',help='some string')
args = parser.parse_args()
BasePath=args.BasePath
RA=float(args.RA)
DEC=float(args.DEC)
#print("Cube prep", RA,DEC)

CubeConversions(BasePath,RA,DEC)


