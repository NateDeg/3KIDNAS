import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs


def CubeAnalysisFuncDict():
    CubeFuncs={'SoFiA_CubeAnalysis':SoFiA_CubeAnalysis,'TotalCubeHI':TotalCubeHI,'ConstructPVDiagram':ConstructPVDiagram,'MassCalcFromCube':MassCalcFromCube,'GetVelsFromHeader':GetVelsFromHeader,'GetStartPixelLoc':GetStartPixelLoc,'BasicCubeLoad':BasicCubeLoad,'BasicCubeAnalysis':BasicCubeAnalysis}
    
    return CubeFuncs


    
def BasicCubeLoad(CubeName):
    #   Open up the cube
    CubeHDU=fits.open(CubeName)
    #   Load in the header
    CubeHeader=CubeHDU[0].header
    #   Load in the data
    CubeData=CubeHDU[0].data
    #   Get the wcs values
    CubeWCS = wcs.WCS(CubeHeader)
    #   Close the cube
    CubeHDU.close()
    return CubeHeader,CubeData,CubeWCS
    
def BasicCubeAnalysis(CubeName):
    #   Load in the cube header
    CubeHeader,CubeData,CubeWCS=BasicCubeLoad(CubeName)
    #   Get the cube velocities
    CubeVels=np.array(GetVelsFromHeader(CubeHeader))
    #   Get the start location of the pixels from the header
    CubeHeader=GetStartPixelLoc(CubeWCS,CubeHeader)
    #   Store all the info into a dictionary and return in
    CubeInfo={'Data':CubeData,'CubeHeader':CubeHeader,'CubeWCS':CubeWCS,'CubeVels':CubeVels}
    return CubeInfo
    
def GetStartPixelLoc(CubeWCS,CubeHeader):
    pixcrd=[[0,0,0]]
    RealCoord=CubeWCS.wcs_pix2world(pixcrd,0)
    print("Get pixel start location, real coords--", RealCoord)
    CubeHeader['CRLOC1']=RealCoord[0,0]
    CubeHeader['CRLOC2']=RealCoord[0,1]
    return CubeHeader

    
    
def TotalCubeHI(cube,PixSize,BeamSize):
    cubeTot=np.nansum(cube)         #Get the total signal in cube units (Jy/beam)
    BeamArea=np.pi*BeamSize**2./(4.*np.log(2.)) #Get the beam area in degrees
    pixperbeam=BeamArea/PixSize**2. #Find the number of pixels per beam
    HISignal=cubeTot/pixperbeam #Get the integrated HI signal in Jy
    return HISignal


def MassCalcFromCube(cube,Distance,ChannelSize,BeamSize,PixSize):
    #   Get the total signal in Jy/beam
    Signal_JyBeam=np.nansum(cube)
    #       Get the beam size in degrees
    beamarea=(np.pi*BeamSize**2.)/(4.*np.log(2.))
    #   Ge the number of pixels/beam
    pixperbeam=beamarea/(abs(PixSize)*abs(PixSize))
    #   Get the signal in Jy
    Signal_Jy=Signal_JyBeam/pixperbeam
    #print("Signal in Jy", Signal_Jy,ChannelSize)
    #   Get the mass
    Mass=0.236*(Distance*1000.)**2.*Signal_Jy*abs(ChannelSize)  #Original version
    return Mass


def GetVelsFromHeader(CubeHeader):
    Vels=[]
    for i in range(CubeHeader['NAXIS3']):
        VTemp=float(i-CubeHeader['CRPIX3']+1)*CubeHeader['CDELT3']+CubeHeader['CRVAL3']
        Vels.append(VTemp)
    return Vels


def ConstructPVDiagram(CubeData,Angle,PixelCenter,BeamSize,Vels):
    AngUse=(Angle-90.)*np.pi/180.
    if Angle+90.>360.:
        AngUse=AngUse-2.*np.pi
    if Angle-90.<0.:
        AngUse=AngUse+2.*np.pi
        
    #   Figure out the max size of the PV diagram
    CubeDim=GetPVDiagramDimensions(CubeData,PixelCenter,Vels)
    
    PV=np.zeros(CubeDim)
    

    for i in range(np.shape(CubeData)[1]):
        for j in range(np.shape(CubeData)[2]):
            X=(j-PixelCenter[0])
            Y=(i-PixelCenter[1])
            XP=X*np.cos(-AngUse)-Y*np.sin(-AngUse)
            YP=X*np.sin(-AngUse)+Y*np.cos(-AngUse)
            
            k=int(round(XP+PixelCenter[0]))
            l=int(round(YP+PixelCenter[1]))

            if k >= 0 and k < CubeDim[0]:
                if abs(YP) <= BeamSize/2.:
                    for m in range(np.shape(CubeData)[0]):
                        if np.isnan(CubeData[m,i,j]):
                            PV[k,m]+=0.
                        else:
                            PV[k,m]+=CubeData[m,i,j]
                            
    #print("PVSum check", np.sum(PV))
    return PV

def GetPVDiagramDimensions(CubeData,PixelCenter,Vels):
    #print("Getting PV diagram size")
    Size=0.
    #   Loop through x and y dimensions and figure out the largest difference in pixels from the center to the edge
    for i in range(2):
        #   The j index is due to cube data dimensions being inverted
        j=2-i
        Size=np.max([Size,PixelCenter[i],np.shape(CubeData)[j]-PixelCenter[i]])
    PixSize=int(2*Size)+1
    CubeDim=(PixSize,np.shape(Vels)[0])
    return CubeDim
