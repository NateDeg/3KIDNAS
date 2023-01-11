#!/usr/bin/env python3
import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

def CubeAnalysis(CatID,StrDict,WallCat):
    Name=WallCat.name[CatID]
    NumStr=Name.split(' ')[1]
    CubeName=StrDict['VelCubePath'][0]+"/"+StrDict['BaseFileName']+"_"+NumStr+"_cube_Vel.fits"
    MaskName=StrDict['DataFolder']+StrDict['BaseFileName']+"_source_products/"+StrDict['BaseFileName']+"_"+NumStr+"_mask.fits"
    CubeHDU=fits.open(CubeName)
    MaskHDU=fits.open(MaskName)
    
    CubeHeader=CubeHDU[0].header
    CubeData=CubeHDU[0].data
    MaskData=MaskHDU[0].data
    
    CubeWCS = wcs.WCS(CubeHeader)
    
    MaskedCubeData=CubeData*MaskData
    
    print("Data cube pixel size", CubeHeader['CDELT1'])
    
    BeamSize=CubeHeader['BMAJ']
    PixSize=CubeHeader['CDELT2']
    ChannelSize=CubeHeader['CDELT3']
    ChannelSize=ChannelSize/1000.
    
    
    Distance=GetDistance(CatID,WallCat)
    IntegratedHI=TotalHI(MaskedCubeData,PixSize,BeamSize)
    #CubeMass=MassCalc(IntegratedHI,Distance,ChannelSize)
    CubeMass=MassCalc(MaskedCubeData,Distance,ChannelSize,BeamSize,PixSize)

    CubeHeader=GetStartPixelLoc(CubeWCS,CubeHeader)
    
    print("kin pa",WallCat.kin_pa[CatID]+180.)
    
    CubeVels=np.array(GetVels(CubeHeader))
    
    CenterPos=[[WallCat.ra[CatID],WallCat.dec[CatID],0]]
    CentPix=CubeWCS.wcs_world2pix(CenterPos,0)
    CentPix=[CentPix[0,0],CentPix[0,1]]
    
    CubePV=CalcPVDiagram(WallCat.kin_pa[CatID],BeamSize/PixSize,MaskedCubeData,CubeVels,CentPix)

    CubeHDU.close()
    MaskHDU.close()
    
    RMS=WallCat.rms[CatID]

    CubeInfo={'Mask':MaskData,'MaskedCube':MaskedCubeData\
        ,'CubeHeader':CubeHeader, 'Distance':Distance,'IntegratedHI':IntegratedHI\
        ,'Mass':CubeMass,'CubeWCS':CubeWCS,'Noise':RMS,'CubeVels':CubeVels,'PV':CubePV
        }
    return CubeInfo

def GetDistance(CatID,WallCat):
    H0=70
    CatVel=WallCat.cz[CatID]
    Distance=CatVel/H0
    return Distance

def TotalHI(cube,PixSize,BeamSize):
    cubeTot=np.nansum(cube)         #Get the total signal in cube units (Jy/beam)
    BeamArea=np.pi*BeamSize**2./(4.*np.log(2.)) #Get the beam area in degrees
    pixperbeam=BeamArea/PixSize**2. #Find the number of pixels per beam
    HISignal=cubeTot/pixperbeam #Get the integrated HI signal in Jy
    return HISignal


def MassCalc(cube,Distance,ChannelSize,BeamSize,PixSize):
    #   Get the total signal in Jy/beam
    Signal_JyBeam=np.nansum(cube)
    #       Get the beam size in degrees
    beamarea=(np.pi*BeamSize**2.)/(4.*np.log(2.))
    #   Ge the number of pixels/beam
    pixperbeam=beamarea/(abs(PixSize)*abs(PixSize))
    #   Get the signal in Jy
    Signal_Jy=Signal_JyBeam/pixperbeam
    #   Get the mass
    Mass=0.236*(Distance*1000.)**2.*Signal_Jy*abs(ChannelSize)  #Original version

    return Mass

def GetStartPixelLoc(CubeWCS,CubeHeader):
    pixcrd=[[0,0,0]]
    RealCoord=CubeWCS.wcs_pix2world(pixcrd,0)
    print("Cube Informatio n- Start Pixel", RealCoord)
    CubeHeader['CRLOC1']=RealCoord[0,0]
    CubeHeader['CRLOC2']=RealCoord[0,1]
    return CubeHeader

def GetVels(CubeHeader):
    Vels=[]
    for i in range(CubeHeader['NAXIS3']):
        VTemp=float(i-CubeHeader['CRPIX3']+1)*CubeHeader['CDELT3']+CubeHeader['CRVAL3']
        Vels.append(VTemp)
    return Vels


def CalcPVDiagram(Angle,BeamSize,CubeData,Vels,PixCenter):
    print(Angle,BeamSize)
    AngUse=(Angle+90.)*np.pi/180.
    PV=np.zeros([np.shape(CubeData)[2],np.shape(CubeData)[0]])
  
    for i in range(np.shape(CubeData)[1]):
        for j in range(np.shape(CubeData)[2]):
            X=(j-PixCenter[0])
            Y=(i-PixCenter[1])
            XP=X*np.cos(-AngUse)-Y*np.sin(-AngUse)
            YP=X*np.sin(-AngUse)+Y*np.cos(-AngUse)
            
            k=int(round(XP+PixCenter[0]))
            l=int(round(YP+PixCenter[1]))
            
            #print(i,j,k,l,XP,YP,np.nansum(CubeData[:,i,j]))
            if k >= 0 and k < np.shape(CubeData)[2]:
                #print("inloop",i,j,k,l,XP,YP,np.nansum(CubeData[:,i,j]))
                #print(i,j,k,YP)
                if abs(YP) <= BeamSize/2.:
                    # print(i,j,X,Y,XP,YP)
                    for m in range(np.shape(CubeData)[0]):
                        #print("NanCheck",np.isnan(CubeData[m,i,j]))
                        if np.isnan(CubeData[m,i,j]):
                            PV[k,m]+=0.
                        else:
                            PV[k,m]+=CubeData[m,i,j]
        #print(i,j,k,m,PV[k,m],CubeData[m,i,j],CubeData[m,j,i])
    return PV



def MakeModelPV(CubeFile,ObsDict,WallCat,CatID):
    CubeHDU=fits.open(CubeFile)
    CubeHeader=CubeHDU[0].header
    CubeData=CubeHDU[0].data
    
    BeamSize=CubeHeader['BMAJ']
    PixSize=CubeHeader['CDELT2']

    
    CubeVels=np.array(GetVels(CubeHeader))
    #    print(CubeVels)
    try:
        if CubeHeader['CUNIT3']=='km/s':
            CubeVels=CubeVels
        else:
            CubeVels=CubeVels/1000.
#       If no unit type, assume m/s and switch to km/s
    except:
        CubeVels=CubeVels/1000.
  
#    print(CubeVels)
    
    CenterPos=[[WallCat.ra[CatID],WallCat.dec[CatID],0]]
    CentPix=ObsDict['CubeWCS'].wcs_world2pix(CenterPos,0)
    CentPix=[CentPix[0,0],CentPix[0,1]]
    
    
    CubePV=CalcPVDiagram(WallCat.kin_pa[CatID],BeamSize/PixSize,CubeData,CubeVels,CentPix)
    ModelCube={'CubeHeader':CubeHeader,'CubeData':CubeData,'PV':CubePV,'CubeVels':CubeVels}
    
    CubeHDU.close()
    
    return ModelCube







