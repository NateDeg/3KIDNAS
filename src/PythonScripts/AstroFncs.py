import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

def FuncDict():
    AstronFncDict={'RedShiftConv':RedShiftConv,'DistEst':DistEst,'MassFromFlux':MassFromFlux,'CalcRHI':CalcRHI,'ProjectedSize_ArcSecond':ProjectedSize_ArcSecond,'CalcRA_Dec_FromCube':CalcRA_Dec_FromCube_And_Center,'CalcCenter_FromCube':CalcCenter_FromCube_And_RADEC}
    return AstronFncDict

def RedShiftConv(Freq,RestFreq):
    z=(RestFreq-Freq)/Freq
    Vel=z*2.9979245e8
    return Vel
    
def DistEst(H0,Vel):
    Distance=Vel/H0
    return Distance
    
def MassFromFlux(Distance,Flux,dF):
    #   This assumes that the flux is in units of Jy freq
    #   Need dF in velocity terms
    RestFreq=1.42040575179E+09
    FreqT=RestFreq+dF
    dV=RedShiftConv(FreqT.values,RestFreq)/1000.
    print(dV)
    
    
    Mass=0.236*(Distance*1000.)**2.*Flux/abs(dF)
    MassT=Mass.astype(float)*np.abs(dV)
    Mass=np.log10(Mass.astype(float))
    print(Mass,np.log10(MassT))
    return Mass
   
def CalcRHI(logMHI):
#Use the HI mass - diameter relation from
#Wang+2016 to estimate RHI from logMHI.
#-->INPUT: log(MHI/Mo)
#-->OUTPUT: RHI in kpc
#--------------------
    slope = 0.506
    intercept = -3.293
    logDHI = slope*logMHI + intercept
    RHI = (10**logDHI)/2.
    return RHI
   
   
def ProjectedSize_ArcSecond(Size,Distance):
    #   Note that this assumes the size is in kpc and the distance is in Mpc
    Size_P=Size/(Distance*1000) * 206265
    return Size_P



def CalcRA_Dec_FromCube_And_Center(X,Y,CubeDict):
    w=CubeDict['CubeWCS']
    RA=[]
    DEC=[]
    for i in range(np.shape(X)[0]):
        pixcrd=[[X[i], Y[i],0]]
        CentCoord=w.wcs_pix2world(pixcrd,0)
        RA.append(CentCoord[0,0])
        DEC.append(CentCoord[0,1])

    return RA,DEC

def CalcCenter_FromCube_And_RADEC(RA,DEC,CubeDict):
    w=CubeDict['CubeWCS']
    X=[]
    Y=[]
    for i in range(len(RA)):
        RealCoords=[[RA[i],DEC[i],0]]
        PixCoords=w.wcs_world2pix(RealCoords,0)
        X.append(PixCoords[0,0])
        Y.append(PixCoords[0,1])

    return X,Y
