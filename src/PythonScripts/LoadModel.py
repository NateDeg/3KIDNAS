import numpy as np
from decimal import Decimal

import os.path
from os import path
import re

def LoadBestFitModelFile(FileName):
    print("Loading Avg Model")
    if path.isfile(FileName) == False:
        FitAchieved='False'
        AvgDict={'FITAchieved':FitAchieved}
        return AvgDict
    else:
        FitAchieved='True'
        #   Read the file
        Fit=open(FileName,"r")
        Lines=Fit.readlines()
        #   Figure out the number of profile values
        ProfileStartLine=18
        nR=len(Lines)-ProfileStartLine
        #   Load in the geometric parameters
        XCent,XErr=GeoLineAssign(Lines[7],nR)
        YCent,YErr=GeoLineAssign(Lines[8],nR)
        RA,RAErr=GeoLineAssign(Lines[9],nR)
        DEC,DECErr=GeoLineAssign(Lines[10],nR)
        Inc,IncErr=GeoLineAssign(Lines[11],nR)
        PA,PAErr=GeoLineAssign(Lines[12],nR)
        VSys,VSysErr=GeoLineAssign(Lines[13],nR)
        #print("Geo load test", XCent[0],YCent[0],RA[0],DEC[0],Inc[0],PA[0],VSys[0])

        #   Initialize the radial profiles
        R=np.zeros(nR)
        VRot=np.zeros(nR)
        VErr=np.zeros(nR)
        VErrInc=np.zeros(nR)
        SD=np.zeros(nR)
        SDErr=np.zeros(nR)
        #   Load in the radial profile values
        for i in range(nR):
            j=i+ProfileStartLine
            R[i],VRot[i],VErr[i],VErrInc[i],SD[i],SDErr[i]=ProfileLineAssign(Lines[j])
            #print(R[i],VRot[i])
        #   Close the file
        Fit.close()
        #   Assign values to a standard TR dictionary
        AvgDict={'R':R,'R_SD':R,'XCENTER':XCent, 'XCENTER_ERR':XErr\
            ,'YCENTER':YCent,'YCENTER_ERR':YErr\
            ,'INCLINATION':Inc,'INCLINATION_ERR':IncErr\
            ,'POSITIONANGLE':PA,'POSITIONANGLE_ERR':PAErr\
            ,'VSYS':VSys,'VSYS_ERR':VSysErr\
            ,'VROT':VRot,'VROT_ERR':VErr\
            ,'SURFDENS':SD,'SURFDENS_ERR':SDErr\
            ,'RA':RA,'RA_ERR':RAErr
            ,'DEC':DEC,'DEC_ERR':DECErr
            ,'FITAchieved':FitAchieved,'CHI2':-1}
        return AvgDict

def GeoLineAssign(Line,nR):
    S=Line.split()
    Val=np.full(nR,float(S[2].strip()))
    Err=np.full(nR,float(S[3].strip()))
    return Val,Err
    
def ProfileLineAssign(Line):
    #S=re.split('\t\t\t|\t\t|\t',Line)
    S=Line.split()
    R=float(S[0].strip())
    V=float(S[1].strip())
    VErr=float(S[2].strip())
    VErrInc=float(S[3].strip())
    SD=float(S[4].strip())
    SDErr=float(S[5].strip())
    return R,V,VErr,VErrInc,SD,SDErr
