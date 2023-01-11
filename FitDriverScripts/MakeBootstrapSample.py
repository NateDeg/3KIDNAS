import sys as sys
import os as os
import copy as copy
import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs


    
def SetBootstrapVariableValues(GalaxyDict,Step):
    GalaxyDict['ObjNameU']=GalaxyDict['ObjName']+"_Bootstrap_"+str(Step)
    
    return GalaxyDict
    
    
def MakeBootstrapSample(GeneralDict,GalaxyDict,Step):
    
    #   Start by writing the bootstrap input file
    GalaxyDict['BSName'],GalaxyDict['BSFile']=WriteBootstrapFile(GeneralDict,GalaxyDict,Step)
    #   And make the output folder for the bootstrapped cube
    os.makedirs(GalaxyDict['BootstrapFolder'], exist_ok=True)
    #   Now we can run the bootstrap sample
    BSCmd=GeneralDict['BootstrapExecPath']+" "+GalaxyDict['BSFile']
    os.system(BSCmd)
    #   If everything had gone well, the cube to be used for the SoFiA and WRKP Bootstrap analysis is the output bootstrap cube
    #   Before moving the bootstrap cube, make sure the header has integers for the NAXIS key words
    AdjustBootstrapHeader(GalaxyDict)
    #       Start by putting the bootstrap cube in the correct folder
    MVCmd="mv "+GalaxyDict['BSName']+".fits "+GalaxyDict['BootstrapFolder']+"."
    os.system(MVCmd)
    GalaxyDict['CubeNameU']=GalaxyDict['BootstrapFolder']+GalaxyDict['BSName']+".fits"
    #   Clean up by getting rid of the bootstrap file
    ClnCmd="rm "+GalaxyDict['BSFile']
    os.system(ClnCmd)

    return GalaxyDict
    
def WriteBootstrapFile(GeneralDict,GalaxyDict,Step):

    
    BootstrapName=GalaxyDict['ObjNameU']
    BootstrapFileName=GalaxyDict['TargFolder']+BootstrapName+".txt"
   
    
   
    f=open(BootstrapFileName,'w')
    ExplanatoryStr="#    Name of the input data cube\n"
    f.write(ExplanatoryStr)
    f.write(GalaxyDict['CubeName']+"\n")
    
    ExplanatoryStr="#    Name of the model data cube\n"
    f.write(ExplanatoryStr)
    f.write(GalaxyDict['BestFitModel']['ModelCube']+"\n")
    
    ExplanatoryStr="#    The base name for the bootstrap cube\n"
    f.write(ExplanatoryStr)
    f.write(BootstrapName+"\n")
    
    ExplanatoryStr="#    Size of the resampling blocks in beams and channels\n"
    f.write(ExplanatoryStr)
    f.write("1\n")
    
    
    Model=GalaxyDict['BestFitModel']
    GeoStr=str(Model['XCENTER'][0])+"\t"
    GeoStr+=str(Model['YCENTER'][0])+"\t"
    
    
    RefVel=GalaxyDict['CubeHeader']['CRVAL3']
    RefChan=GalaxyDict['CubeHeader']['CRPIX3']
    dV=GalaxyDict['CubeHeader']['CDELT3']/1000.
    
    DeltaV=Model['VSYS'][0]-RefVel/1000.
    #print(Model['VSYS'][0],RefVel/1000.,DeltaV,dV)
    VCenter=DeltaV/dV+RefChan
    #print(VCenter)
    
    GeoStr+=str(VCenter)+"\t"
    GeoStr+=str((Model['POSITIONANGLE'][0]+90.)*np.pi/180.)+"\t"
    GeoStr+=str(Model['INCLINATION'][0]*np.pi/180.)+"\n"
    
    ExplanatoryStr="#    The geometry used for bootstrap resampling (centre pt. in pixels plus PA & INC in radians) \n"
    f.write(ExplanatoryStr)
    f.write(GeoStr)
    
    
    f.close()
 
    return BootstrapName,BootstrapFileName

def AdjustBootstrapHeader(GalaxyDict):
    BSCubeName=GalaxyDict['BSName']+".fits"
    Cube=fits.open(BSCubeName)
    
    keys=['NAXIS1','NAXIS2','NAXIS3']
    for key in keys:
        OriKey=Cube[0].header[key]
        NewKey=int(OriKey)
        print(type(NewKey))
    
