#!/usr/bin/env python3
import numpy as np
import os

def MakeMCGModel(ParamDict,DataHeader,):
    MCG_DatacubeHeaderName="TempMCG_DCHeader.in"
    MCG_TiltedRingFileName="TempMCG_TR.in"
    MCG_MainInFile="TempMCG.in"
    MCG_ModelName="TempModel"
    
    FinalModelCube=MCG_ModelName+"/"+MCG_ModelName+"_ConvolvedSourceCube.fits"
    
    MCGNames={'OutputName':MCG_ModelName,'MainFile':MCG_MainInFile\
        ,'DatacubeFile':MCG_DatacubeHeaderName,'TiltedRingFile':MCG_TiltedRingFileName}
    
    
    MakeMCG_MainFile(MCGNames)
    
    MakeMCGDataCubeHeader(MCG_DatacubeHeaderName,DataHeader)

    MakeMCG_TiltedRingInput(MCG_TiltedRingFileName,ParamDict)

    LocMCGPath="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT14/MCGSuite/Programs/MockCubeGenerator"

    MakeCubeCmd=LocMCGPath+" " + MCG_MainInFile
    os.system(MakeCubeCmd)
    CleanCmd="rm " + MCG_MainInFile +" " + MCG_DatacubeHeaderName + " " + MCG_TiltedRingFileName
#os.system(CleanCmd)

    return FinalModelCube


def MakeMCG_MainFile(MCGNames):
    file=open(MCGNames['MainFile'],"w")
    
    file.write("#    Base name for the output folder and all file names \n")
    file.write(MCGNames['OutputName']+"\n")
 
    file.write("#    Minimal, Moderate, of Full Outputs (0,1, and 2 respectively)\n")
    file.write("2 \n")
     
    file.write("#    Name of the file containing the cube and beam definitions\n")
    file.write(MCGNames['DatacubeFile']+"\n")
     
    file.write("#    Name of the file containing the underlying ring model\n")
    file.write(MCGNames['TiltedRingFile']+"\n")
     
    file.write("#    Noise Unit Switch (0=mJy/beam)\n")
    file.write("0 \n")
     
    file.write("#    Noise Value\n")
    file.write("1.6 \n")
     
    file.write("#    Random Seed\n")
    file.write("1 \n")

    file.close()

def MakeMCGDataCubeHeader(MCGFileName,CubeHeader):

    file=open(MCGFileName,"w")
    file.write("#    number of pixels and channel\n")
    file.write(str(CubeHeader['NAXIS1'])+"\t" +str(CubeHeader['NAXIS2']) \
               +"\t" +str(CubeHeader['NAXIS3']) +"\n")
    file.write("#    Pixel Size Units (0=degrees, 1=arcsec)  - Channel Size Units (0=m/s, 1=km/s) - Reference Pixel Units (0=degrees,1=arcsec) - Reference Channel Units (0=m/s, 1=km/s)  — Beam axis Units (0=arcsec) — Beam Rotation Angle Units (0=degrees)\n")
    file.write("0\t 0 \t 0 \t 0 \t 1 \t 0\n")

    file.write("#    Pixel dimensions and channel size\n")
    file.write(str(CubeHeader['CDELT1']) +"\t" +str(CubeHeader['CDELT2'])\
               +"\t"+str(CubeHeader['CDELT3']) +"\t"+" \n")

    file.write("#    Reference Pixel(channel) in each dimension (CRPIX values)\n")
    file.write("1" +"\t" +"1"\
               +"\t"+str(CubeHeader['CRPIX3']) +"\t"+" \n")
    
    file.write("#    Reference value in each dimension (CRVAL values)\n")
    file.write(str(CubeHeader['CRLOC1']) +"\t" +str(CubeHeader['CRLOC2'])\
               +"\t"+str(CubeHeader['CRVAL3']) +"\t"+" \n")

    file.write("#    Beam dimensions (major, minor, position angle)\n")
    file.write(str(CubeHeader['BMAJ']*3600.) +"\t" +str(CubeHeader['BMIN']*3600.)
               +"\t"+str(CubeHeader['BPA']) +"\t"+" \n")
    
    
    file.write("#    number of sigma lengths to reach with the beam\n")
    file.write("5. \n")
    
    file.write("#    Type of velocity smoothing to use (0=none, 1=Gaussian)\n")
    file.write("0 \n")
    
    file.write("#    The velocity smoothing sigma if using Gaussian smoothing (km/s)\n")
    file.write("10. \n")

    file.close()

def MakeMCG_TiltedRingInput(MCGFileName,ParamDict):
    
    nRings=np.shape(ParamDict['R'])[0]
    Rwidth=ParamDict['R'][1]-ParamDict['R'][0]
    file=open(MCGFileName,"w")
    
    
    file.write("#        The format of the input file — 1=row, 2= column\n")
    file.write("2 \n")
    
    file.write("#   Number of Rings in the model\n")
    file.write(str(nRings)+"\n")


    file.write("#    The cloud mode you are using\n")
    file.write("0 \n")

    file.write("#    The base cloud surface density\n")
    file.write("10. \n")
    
    file.write("#    Central Position Switch (0=degrees, 1=arcsec), Inc/PA Unit switch (0=degrees, 1=arcsec), Velocity Units (0=m/s, 1=km/s), Brightness Units (0=Jy km/s arcsec^-2)\n")
    file.write("0 \t 0 \t 1 \t  0 \n")
    
    file.write("#    The parameters in each radial bin\n")
    file.write("#    Rmid    Rwidth    Xcent    Ycent    Inc    PA    VSys    VRot    VRad    Vvert    VDisp    dvdz    Sigma        z0    zGradStart\n")
    
    #    print(len(ParamDict['R']), len(ParamDict['SURFDENS']))
    for i in range(nRings):
        #print(i)
        file.write(str(ParamDict['R'][i])+"\t"+str(Rwidth) +"\t" + str(ParamDict['RA'][i]) \
                   + "\t"+ str(ParamDict['DEC'][i]) +"\t" + str(ParamDict['INCLINATION'][i])\
                   + "\t" +str(-ParamDict['POSITIONANGLE'][i]+180) +"\t"\
                   + str(ParamDict['VSYS'][i])+ "\t" + str(ParamDict['VROT'][i]) \
                   +"\t" + str(ParamDict['VRAD'][i]) + "\t" + "0." +"\t" \
                   + str(ParamDict['VDISPERSION'][i]) + "\t"+ "0." + "\t" \
                   + str(ParamDict['SURFDENS'][i])+"\t" + str(ParamDict['Z0'][i])\
                   + "\t" + str(5.*ParamDict['Z0'][i]) +   "\n")

    file.close()
