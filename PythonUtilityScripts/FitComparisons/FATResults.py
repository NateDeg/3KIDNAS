#!/usr/bin/env python3
import numpy as np

from . import MCGGenerator as MCG
from . import CubeComparison as CC
from . import CubeInformation as CI



def GetFATFit(CatID,WallCat,Path,BaseName,ObsDict,FolderNameSwitch):
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    if FolderNameSwitch==0:
        OutputDir=Path+"/"+BaseName+"_"+NameSuffix+"/"
    elif FolderNameSwitch==1:
        OutputDir=Path+"/"+NameSuffix+"/"
    FinalParamsFile=OutputDir+"Finalmodel/FinalModel.def"
    CubeFile=OutputDir+"Finalmodel/FinalModel.fits"

    print(FinalParamsFile)
    Disk1,Disk2=ReadFATParamsFile(FinalParamsFile)
    
    Disk1=FATCentConversion(Disk1,ObsDict)
    Disk2=FATCentConversion(Disk2,ObsDict)

    if Disk1['FITAchieved']=='True':
        MCG_FATModelCube=MCG.MakeMCGModel(Disk1,ObsDict['CubeHeader'])
        ResidTot,chi2=CC.CubeCompare(ObsDict,MCG_FATModelCube)
        Disk1['RESID']=ResidTot
        Disk1['CHI2']=chi2
        ModelCube=CI.MakeModelPV(CubeFile,ObsDict,WallCat,CatID)
        Disk1['CUBE']=ModelCube
        Disk1['R_SD']=Disk1['R']

    if Disk2['FITAchieved']=='True':
        MCG_FATModelCube=MCG.MakeMCGModel(Disk2,ObsDict['CubeHeader'])
        ResidTot,chi2=CC.CubeCompare(ObsDict,MCG_FATModelCube)
        Disk2['RESID']=ResidTot
        Disk2['CHI2']=chi2
        Disk2['CUBE']=Disk1['CUBE']
        Disk2['R_SD']=Disk1['R']

    #ResidTot,chi2=BBaroloCubeCompare(ObsDict,MCG_BBaroloModelCube)

    return Disk1,Disk2


def ReadFATParamsFile(FileName):
    try:
        f=open(FileName,"r")
    except:
        print("No Final FAT Model")
        Disk1=BadFATFit()
        Disk2=BadFATFit()
        return Disk1,Disk2


    FullFile=f.readlines()  #   Read file into lines
    #print(np.shape(FullFile))
    #       Line 16 is the number of radial bins
    #print(FullFile[16])
    
    SpaceChar=(" ","   ","  ")
    #    nR=ExtractLineValues(FullFile[16],SpaceChar[0])
    #nR=int(nR)

    Disk1={'FITAchieved':'True'}
    for i in range(17,31,1):
        Disk1=ExtractLineValues(FullFile[i],SpaceChar,Disk1)
        Disk1['VRAD']=Disk1['R']*0.

    Disk2={'FITAchieved':'True'}
    for i in range(32,48,1):
        Disk2=ExtractLineValues(FullFile[i],SpaceChar,Disk2)
        Disk2['VRAD']=Disk1['R']*0.

    Disk2['R']=Disk1['R']
# print(Disk1['SURFDENS'])
#print(Disk2['SURFDENS'])

    f.close()
    return Disk1,Disk2

def ExtractLineValues(LineStr,Spacer,FATDict):
    #print(LineStr)
#       Split by equals line
    EqualSplit=LineStr.split("=")
    #print("EqualSplit", EqualSplit[0])
    
    #       Get the keyword and the line formating
    ReadFormat,VarName=ReadFATKeyWord(EqualSplit[0])
    if ReadFormat == -1:
        #       If it's not one of the recognized keywords, don't do anything
        return FATDict
    else:
        #       Remove the \n character from the line
        values=EqualSplit[1].split("\n")[0]
        #       Split into an array according to the spacer and format for the variable
        #        Arr=values.split(Spacer[ReadFormat])
        Arr=values.split()

#print(VarName,ReadFormat,Arr)
    ArrFloat=np.array(list(map(float,Arr)))
    FATDict[VarName]=ArrFloat
    return FATDict

def FATCentConversion(Disk,ObsDict):
    w=ObsDict['CubeWCS']
    if Disk['FITAchieved']=='True':
        X=[]
        Y=[]
        for i in range(np.shape(Disk['RA'])[0]):
            RealCoords=[[Disk['RA'][i],Disk['DEC'][i],0]]
            PixCoords=w.wcs_world2pix(RealCoords,0)
            X.append(PixCoords[0,0])
            Y.append(PixCoords[0,1])
        Disk['XCENTER']=X
        Disk['YCENTER']=Y
    else:
        Disk['XCENTER']=[]
        Disk['YCENTER']=[]
    return Disk


def ReadFATKeyWord(KeyStr):
    TestStr=KeyStr.strip('#')
    TestStr=TestStr.strip()
    VarName=[]

    ReadFormat=-1
    if TestStr == 'NUR':
        ReadFormat=0
        VarName='nR'
    elif TestStr == 'RADI':
        ReadFormat=0
        VarName='R'
    elif TestStr == 'VROT':
        ReadFormat=0
        VarName='VROT'
    elif TestStr == 'VROT_ERR':
        ReadFormat=2
        VarName='VROT_ERR'
    elif TestStr == 'Z0':
        ReadFormat=0
        VarName='Z0'
    elif TestStr == 'SBR':
        ReadFormat=0
        VarName='SURFDENS'
    elif TestStr == 'INCL':
        ReadFormat=0
        VarName='INCLINATION'
    elif TestStr == 'INCL_ERR':
        ReadFormat=2
        VarName='INC_ERR'
    elif TestStr == 'PA':
        ReadFormat=0
        VarName='POSITIONANGLE'
    elif TestStr == 'PA_ERR':
        ReadFormat=2
        VarName='PA_ERR'
    elif TestStr == 'XPOS':
        ReadFormat=0
        VarName='RA'
    elif TestStr == 'YPOS':
        ReadFormat=0
        VarName='DEC'
    elif TestStr == 'VSYS':
        ReadFormat=0
        VarName='VSYS'
    elif TestStr == 'SDIS':
        ReadFormat=0
        VarName='VDISPERSION'
    elif TestStr == 'SDIS_ERR':
        ReadFormat=2
        VarName='VDISP_ERR'


    elif TestStr == 'VROT_2':
        ReadFormat=0
        VarName='VROT'
    elif TestStr == 'VROT_ERR':
        ReadFormat=2
        VarName='VROT_2_ERR'
    elif TestStr == 'Z0_2':
        ReadFormat=0
        VarName='Z0'
    elif TestStr == 'SBR_2':
        ReadFormat=0
        VarName='SURFDENS'
    elif TestStr == 'INCL_2':
        ReadFormat=0
        VarName='INCLINATION'
    elif TestStr == 'INCL_2_ERR':
        ReadFormat=2
        VarName='INC_ERR'
    elif TestStr == 'PA_2':
        ReadFormat=0
        VarName='POSITIONANGLE'
    elif TestStr == 'PA_2_ERR':
        ReadFormat=2
        VarName='PA_ERR'
    elif TestStr == 'XPOS_2':
        ReadFormat=0
        VarName='RA'
    elif TestStr == 'YPOS_2':
        ReadFormat=0
        VarName='DEC'
    elif TestStr == 'VSYS_2':
        ReadFormat=0
        VarName='VSYS'
    elif TestStr == 'SDIS_2':
        ReadFormat=0
        VarName='VDISPERSION'
    elif TestStr == 'SDIS_2_ERR':
        ReadFormat=2
        VarName='VDISP_ERR'

    return ReadFormat,VarName

def BadFATFit():
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':'False'}
    return Disk

