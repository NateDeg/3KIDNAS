import numpy as np
from decimal import Decimal

import os.path
from os import path
import re

def LoadFATModelFile(FileName):
    Disk1,Disk2=ReadFATParamsFile(FileName)
    
    
    if Disk1['FITAchieved']:
        Disk1['R_SD']=Disk1['R']
        Disk1=FATErrAdjust(Disk1)
        Disk1['FitForAveraging']=True
        Disk1['SURFDENS_FACEON']=Disk1['SURFDENS']  #The surface density is already face-on, but needs to be set for MCG later on
        
    if Disk2['FITAchieved']:
        Disk2['R_SD']=Disk2['R']
        Disk2=FATErrAdjust(Disk2)
        Disk2['FitForAveraging']=False
        Disk2['SURFDENS_FACEON']=Disk2['SURFDENS']
    
    Disk1['SURFDENS_FACEON2']=Disk2['SURFDENS_FACEON']
    return Disk1
        
        
def ReadFATParamsFile(FileName):
    try:
        f=open(FileName,"r")
    except:
        print("No Final FAT Model")
        Disk1=BadFATFit()
        Disk2=BadFATFit()
        return Disk1,Disk2


    FullFile=f.readlines()  #   Read file into lines

    
    SpaceChar=(" ","   ","  ")

    Disk1={'FITAchieved':True}
    Disk2={'FITAchieved':True}
    for i in range(17,55,1):
        Disk1,Disk2=ExtractLineValues(FullFile[i],SpaceChar,Disk1,Disk2)
    #Disk1['VRAD']=Disk1['R']*0. # A radial velocity value is needed for MCG


    #for i in range(32,48,1):
    #    Disk2=ExtractLineValues(FullFile[i],SpaceChar,Disk2)
    #Disk2['VRAD']=Disk1['R']*0.

    Disk1['VRAD']=Disk1['R']*0. # A radial velocity value is needed for MCG
    Disk2['VRAD']=Disk1['R']*0.
    Disk2['R']=Disk1['R']

    f.close()
    return Disk1,Disk2
    
    
def FATErrAdjust(Disk):
    Disk['VROT_ERR']=[Disk['VROT_ERR'],Disk['VROT_ERR']]
    Disk['POSITIONANGLE_ERR']=[Disk['PA_ERR'],Disk['PA_ERR']]
    Disk['INCLINATION_ERR']=[Disk['INC_ERR'],Disk['INC_ERR']]
    Disk['SURFDENS_ERR']=[Disk['R']*0.,Disk['R']*0.]
    return Disk

def ExtractLineValues(LineStr,Spacer,FATDict1,FATDict2):
#       Split by equals line
    EqualSplit=LineStr.split("=")
    #print("Line value", EqualSplit)
    #       Get the keyword and the line formating
    ReadFormat,VarName,DiskSwitch=ReadFATKeyWord(EqualSplit[0])
    if ReadFormat == -1:
        #       If it's not one of the recognized keywords, don't do anything
        return FATDict1,FATDict2
    else:
        #       Remove the \n character from the line
        values=EqualSplit[1].split("\n")[0]
        #       Split into an array according to the spacer and format for the variable
        Arr=values.split()
    ArrFloat=np.array(list(map(float,Arr)))
    if DiskSwitch==0:
        FATDict1[VarName]=ArrFloat
    else:
        FATDict2[VarName]=ArrFloat
    
    #FATDict[VarName]=ArrFloat
    return FATDict1,FATDict2


def ReadFATKeyWord(KeyStr):
    TestStr=KeyStr.strip('#')
    TestStr=TestStr.strip()
    VarName=[]
    DiskSwitch=0

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
        DiskSwitch=1
    elif TestStr == 'VROT_2_ERR':
        ReadFormat=2
        VarName='VROT_ERR'
        DiskSwitch=1
    elif TestStr == 'Z0_2':
        ReadFormat=0
        VarName='Z0'
        DiskSwitch=1
    elif TestStr == 'SBR_2':
        ReadFormat=0
        VarName='SURFDENS'
        DiskSwitch=1
    elif TestStr == 'INCL_2':
        ReadFormat=0
        VarName='INCLINATION'
        DiskSwitch=1
    elif TestStr == 'INCL_2_ERR':
        ReadFormat=2
        VarName='INC_ERR'
        DiskSwitch=1
    elif TestStr == 'PA_2':
        ReadFormat=0
        VarName='POSITIONANGLE'
        DiskSwitch=1
    elif TestStr == 'PA_2_ERR':
        ReadFormat=2
        VarName='PA_ERR'
        DiskSwitch=1
    elif TestStr == 'XPOS_2':
        ReadFormat=0
        VarName='RA'
        DiskSwitch=1
    elif TestStr == 'YPOS_2':
        ReadFormat=0
        VarName='DEC'
        DiskSwitch=1
    elif TestStr == 'VSYS_2':
        ReadFormat=0
        VarName='VSYS'
        DiskSwitch=1
    elif TestStr == 'SDIS_2':
        ReadFormat=0
        VarName='VDISPERSION'
        DiskSwitch=1
    elif TestStr == 'SDIS_2_ERR':
        ReadFormat=2
        VarName='VDISP_ERR'
        DiskSwitch=1

    return ReadFormat,VarName,DiskSwitch



def BadFATFit():
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':False,'FitForAveraging':False}
    return Disk
