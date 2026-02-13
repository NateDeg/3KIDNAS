import numpy as np
import pandas as pd
import os.path
from os import path
import os as os

    
def GetGalaxyDictNames(step,Cat,RTDict):
    #   Initialize a galaxy dictionary object
    GalaxyDict={}
    #   Set the name of the object in the galaxy dictionary
    GalaxyDict['name']=Cat['name'][step]
    #   Get an object name with underscores
    GalaxyDict['name_underscore']=GalaxyDict['name'].replace(' ','_')
    #   Get the name of the folder where the source data lives
    GalaxyDict['cubefolder']=RTDict['SourceFolder']+GalaxyDict['name_underscore']+"/"
    #   Get the name of the frequency cube -- use a general method due to possible naming convention changes
    SourceFiles=os.listdir(GalaxyDict['cubefolder'])
    subs = 'cube.fits'
    GalaxyDict['FreqCubeName']=GalaxyDict['cubefolder']+ [i for i in SourceFiles if subs in i][0]
    subs = 'mask.fits'
    GalaxyDict['MaskName']=GalaxyDict['cubefolder']+ [i for i in SourceFiles if subs in i][0]
    #   Set the name of the processed velocity cube
    GalaxyDict['VelCubeName']=GalaxyDict['name_underscore']+"_VelCube.fits"
    #   Set the name of the galaxy fit parameter file
    GalaxyDict['FitParameterFile']=GalaxyDict['name_underscore']+"_RTParameters.py"
    
    #   Set the size of the galaxy beam in pixels
    GalaxyDict['Beam_Pix']=6.
    
    return GalaxyDict

def WriteSingleGalaxyIni(GalaxyDict,RTDict,SGDict):
    f=open(GalaxyDict['FitParameterFile'],'w')
    
    Str="CubeName='"+GalaxyDict['VelCubeName']+"'\n"
    Str+="MaskName='"+GalaxyDict['MaskName']+"'\n"
    Str+="ObjName='"+GalaxyDict['name_underscore']+"'\n"
    Str+="TargFolder='"+RTDict['TargFolder']+"'\n"
    Str+="nBootstraps="+str(RTDict['nBootstraps'])+"\n"
    Str+="nProcessors_Bootstraps="+str(RTDict['nProcessors_Bootstrap'])+"\n"
    Str+="PA_Estimate="+str(GalaxyDict['PA_Estimate'])+"\n"
    Str+="Inc_Estimate="+str(GalaxyDict['IncEst'])+"\n"
    Str+="\n"
    for key in SGDict.keys():
        print(key)
        Str+=key+"='"+SGDict[key]+"'\n"
    f.write(Str)
    f.close()
