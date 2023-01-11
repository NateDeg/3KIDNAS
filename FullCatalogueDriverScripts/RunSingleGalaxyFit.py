import numpy as np
import pandas as pd
import os.path
from os import path

    
def GetGalaxyDictNames(step,Cat,RTDict):
    #   Initialize a galaxy dictionary object
    GalaxyDict={}
    #   Set the name of the object in the galaxy dictionary
    GalaxyDict['name']=Cat['name'][step]
    #   Get an object name with underscores
    GalaxyDict['name_underscore']=GalaxyDict['name'].replace(' ','_')
    #   Get the name of the folder where the source data lives
    GalaxyDict['cubefolder']=RTDict['SourceFolder']+GalaxyDict['name_underscore']+"/"
    #   Get the name of the frequency cube
    GalaxyDict['FreqCubeName']=GalaxyDict['cubefolder']+GalaxyDict['name_underscore']+"_cube.fits"
    #   Get the name of the mask
    GalaxyDict['MaskName']=GalaxyDict['cubefolder']+GalaxyDict['name_underscore']+"_mask.fits"
    #   Set the name of the processed velocity cube
    GalaxyDict['VelCubeName']=GalaxyDict['name_underscore']+"_VelCube.fits"
    #   Set the name of the galaxy fit parameter file
    GalaxyDict['FitParameterFile']=GalaxyDict['name_underscore']+"_RTParameters.py"
    
    return GalaxyDict
    
def GetGeometryEstimates(step,Cat,GalaxyDict):
    GalaxyDict['PA_Estimate']=Cat['kin_pa'][step]
    IncEst=np.arccos(Cat['ell_min'][step]/Cat['ell_maj'][step])*180./np.pi
    GalaxyDict['IncEst']=IncEst

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

def RunSingleGalaxyFit(GalaxyDict,RTDict):
    print(RTDict.keys())
