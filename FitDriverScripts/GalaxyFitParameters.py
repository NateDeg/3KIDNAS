"""
GalaxyFitParameters.py
=======

**Author:** Nathan Deg

**Description:**
This module contains functions that both define and collect key parameters for running 3KIDNAS.
"""

import sys as sys
import os as os
import copy as copy

def GetRuntimeArguments():
    """
    This function gets the name of the python file containing the runtime arguments.
    Returns
        -------
        ParamFile : string
            The name of the python formatted file containing the required runtime arguments
    """

    #   Grab the command line arguments
    Commands=sys.argv
    #   Make sure the user has supplied a python parameter file
    if len(Commands)==1:
        print("You must supply a python parameter file with runtime arguments")
        exit()
    #   Check that the python parameter file exists
    FileCheck=os.path.isfile(Commands[1])
    if FileCheck == False:
        print("The supplied parameter file does not exist")
        print(Commands[1])
        exit()
    #   Finally store the name of the parameter file and return it.
    ParamFile=Commands[1]
    return ParamFile
        
    
def ReadParamFile(FileName):
    """
    This function reads the required python runtime input file.  Because this is a python file (which gives some increased flexibility, it will be read in as a module.
    
    Parameters
        ----------
        FileName : string
            The name of the file to be read.
        
    Returns
        -------  
        ModTest : python module
            A module containing all the runtime arguments
    """
    #print("About to import parameter file ", FileName)
    #print("Absolute path to parameter file:",os.path.abspath(FileName))
    #   Convert the path to an absolute path
    AbsPath=os.path.abspath(FileName)
    #   Figure out the path and the module names by splitting the absolute path
    Package=AbsPath.rsplit("/",1)[0]
    ModName=AbsPath.rsplit("/",1)[1]
    #   Remove the .py from the ModName
    ModName=ModName.rsplit(".",1)[0]
    #   Now import Sys
    import sys
    #   And add the path to the parameter file to sys.path
    sys.path.append(Package)
    #   Finally use the importlib module to import the module
    import importlib.util
    #print("Package", Package)
    #print("ModName", ModName)
    #   Start by setting the instance for loading in a module
    spec = importlib.util.spec_from_file_location(name="RTParams",location=AbsPath)
    #   And now set it so that the module can be loaded.
    ModTest=importlib.util.module_from_spec(spec)
    #   Finally load the module
    spec.loader.exec_module(ModTest)
    #   Return the module contents
    return ModTest

def OverwriteDefaults(GeneralDict,RTParams,KeyRTParams):
    #   First loop through all the default key values and see if they are re-set in the runtime parameters
    #print(GeneralDict.keys())
    for key in GeneralDict.keys():
        try:
            RTVal=vars(RTParams)[key]
            print("Changing value for parameter", key, "to", RTVal)
            GeneralDict[key]=RTVal
        except:
            pass
    #   We also want to check the fitting options
    for key in GeneralDict['FittingOptsDict']['MethodParams']:
        try:
            RTVal=vars(RTParams)[key]
            print("Changing value for parameter", key, "to", RTVal)
            GeneralDict['FittingOptsDict'][key]=RTVal
        except:
            pass
            
    #       The fitting options need to be split into constant and fixed keys
    for key in GeneralDict['FittingOptsDict']['FittingParams']:
        try:
            RTVal=vars(RTParams)[key]
            RTVal=RTVal.split()
            print("Changing value for parameter", key, "to", RTVal)
            GeneralDict['FittingOptsDict'][key]['Constant']=RTVal[0]
            GeneralDict['FittingOptsDict'][key]['Fixed']=RTVal[1]
            
        except:
            pass
    
    #   Now set all the key parameters from the runtime file and story them in a galaxy specific dictionary.  This will make it easier to generalize later functions to take in the general dictionary and a specific galaxy dictionary.
    GalaxyDict={}
    KeyList=KeyRTParams[0][:]
    for key in KeyRTParams:
        keyU=key[0]
        try:
            RTVal=vars(RTParams)[keyU]
            GalaxyDict[keyU]=RTVal
        except:
            print("Input file is missing value for: ", keyU)
            print("Stopping run here")
            exit()
    return GeneralDict,GalaxyDict
    
def CheckParamTypes(GalaxyDict,KeyRTParams):
    print("Checking to make sure all variables have a correct type/format")
    for key in KeyRTParams:
        keyU=key[0]
        TargType=key[1]
        if type(GalaxyDict[keyU]) != TargType:
            print("Wrong variable type for ", keyU)
            print("Current entry is ", GalaxyDict[keyU])
            print("And has type ", type(GalaxyDict[keyU]))
            print("Stopping run here")
            exit()


def SetCurrentRunVariables(GalaxyDict):
    WorkingVariableList=["CubeNameU","MaskNameU","ObjNameU","Inc_EstimateU","PA_EstimateU","TargFolderU"]
    OriginalVariableList=["CubeName","MaskName","ObjName","Inc_Estimate","PA_Estimate","TargFolder"]
    for i in range(len(WorkingVariableList)):
        x=WorkingVariableList[i]
        y=OriginalVariableList[i]
        GalaxyDict[x]=copy.deepcopy(GalaxyDict[y])
    return GalaxyDict


def DefaultRuntimeOptions():
    """
    This function sets a variety of default runtime options
    Returns
        -------
        DefaultOpts : dictionary
            A dictionary containing the default runtime parameter values
    """

    DefaultOpts={}
    DefaultOpts['FittingAlgorithm']=1
    DefaultOpts['LikeFnc']=1
    DefaultOpts['ParamConversion']=1
    DefaultOpts['MomMapCalc']=3
    DefaultOpts['CentMethod']=0
    DefaultOpts['ShapeMethod']=0
    DefaultOpts['SizeMethod']=0
    DefaultOpts['VProfileMethod']=0
    DefaultOpts['SD_LogSwitch']=0
    DefaultOpts['CloudMode']=0
    DefaultOpts['cdens']=100
    DefaultOpts['BeamSmear']=2.5
    DefaultOpts['SDLim_NoiseFactor']=1.0
    DefaultOpts['NRings']=-1
    DefaultOpts['RingPerBeam']=2
    DefaultOpts['XFitting']={'Constant':"T",'Fixed':"F"}
    DefaultOpts['YFitting']={'Constant':"T",'Fixed':"F"}
    DefaultOpts['IncFitting']={'Constant':"T",'Fixed':"F"}
    DefaultOpts['PAFitting']={'Constant':"T",'Fixed':"F"}
    DefaultOpts['VSysFitting']={'Constant':"T",'Fixed':"F"}
    DefaultOpts['VRotFitting']={'Constant':"F",'Fixed':"F"}
    DefaultOpts['VRadFitting']={'Constant':"T",'Fixed':"T"}
    DefaultOpts['VDispFitting']={'Constant':"T",'Fixed':"T"}
    DefaultOpts['VVertFitting']={'Constant':"T",'Fixed':"T"}
    DefaultOpts['dvdzFitting']={'Constant':"T",'Fixed':"T"}
    DefaultOpts['SDFitting']={'Constant':"F",'Fixed':"F"}
    DefaultOpts['ZHeightFitting']={'Constant':"T",'Fixed':"T"}
    DefaultOpts['ZGradFitting']={'Constant':"T",'Fixed':"T"}

    DefaultOpts['LineKeyMap']={'FittingAlgorithm':1,'LikeFnc':3,'ParamConversion':5,'MomMapCalc':7,'CentMethod':12,'ShapeMethod':14,'SizeMethod':16,'VProfileMethod':18,'SD_LogSwitch':20,'CloudMode':22,'cdens':24,'BeamSmear':26,'SDLim_NoiseFactor':28,'NRings':30,'RingPerBeam':32,'XFitting':34,'YFitting':36,'IncFitting':38,'PAFitting':40,'VSysFitting':42,'VRotFitting':44,'VRadFitting':46,'VDispFitting':48,'VVertFitting':50,'dvdzFitting':52,'SDFitting':54,'ZHeightFitting':56,'ZGradFitting':58}
    
    DefaultOpts['FittingParams']=['XFitting','YFitting','IncFitting','PAFitting','VSysFitting','VRotFitting','VRadFitting','VDispFitting','VVertFitting','dvdzFitting','SDFitting','ZHeightFitting','ZGradFitting']
    DefaultOpts['MethodParams']=['FittingAlgorithm','LikeFnc','ParamConversion','MomMapCalc','CentMethod','ShapeMethod','SizeMethod','VProfileMethod','SD_LogSwitch','CloudMode','cdens','BeamSmear','SDLim_NoiseFactor','NRings','RingPerBeam']

    return DefaultOpts
