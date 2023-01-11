import sys as sys
import os as os
import copy as copy

def GetRuntimeArguments():
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
    ParamFile=Commands[1]
    return ParamFile
        
    
def ReadParamFile(FileName):
    print("About to import parameter file ", FileName)
    print("Absolute path to parameter file:",os.path.abspath(FileName))
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
    #import importlib as il
    import importlib.util
    print("Package", Package)
    print("ModName", ModName)
    spec = importlib.util.spec_from_file_location(name="RTParams",location=AbsPath)
    #ModTest=il.import_module(ModName,package=Package)
    ModTest=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ModTest)
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


