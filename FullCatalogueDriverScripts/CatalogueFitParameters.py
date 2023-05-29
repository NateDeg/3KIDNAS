import sys as sys
import os as os
import copy as copy
import numpy as np

def SetKeyParams():
    #   Set the list of runtime parameters that must be set in the RT Params  file
    #       It is necessary to list each variable with it's type i.e. [Variable Name, Type] as the types will be checked in GalaxyFitParamters.CheckParamTypes
    KeyRTParamsList=np.array([['CatName',str]
    ,['SourceFolder',str]
    ,['nBootstraps',int]
    ,['TargFolder',str]
    ,['nProcessors',int]
    ,['nProcessors_Bootstrap',int]
    ,['KinTR',str]
    ])
    return KeyRTParamsList
    
    
def CatalogueRTImport():
    #   Get a list of required keywords in the catalogue RT parameter file
    KeyRTParams=SetKeyParams()
    #   Get the name of the runtime parameter file
    RTParamFile=GetRuntimeArguments()
    #   Parse the runtime parameter file
    RTParams=ReadParamFile(RTParamFile)
    #   Now put the variables into either that catalogue runtime dictionary or a dictionary that will be used and parsed on each individual galaxy fit
    RTDict,SGDict=SetupRTDictionaries(RTParams,KeyRTParams)
    #   Return the two dictionaries
    return RTDict,SGDict
    
    

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
    ModName=ModName.split(".")[0]
    #   Now import Sys
    import sys
    #   And add the path to the parameter file to sys.path
    sys.path.append(Package)
    #   Finally use the importlib module to import the module
    import importlib as il
    ModTest=il.import_module(ModName,package=Package)

    return ModTest

def SetupRTDictionaries(RTParams,KeyRTParams):
    #   Initialize the catalogue dictionary and the galaxy dictionary
    CatDict={}
    NewGalaxyRTVals={}
    #   Get the list of all required keywords
    KeyList=KeyRTParams[:,0]
    #       Loop through all the keys in the RTParams imported file
    for x in dir(RTParams):
        #   Only grab the actual variables
        if not x.startswith('__'):
            #   Check if the variable is a required catalogue variable
            if x in KeyList:
                #   if it is, get the variable value
                VarUse=vars(RTParams)[x]
                #   Also get the target variable type
                ListLoc=np.where(KeyRTParams[:,0]==x)
                j=ListLoc[0][0]
                VarType=KeyRTParams[j][1]
                #   Check whether the type is correct
                if VarType == type(VarUse):
                    CatDict[x]=vars(RTParams)[x]
                #   If it isn't, exit the program
                else:
                    print("Required keyword", x, " is the wrong type")
                    print("The target type is", VarType)
                    print("exiting program")
                    exit()
            #   If the variable isn't required for that catalogue, it will be saved and used in the single galaxy fitting program
            else:
                NewGalaxyRTVals[x]=vars(RTParams)[x]
    #   Once done parsing everything, double check that all keywords are accounted for
    
    for key in KeyRTParams:
        keyU=key[0]
        try:
            RTVal=vars(RTParams)[keyU]
        except:
            print("Input file is missing required value for: ", keyU)
            print("Stopping run here")
            exit()

    return CatDict,NewGalaxyRTVals

def SupplementalIni(RTDict):
    #   Hardcode in the version of the code
    RTDict['KinVer']="WRKP V1"
    #   Set the name of the accepted model catalogue file
    RTDict['AcceptedModelCatalogueFile']=RTDict['KinTR']+"_KinematicModels.csv"
    return RTDict
