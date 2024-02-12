import os as os

def GetFileAndFolderLocationAndNames():
    #   Get the current directory for future writing
    CurrDir=os.getcwd()
    #   Get the WRKP directory
    WRKPDir= os.path.dirname(os.path.realpath(__file__))
    #   Adjust the WRKP directory path slightly to point one directory up due to this file being located in WRKP/FitDriverScripts
    WRKPDir=WRKPDir.rsplit('/', 1)[0]
    #   Set the SingleGalaxyFitter and BootstrapSampler paths
    FitterExecPath=WRKPDir+"/Programs/SingleGalaxyFitter"
    BootstrapExecPath=WRKPDir+"/Programs/BootStrapSampler"
    #   Set the SoFiA path to the third_party install
    SoFiAExecPath=WRKPDir+"/third_party/SoFiA-2-master_2_5_1/sofia"
    #   Set the default SoFiA template file
    SoFiATemplateFile=WRKPDir+"/third_party/SoFiA-2-master_2_5_1/template_par_file.par"
    #   Set the default WRKP Main file
    WRKP_GeneralMainIn=WRKPDir+"/Inputs/SingleFitInput_Base.in"
    #   Set the default WRKP Options file
    WRKP_GeneralOptionsIn=WRKPDir+"/Inputs/SingleGalaxyTestFittingOptions_Base.txt"
    #   Set the general file dictionary as the set of local variables thus far
    FileDict=locals()
    #   Set the list of runtime parameters that must be set in the RT Params  file
    #       It is necessary to list each variable with it's type i.e. [Variable Name, Type] as the types will be checked in GalaxyFitParamters.CheckParamTypes
    KeyRTParamsList=[['CubeName',str]
    ,['MaskName',str]
    ,['PA_Estimate',float]
    ,['Inc_Estimate',float]
    ,['nBootstraps',int]
    ,['ObjName',str]
    ,['TargFolder',str]
    ,['nProcessors_Bootstraps',int]]
    return FileDict,KeyRTParamsList
    
