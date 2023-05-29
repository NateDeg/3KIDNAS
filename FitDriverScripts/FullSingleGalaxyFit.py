import sys as sys
import os as os
import multiprocessing as mp
from multiprocessing import freeze_support

import time

from . import SoFiA_Driver as SD
from . import SetFileLocations as SFL
from . import GalaxyFitParameters as GFP
from . import RunWRKP as RW
from . import ReadWRKPFit as RWF
from . import MakeBootstrapSample as MBS
from . import Bootstrap_Error_Analysis as BEA
from . import Bootstrap_Outputs as BO
from . import GeometryCorrection as GC
from . import BootstrapModelPlot as BMP



def GalaxyFit():
    print("This is the WRKP Fit Driver Script")
    
    #   Get the various program and file locations
    GeneralDict,KeyRTParams=SFL.GetFileAndFolderLocationAndNames()
    #   Figure out the input parameter file
    ParamFile=GFP.GetRuntimeArguments()
    #   Read in the runtime parameters
    RTParams=GFP.ReadParamFile(ParamFile)
    #   Overwrite possible default parameters with the runtime parameters
    GeneralDict,GalaxyDict=GFP.OverwriteDefaults(GeneralDict,RTParams,KeyRTParams)
    #   Make sure all the variables are formatted correctly for use
    GFP.CheckParamTypes(GalaxyDict,KeyRTParams)
    #   Now use the default WRKP file locations that have already been specified and load them in
    GeneralDict=RW.LoadDefaultWRKPFiles(GeneralDict)
    
    #   Assuming that both the general dictionaries and galaxy dictionaries have the correct entries, we can run set up the WRKP input files
    #       Before beginning the analysis we'll copy the galaxy dictionary names to variables that will be used for a specific run.  This allows the running commands to be more general.
    GalaxyDict=GFP.SetCurrentRunVariables(GalaxyDict)
    
    #       The first step is create a 'catalogue' file with the initial estimates for the inclination and position angle
    GalaxyDict['SoFiAShapeFile']=SD.WriteSoFiACatFileForWRKP(GalaxyDict)
    #   Now run WRKP
    start = time.time() #   Time how long the fit takes
    GalaxyDict=RW.RunWRKP(GeneralDict,GalaxyDict,0)
    end = time.time()
    InitialFitTime=(end-start)
    #   Write out the total time of the fit
    TimeFile=GalaxyDict['TargFolderU']+GalaxyDict['ObjNameU']+"/"+GalaxyDict['ObjNameU']+"_FitTimeCheck.txt"
    f=open(TimeFile,'w')
    Str="Initial Fit Runtime (in seconds) = "+str(InitialFitTime)+"\n"
    f.write(Str)
    f.close()
    
    
    #   With WRKP now completed, the fit must be ingested
    #       It is possible that the initial fit will have failed.  In that case, we won't be able to read the output file, so we should end the program here
    try:
        GalaxyDict['BestFitModel']=RWF.ReadWRKPOutputFile(GeneralDict,GalaxyDict)
    except:
        print("No best fit model made")
        exit()
    
    
    #   The next step is to estimate the uncertainties for the model parameters by running a number of bootstrap fits.
    #   Before we loop through all the bootstraps, we need set a few naming conventions
    GalaxyDict=BEA.SetErrorAnalysisVariables(GalaxyDict)
    #       It's also necessary to load in the SoFiA template file
    GeneralDict['SoFiATemplate']=SD.LoadSoFiATemplate(GeneralDict['SoFiATemplateFile'])

    #       Now the loop over all bootstraps can be started
    BootstrapModels=[None]*GalaxyDict['nBootstraps']
    #           Set the number of processors to use for bootstrapping
    nProcessors=GalaxyDict['nProcessors_Bootstraps']
    #           Split the bootstrap fits over some number of parameters
    pool=mp.Pool(processes=nProcessors)
    BootstrapModels=pool.starmap(BootstrapRunStep, [(i,GeneralDict,GalaxyDict) for i in range(GalaxyDict['nBootstraps'])],chunksize=1)
    #for i in range(GalaxyDict['nBootstraps']):
    #    BootstrapModels[i]=BootstrapRunStep(i,GeneralDict,GalaxyDict)

 
    #   Get the uncertaintites on each parameter
    BEA.EstimateUncertaintiesFromBootstraps(GeneralDict,GalaxyDict,BootstrapModels)
    #   Adjust the PA to reach the global PA
    GalaxyDict['BestFitModel']=GC.GetGlobalPositionAngle(GalaxyDict)
    
    #   Output a text file with the errors and uncertainties
    BO.WriteBootstrappedFitOutputFile_Text(GalaxyDict)
    
    #   Make diagnostic plot from the bootstrapped model
    BMP.MakeBootstrapModelPlot(GalaxyDict,GeneralDict)

    #   Keep track of the total runtime
    end = time.time()
    TotalFitTime=(end-start)
    f=open(TimeFile,'a')
    Str="Total Fit Runtime (in seconds) = "+str(TotalFitTime)+"\n"
    f.write(Str)
    f.close()
    
    #   Once the bootstrap run is done, remove the bootstrap and SoFiA folders
    ClnCmd="rm -r "+GalaxyDict['BootstrapFolder']+" "+GalaxyDict['SoFiAFolder']
    os.system(ClnCmd)

def BootstrapRunStep(step,GeneralDict,GalaxyDict):
    BootstrapModel=BEA.GetBootstrapModel(GeneralDict,GalaxyDict,step)
    return BootstrapModel
