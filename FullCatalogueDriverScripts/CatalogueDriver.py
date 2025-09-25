import numpy as np
import pandas as pd
import os.path
from os import path
import multiprocessing as mp


from . import CatalogueFitParameters as CFP
from . import CatalogueInput as CI
from . import RunSingleGalaxyFit as RSGF
from . import CubePreprocessing as CP
from . import GenerateAcceptedModelCatalogue as GAMC




def CatalogueDriverMain():
    print("Doing full field analysis")
    
    #   Start by getting a dictionary of required catalogue variables and extra entries from a user-supplied python input file (set in the terminal)
    RTDict,SGDict=CFP.CatalogueRTImport()
    #   Do a few supplemental initializations
    RTDict=CFP.SupplementalIni(RTDict)
    #   Now load in the full source catalogue
    Cat=CI.LoadCatalogue(RTDict['CatName'])
    #   Sort the catalogue by size
    SizeSort=np.argsort(Cat['ell_maj'])
    Cat['SizeIndx']=SizeSort
    #   Now run the loop for the fitting
    i=0
   

    #nTot=1
    nTot=len(Cat)
    nProcessors=RTDict['nProcessors']
    pool=mp.Pool(processes=nProcessors,maxtasksperchild=1)
    pool.starmap(RunGalaxyFit, [(i,Cat,RTDict,SGDict) for i in range(nTot)],chunksize=1)
    
    #RunGalaxyFit(i,Cat,RTDict,SGDict)

    GAMC.GenerateAcceptedModelOutpouts(Cat,RTDict)
 
    


def RunGalaxyFit(step,Cat,RTDict,SGDict):
    #print(step)
    print("Process",mp.current_process(),step)

    #for i in range(len(Cat)):
    #    print(i,Cat['name'][i])
    Indx=Cat['SizeIndx'][step]

    #   Set up all the names that will be needed to write the input file
    GalaxyDict=RSGF.GetGalaxyDictNames(Indx,Cat,RTDict)
    #   Now do the cube pre-processing to get a velocity cube
    CP.CubePreprocessing(GalaxyDict)
    #   With the cube preprocessing completed, move on to getting the inclination and postion angle estimates for the cube
    GalaxyDict=RSGF.GetGeometryEstimates(Indx,Cat,GalaxyDict)
    #   With these in tow, the parameter file needed for the galaxy fit can be writtent
    RSGF.WriteSingleGalaxyIni(GalaxyDict,RTDict,SGDict)
    #   The single galaxy fitting code should be one directory up from this file so, figure out that path first
    CurrFile=__file__
    DriverPath=CurrFile.rsplit("/",2)[0]
    #   Use this path to get the fit driver
    FitDriver=DriverPath+"/WRKP_GalaxyFitDriver.py"
    #   Now run the fitter
    FitCmd="python3.9 "+FitDriver+" "+GalaxyDict['FitParameterFile']
    os.system(FitCmd)

        #   Finally clean things up
    MvCmd="mv "+GalaxyDict['VelCubeName'] + " "+RTDict['TargFolder']+GalaxyDict['name_underscore']+"/."
    os.system(MvCmd)
    MvCmd="mv "+GalaxyDict['FitParameterFile'] + " "+RTDict['TargFolder']+GalaxyDict['name_underscore']+"/."
    os.system(MvCmd)
   
    #print(ResultsFile,os.path.isfile(ResultsFile))
        #   If the fit was successful add the provenance values to the keywords
    ResultsFile=RTDict['TargFolder']+GalaxyDict['name_underscore']+"/"+GalaxyDict['name_underscore']+"_BSModel.txt"
    if os.path.isfile(ResultsFile):
        GAMC.AddProvenanceKeywords(GalaxyDict,RTDict)
    #else:
    #    RmCmd="rm "+GalaxyDict['VelCubeName']+" "+GalaxyDict['FitParameterFile']


