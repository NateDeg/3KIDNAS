import numpy as np
import pandas as pd
import os.path
from os import path
import sys as sys
import multiprocessing as mp

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from datetime import date


from . import ModelAcceptance as MA



def GetModelNames(FileDict,GalName):

    #   Set the name of the galaxy
    ModelNames={'GalName':GalName}
    #   And the name of the galaxy with underscores
    ModelNames['GalName_Underscore']=GalName.replace(' ','_')
    #   Set the name of the output folder
    ModelNames['ResultsFolder']=FileDict['TargFolder']+ModelNames['GalName_Underscore']+"/"
    #   Set the name of the model file
    ModelNames['ModelFile']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_BSModel.txt"
        #   Set the name of the model plot
    ModelNames['ModelPlot']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_BSModel.png"
        #   Set the name of the time check file
    ModelNames['TimeFile']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_FitTimeCheck.txt"
    #   Set the name of the cube file
    ModelNames['CubeFile']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_AverageModel_v1.fits"
    #   Set the name of flags file
    ModelNames['FlagFile']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_Flags_v1.txt"
    #   Set the name of the difference cube
    ModelNames['DiffCube']=ModelNames['ResultsFolder']+"DifferenceCube.fits"
    #   Set the name of the processed cube
    ModelNames['ProcCube']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_VelCube.fits"
    #   Set the name of the PVMaps
    ModelNames['PVMap1']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_PVMajor_Data.fits"
    ModelNames['PVMap2']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_PVMinor_Data.fits"
    ModelNames['PVMap3']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_PVMajor_Model.fits"
    ModelNames['PVMap4']=ModelNames['ResultsFolder']+ModelNames['GalName_Underscore']+"_PVMinor_Model.fits"
    
    StoreKeys=['CubeFile','ModelFile','FlagFile','DiffCube','ModelPlot','ProcCube','PVMap1','PVMap2','PVMap3','PVMap4']
    TargKeys=['TargCubeFile','TargModelFile','TargFlagFile','TargDiffCube','TargModPlot','TargProcCube','TargPVMap1','TargPVMap2','TargPVMap3','TargPVMap4']
    TargNames=['ModCube.fits','AvgMod.txt','Flag.txt','DiffCube.fits','DiagnosticPlot.png','ProcData.fits','PVMajor_Data.fits','PVMinor_Data.fits','PVMajor_Model.fits','PVMinor_Model.fits']
    j=0
    for key in TargKeys:
        ModelNames[key]=ModelNames['GalName_Underscore']+"_"+FileDict['KinTR']+"_"+TargNames[j]
        j+=1
    ModelNames['StoreKeys']=StoreKeys
    ModelNames['TargKeys']=TargKeys
    
    ModelNames['SourceFolder']=FileDict['SourceFolder']+ModelNames['GalName_Underscore']+"/"
    SourceFiles=os.listdir(ModelNames['SourceFolder'])
    subs = 'cube.fits'
    ModelNames['OriCube']=ModelNames['SourceFolder']+ [i for i in SourceFiles if subs in i][0]
    subs = 'mask.fits'
    ModelNames['OriMask']=ModelNames['SourceFolder']+ [i for i in SourceFiles if subs in i][0]

    return ModelNames
    
        
    
def LoadModelOutputs(ModelNames,FileDict):
    #   Check if the model finished
    ModelNames['ModelFitAchieved']=os.path.exists(ModelNames['ModelFile'])
    #print(ModelNames['ModelFile'],os.path.exists(ModelNames['ModelFile']))
    #   Read the time file
    try:
        ModelNames=ReadTimeFile(ModelNames)
    except:
        ModelNames['T_IniFit']=-1

    if ModelNames['ModelFitAchieved']:
        import ReadWRKPFit as RPF
        #   Read the model file
        ModelNames['Model']=RPF.LoadBootstrappedFit(ModelNames['ModelFile'])
    #print(ModelNames['ModelFitAchieved'],ModelNames['ModelFile'])
    return ModelNames
    
    
def CopySuccessfulModel(ModelNames,FileDict):
    #   First make the kinematic model folder
    os.makedirs(FileDict['KinFolder'], exist_ok=True)
    #   Also make a folder for all accepted kinematic model diagnostic plots
    os.makedirs(FileDict['PlotFolder'], exist_ok=True)
    #   First copy the model plot to the overarching model folder
    CpCmd="cp "+ModelNames['ModelPlot'] + " "+FileDict['PlotFolder']+"."
    os.system(CpCmd)
    #   Now make a directory for the model results
    AcceptedModelFolder=FileDict['KinFolder']+ModelNames['GalName_Underscore']+"/"
    os.makedirs(AcceptedModelFolder, exist_ok=True)
    #   Set the keys of files to copy
    j=0
    for key in ModelNames['StoreKeys']:
        CpCmd="cp "+ModelNames[key]+" "+AcceptedModelFolder+ModelNames[ModelNames['TargKeys'][j]]
        j+=1
        os.system(CpCmd)
        
def IniResultsDict(SrcCat,KeyParams,SuccesParams,ProfParams):
    ResultsDict={}
    for x in KeyParams:
        ResultsDict[x]=[None]*len(SrcCat)
    for x in SuccesParams:
        ResultsDict[x]=[None]*len(SrcCat)
    for x in ProfParams:
        ResultsDict[x]=[None]*len(SrcCat)
    print(ResultsDict.keys())
    
    return ResultsDict

def GetDriverFolder():
    ModFile=(os.path.realpath(__file__))
    CatDriverDir=ModFile.rsplit('/', 1)[0]
    PipelineDir=CatDriverDir.rsplit('/', 1)[0]
    GalFitterDir=PipelineDir+"/FitDriverScripts/"
    return locals()
    

def GenerateAcceptedModelOutpouts(Cat,RTDict):
    print("Checking Results")

    #   Get the Cut limits
    CutLimits=MA.SetCutLimits()
    #   Set the kinematic models folder
    #RTDict['KinTR']=Cat['team_release'][0]+"_KinTR1"
    KinName=RTDict['TargFolder'].strip("/")+"_"+RTDict['KinTR']
    RTDict['KinFolder']=RTDict['TargFolder']+KinName+"/"
    #   Make a folder to store a copy of all the diagnostic plots to make life simpler for checking everything
    RTDict['PlotFolder']=RTDict['KinFolder']+"DiagnosticPlots/"
    #   Also make a folder to store all diagnostic plots regardless of success
    RTDict['AllPlotsFolder']=RTDict['TargFolder']+"DiagnosticPlots/"
    #print("All Plots Folder",RTDict['AllPlotsFolder'])
    #   Make the model and plot folder inside the fitted folder
    os.makedirs(RTDict['PlotFolder'], exist_ok=True)
    #   We'll also want to make the diagnostic plot folder
    os.makedirs(RTDict['AllPlotsFolder'],exist_ok=True)
    
    #   We'll need to figure out the path of the current folder
    ModDict=GetDriverFolder()
    sys.path.append(ModDict['GalFitterDir'])
    #   Now set up an empty dictionary to store the needed parameters


    KeyParams=['name']
    ModelParamNames=['RMS', 'SN_Integrated', 'SN_Peak', 'SN_Avg', 'SN_Median','XCENTER', 'XCENTER_ERR', 'YCENTER', 'YCENTER_ERR', 'INCLINATION', 'INCLINATION_ERR', 'POSITIONANGLE', 'POSITIONANGLE_ERR', 'VSYS', 'VSYS_ERR', 'RA', 'RA_ERR', 'DEC', 'DEC_ERR', 'VDISP', 'VDISP_ERR']
    SuccesParams=['RMS', 'SN_Integrated', 'SN_Peak', 'SN_Avg', 'SN_Median','X_model', 'e_X_model', 'Y_model',
       'e_Y_model', 'Inc_model', 'e_Inc_model', 'PA_model', 'e_PA_model', 'Vsys_model', 'e_Vsys_model', 'RA_model', 'e_RA_model', 'DEC_model', 'e_DEC_model', 'Vdisp_model', 'e_Vdisp_model']
    ProfParams=['Rad', 'Vrot_model', 'e_Vrot_model', 'Rad_SD', 'SD_model', 'e_SD_model']
    ProfParamNames=['R','VROT', 'VROT_ERR', 'R_SD','SURFDENS', 'SURFDENS_ERR']

    
    ResultsDict=IniResultsDict(Cat,KeyParams,SuccesParams,ProfParams)
    #print("Results Dictionary", ResultsDict)
    


    AutoAccepted=0
    for i in range(len(Cat)):
    #for i in range(10):
        #print(i,Cat['name'][i])

        #   Get the names for the model files
        ModelNames=GetModelNames(RTDict,Cat['name'][i])

        #   Load in the model results
        Results=LoadModelOutputs(ModelNames,RTDict)
        #   If the model plot exists, copy it to the all plots folder
        if os.path.isfile(Results['ModelPlot']):
            CpCmd="cp "+Results['ModelPlot']+" "+RTDict['AllPlotsFolder']
            os.system(CpCmd)
        
        #   Get the cube header so that we can get the SoFiA size estimate in beams
        Cube=fits.open(ModelNames['OriCube'])
        CHeader=Cube[0].header
        Cube.close()
        BeamSize_Pixels=np.abs(CHeader['BMAJ']/CHeader['CDELT1'])
        

        #   Add some of the SoFiA estimates to the model dictionary for checking
        Results['Size']=Cat['ell_maj'][i]/BeamSize_Pixels
        Results['Inc_SoFiA']=np.arccos(Cat['ell_min'][i]/Cat['ell_maj'][i])
        Results['PA_SoFiA']=Cat['kin_pa'][i]
        Results['name']=Cat['name'][i]
        Results['ell_maj_SoFiA']=Cat['ell_maj'][i]
        Results['ell_min_SoFiA']=Cat['ell_min'][i]

        #   See whether the model fits our success rules
        Results=MA.DetermineSuccess(Results,CutLimits,ModelNames,BeamSize_Pixels)

        #   If the model is successfully, then we'll want to set up the output folder properly and add it to a dictionary
        if Results['ModelSuccess']==1:
            Results['Model']['KinTR']=RTDict['KinTR']
            k=AutoAccepted
            #   Copy the appropriate files to a new folder
            CopySuccessfulModel(ModelNames,RTDict)
            
            for x in KeyParams:
                ResultsDict[x][k]=Results[x]
            #if Results['ModelSuccess']==1:
            #    Results['Model']['T_FinFit']=np.array([Results['T_FinFit']])

            for j in range(len(SuccesParams)):
                x=SuccesParams[j]
                y=ModelParamNames[j]
                ResultsDict[x][k]=Results['Model'][y][0]
            j=0
            for x in ProfParams:
                y=ProfParamNames[j]
                j+=1
                Arr=Results['Model'][y]
                Str=', '.join([str(Val) for Val in Arr])
                ResultsDict[x][k]=Str

            AutoAccepted+=1
            
    print("Total Automatically Accepted Models",AutoAccepted)

    #   Once done, we want to trim the final dictionary to only the accepted elements
    for key in ResultsDict.keys():
        ResultsDict[key]=np.delete(ResultsDict[key],slice(AutoAccepted,len(ResultsDict[key])))

    
    DF=pd.DataFrame.from_dict(ResultsDict)
    #   Name the CSV file
    RTDict['AcceptedModelCatalogueFile']=RTDict['KinFolder']+RTDict['AcceptedModelCatalogueFile']
    DF.to_csv(RTDict['AcceptedModelCatalogueFile'],index=False)

def AddProvenanceKeywords(GalaxyDict,RTDict):
    #       Fix the processed cube keywords
    #       Start by getting the current folder and cube name
    CurrFolder=RTDict['TargFolder']+GalaxyDict['name_underscore']+"/"
    #   Set the processed cube name
    ProcCubeName=CurrFolder+GalaxyDict['VelCubeName']
    #       Open the processed cube
    ProcCube=fits.open(ProcCubeName,mode='update')
    #   Use the specific element of the cube HDU
    HeaderID=0
    #   Use the header to set the SRCVERS keywords for other fits files
    SRCVers=ProcCube[HeaderID].header['ORIGIN']
    #   Get rid of the SoFiA history here
    try:
        del ProcCube[HeaderID].header['HISTORY']
    except:
        print("No History keyword to delete in cube header")
    #   Specifically add the SRCVERS, KINVERS, and KINTR keywords
    #       Start with setting the provenance keyword dictionary
    ProvDict={"SRCVER":SRCVers,"KINVER":RTDict['KinVer'],"KINTR":RTDict['KinTR']}
    for key in ProvDict:
        ProcCube[HeaderID].header.set(key,ProvDict[key])
    #   Add the date
    now = date.today()
    Date= str(now.year) + "-" + str(now.month) + "-" +  str(now.day)
    ProcCube[HeaderID].header.set('DATE',Date)
    #   Save the current cube
    ProcCube.flush()
    #   Close the current cube
    ProcCube.close()
    
    #   Set the name of the other two cubes to add provenance to
    ModCube=CurrFolder+GalaxyDict['name_underscore']+"_AverageModel_v1.fits"
    DiffCube=CurrFolder+"DifferenceCube.fits"
    #   PVCubes
    PVMap1=CurrFolder+GalaxyDict['name_underscore']+"_PVMajor_Data.fits"
    PVMap2=CurrFolder+GalaxyDict['name_underscore']+"_PVMinor_Data.fits"
    PVMap3=CurrFolder+GalaxyDict['name_underscore']+"_PVMajor_Model.fits"
    PVMap4=CurrFolder+GalaxyDict['name_underscore']+"_PVMinor_Model.fits"
    #   Add provenance keywords model and difference cubes
    CubeNames=[ModCube,DiffCube,PVMap1,PVMap2,PVMap3,PVMap4]
    for Name in CubeNames:
        #   Open the cube
        Cube=fits.open(Name,mode='update')
        #   Add the provenance keywords
        for key in ProvDict:
            Cube[HeaderID].header.set(key,ProvDict[key])
        #   Add the date
        Cube[HeaderID].header.set('DATE',Date)
            #   Save the current cube
        Cube.flush()
        #   Close the current cube
        Cube.close()

