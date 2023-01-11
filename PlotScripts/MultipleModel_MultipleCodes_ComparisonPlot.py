import numpy as np
from decimal import Decimal
import argparse

import LoadModel as LM
import CubeAnalysis as CA
import DiagnosticPlotFncs as DPF
import LoadBaroloModel as LBM
import LoadFATModel as LFM

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator


def FileNames():
    ParentFolder="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT2/WRKP/WALLABY_NoiseTests/"
    BaroloResultsFolder="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/BaroloNoiseTests/"
    FATResultsFolder="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/FAT_NoiseTests/"
    
    
    TrueModelFolder="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/WALLABY_J100342-270137/"
    TrueModelLabel="Model"
    TrueModelFile=TrueModelFolder+"WALLABY_J100342-270137_AvgModel_v2.txt"

    
    TrueModel={'ModelFolder':TrueModelFolder,'ModelFile':TrueModelFile,'Label':TrueModelLabel,'DensFile':[None],'CodeSwitch':0}
    
    nCodes=3
    nModels=9
    
    ModelList=[None]*nModels*nCodes
    ModelList=np.reshape(ModelList,(nModels,nCodes))
    
    
    for i in range(nCodes):
        if i == 0:
            BaseCodeFolder=ParentFolder
            CodeBaseName="WALLABY_J100342-270137_Noise"
            ResultsBase=CodeBaseName
            CodeSwitch=0
            Label="WRKP"
        elif i ==1:
            BaseCodeFolder=BaroloResultsFolder
            CodeBaseName="TempModel_Noise"
            ResultsBase="rings_final1.txt"
            CodeSwitch=1
            Label="Barolo"
        elif i == 2:
            BaseCodeFolder=FATResultsFolder
            CodeBaseName="TempModel_Noise"
            ResultsBase="Finalmodel/Finalmodel.def"
            CodeSwitch=1
            Label="pyFAT"
            
        for j in range(nModels):
            FolderName=BaseCodeFolder+CodeBaseName+str(j+2)+"/"
            if i==0:
                ModelFile=FolderName+ResultsBase+str(j+2)+"_AvgModel_v1.txt"
                DensFile=[None]
            elif i==1:
                ModelFile=FolderName+ResultsBase
                DensFile=FolderName+"densprof.txt"
            elif i ==2:
                ModelFile=FolderName+ResultsBase
                DensFile=[None]
            ModelDict={'ModelFolder':FolderName,'ModelFile':ModelFile,'Label':Label,'DensFile':DensFile,'CodeSwitch':CodeSwitch}
            ModelList[j,i]=ModelDict
        
    FolderDict={'TrueModel':TrueModel,'AllModels':ModelList,'nCodes':nCodes,'nModels':nModels}
    return FolderDict


def Code_ParamComp(ModelsDict):
    BasePlotParams={'font.size': 15,'axes.linewidth':2
        ,'xtick.major.size':6,'xtick.minor.size':3
        ,'xtick.major.width':1,'xtick.minor.width':1
            ,'ytick.major.size':6,'ytick.minor.size':3
            ,'ytick.major.pad':10,'xtick.major.pad':10
            ,'ytick.major.width':1,'ytick.minor.width':1
            ,'xtick.labelsize':15 ,'ytick.labelsize':15
            ,'axes.labelsize': 20
            ,'legend.fontsize': 18
                }

    matplotlib.rcParams.update(BasePlotParams)
    
    
    PlotName="Noise_CodeComp_v1.png"
    
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))

    
    
    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    
    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS_FACEON"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    placement=[left,base+1*(h+buf),w,h]
    Key="INCLINATION"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    placement=[left+1*(w+buf),base+1*(h+buf),w,h]
    Key="POSITIONANGLE"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    placement=[left,base+0*(h+buf),w,h]
    Key="XCENTER"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    placement=[left+1*(w+buf),base+0*(h+buf),w,h]
    Key="YCENTER"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    placement=[left,base-1*(h+buf),w,h]
    Key="VSYS"
    DPF.KeywordPlot_CodeComp(fig,placement,Key,ModelsDict)
    
    
        #   Save the plot
    plt.savefig(PlotName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()
    
def CalcCenter_FromCube(RA,DEC,CubeDict):
    w=CubeDict['CubeWCS']
    X=[]
    Y=[]
    for i in range(len(RA)):
        RealCoords=[[RA[i],DEC[i],0]]
        PixCoords=w.wcs_world2pix(RealCoords,0)
        X.append(PixCoords[0,0])
        Y.append(PixCoords[0,1])
    return X,Y
    
def Main():
    #   Get the set of arguments needed to produce the model cube
    
    FileDict=FileNames()
    #print(FileDict)
   
    #First Load the true model
    TrueModel=LM.LoadBestFitModelFile(FileDict['TrueModel']['ModelFile'])
    TrueModel['SURFDENS_FACEON']=TrueModel['SURFDENS']*np.cos(TrueModel['INCLINATION']*np.pi/180.)
    TrueModel['SURFDENS_FACEON2']=TrueModel['SURFDENS_FACEON']
    
    nCodes=FileDict['nCodes']
    nModels=FileDict['nModels']
    ModelArr=[None]*nCodes*nModels
    ModelArr=np.reshape(ModelArr,(nModels,nCodes))
    
    SDConv=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    
    CubeName="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/TempModel_Noise1/TempModel_Noise1.fits"
    CubeInfo=CA.BasicCubeAnalysis(CubeName)
    
    for i in range(nCodes):
        for j in range(nModels):
            CurrModelFiles=FileDict['AllModels'][j][i]
            print(CurrModelFiles['ModelFile'])
            
            if i ==0:
                ModelArr[j,i]=LM.LoadBestFitModelFile(CurrModelFiles['ModelFile'])
                ModelArr[j,i]['SURFDENS_FACEON']=ModelArr[j,i]['SURFDENS']
                ModelArr[j,i]['SURFDENS_FACEON2']=ModelArr[j,i]['SURFDENS_FACEON']
                #print(j,i,ModelArr[j,i]['SURFDENS_FACEON'])
            elif i ==1:
                ModelArr[j,i]=LBM.LoadBaroloModel(CurrModelFiles['ModelFile'],CurrModelFiles['DensFile'])
                ModelArr[j,i]['SURFDENS_FACEON2']=ModelArr[j,i]['SURFDENS_FACEON']
                #print(j,i,ModelArr[j,i]['SURFDENS_FACEON'])
                #print(ModelArr[j,i]['SURFDENS_FACEON']*SDConv)
            elif i ==2 :
                ModelArr[j,i]=LFM.LoadFATModelFile(CurrModelFiles['ModelFile'])
                ModelArr[j,i]['SURFDENS_FACEON']/=SDConv
                ModelArr[j,i]['SURFDENS_FACEON2']/=SDConv
                
                ModelArr[j,i]['XCENTER'],ModelArr[j,i]['YCENTER']=CalcCenter_FromCube(ModelArr[j,i]['RA'],ModelArr[j,i]['DEC'],CubeInfo)
                #print(j,i,ModelArr[j,i]['SURFDENS_FACEON'],ModelArr[j,i]['SURFDENS_FACEON2'])
                
    ModelsDict={'TrueModel':TrueModel,'Models':ModelArr}
    
    Code_ParamComp(ModelsDict)
    
Main()

