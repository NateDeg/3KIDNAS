import numpy as np
from decimal import Decimal
import argparse

import LoadModel as LM
import CubeAnalysis as CA
import DiagnosticPlotFncs as DPF

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator


def FileNames():
    ParentFolder="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT2/WRKP/WALLABY_NoiseTests/"
    
    nModels=3
    ModelFolders=[None]*nModels
    Labels=[None]*nModels
    IniParamFiles=[None]*nModels
    FirstPassFiles=[None]*nModels
    FinalFiles=[None]*nModels
    
    for i in range(nModels):
        if i ==0:
            #ModelFolders[i]="/Users/nate/Dropbox/WALLABY/DataReleases/Wallaby_NormaDR1V1_HydraDR2V2_KinematicModels/Wallaby_Hydra_DR2_KinematicModels_v2/WALLABY_J100342-270137/"
            ModelFolders[i]="/Users/nate/Dropbox/WALLABY/HydraDR2Analysis/Wallaby_Hydra_DR2_KinematicModels_TempTest_v2/WALLABY_J100342-270137/"
            Labels[i]="Model"
            FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_AvgModel_v2.txt"
            IniParamFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_AvgModel_v2.txt"
            FirstPassFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_AvgModel_v2.txt"
        elif i >=1:
            ModelFolders[i]=ParentFolder+"WALLABY_J100342-270137_Noise"+str(i)+"/"
            Labels[i]="Fit"
            
            IniParamFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_IniEstimate_v1.txt"
            FirstPassFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_FirstModel_v1.txt"
            FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_AvgModel_v1.txt"
            
       
        
        
    FolderDict=locals()
    return FolderDict
    
def ModelParamComp(ModelsDict):
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
    
    
    PlotName="Wallaby_J100342_EffectOfSDEstimate.png"
    
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
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left,base+1*(h+buf),w,h]
    Key="INCLINATION"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left+1*(w+buf),base+1*(h+buf),w,h]
    Key="POSITIONANGLE"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left,base+0*(h+buf),w,h]
    Key="XCENTER"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left+1*(w+buf),base+0*(h+buf),w,h]
    Key="YCENTER"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    placement=[left,base-1*(h+buf),w,h]
    Key="VSYS"
    DPF.KeywordPlot(fig,placement,Key,ModelsDict)
    
    

    #   Save the plot
    plt.savefig(PlotName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()


def Main():
    #   Get the set of arguments needed to produce the model cube
    
    FileDict=FileNames()
    #print(FileDict)
    
    IniModel=[None]*FileDict['nModels']
    FirstModel=[None]*FileDict['nModels']
    FinalModel=[None]*FileDict['nModels']
    
    for i in range(FileDict['nModels']):
        print(FileDict['IniParamFiles'][i])
        IniModel[i]=LM.LoadBestFitModelFile(FileDict['IniParamFiles'][i])
        FirstModel[i]=LM.LoadBestFitModelFile(FileDict['FirstPassFiles'][i])
        FinalModel[i]=LM.LoadBestFitModelFile(FileDict['FinalFiles'][i])
        print("hmmm",i, FinalModel[i]['R'])
    
    
   # print(FinalModel[0]['SURFDENS'])
    FinalModel[0]['SD_Obs']=FinalModel[0]['SURFDENS']
    FinalModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    FirstModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    IniModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    #print(FinalModel[0]['SURFDENS'])
    
    ModelsDict={'IniModel':IniModel,'FirstModel':FirstModel,'FinalModel':FinalModel,'nModel':FileDict['nModels'],'Labels':FileDict['Labels']}
    
    
    ModelParamComp(ModelsDict)
    
    
Main()

