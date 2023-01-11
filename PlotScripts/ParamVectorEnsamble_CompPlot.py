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
    
    nModels=2
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
            
        elif i ==1:
            #ModelFolders[i]=ParentFolder+"WALLABY_J100342-270137_Model_NewIniRange/"
            #Labels[i]="Fit"
            
            #IniParamFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Model_NewIniRange_IniEstimate_v1.txt"
            #FirstPassFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Model_NewIniRange_FirstModel_v1.txt"
            #FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Model_NewIniRange_AvgModel_v1.txt"
            ModelFolders[i]=ParentFolder+"WALLABY_J100342-270137_Noise"+str(i)+"/"
            Labels[i]="Fit"
            
            IniParamFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_IniEstimate_v1.txt"
            FirstPassFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_FirstModel_v1.txt"
            FinalFiles[i]=ModelFolders[i]+"WALLABY_J100342-270137_Noise"+str(i)+"_AvgModel_v1.txt"
        
    EnsambleFile="AllFinalParamVectors.txt"
    IniEnsambleFile="AllInitialParamVectors.txt"
    FolderDict=locals()
    return FolderDict
    
def ModelParamComp(ModelsDict,EnsambleModels):
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
    
    
    PlotName="Wallaby_J100342_NoiseTest1_Ensamble.png"
    
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
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
   
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
    
    placement=[left,base+1*(h+buf),w,h]
    Key="INCLINATION"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
    
    placement=[left+1*(w+buf),base+1*(h+buf),w,h]
    Key="POSITIONANGLE"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
    
    placement=[left,base+0*(h+buf),w,h]
    Key="XCENTER"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
    
    placement=[left+1*(w+buf),base+0*(h+buf),w,h]
    Key="YCENTER"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
    
    placement=[left,base-1*(h+buf),w,h]
    Key="VSYS"
    DPF.KeywordPlot_Ensamble(fig,placement,Key,ModelsDict,EnsambleModels)
   
    

    #   Save the plot
    plt.savefig(PlotName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()

def ParamLineToModel(Line,Model):
    nR=len(Model['R'])
    Vals=Line.split()
    Vals=np.asarray(Vals, dtype=np.float64, order='C')
    X=np.array([Vals[0]]*nR)
    Y=np.array([Vals[1]]*nR)
    Inc=np.array([Vals[2]]*nR)*180./np.pi
    PA=np.array([Vals[3]]*nR)*180./np.pi-90.
    VSys=np.array([Vals[4]]*nR)
    
    Start=5
    Fin=5+nR
    VRot=np.array(Vals[Start:Fin])
    SD=np.array(Vals[Fin:])
    SDConv1=3.93346533/36.
    SDConv2=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    SD=SD*SDConv1/SDConv2
    print(SD)

    
    #print(nR, Model['R'])
    #print(VRot,len(VRot))
    #print(SD,len(SD))
    
    

    Model['XCENTER']=X
    Model['YCENTER']=Y
    Model['INCLINATION']=Inc
    Model['POSITIONANGLE']=PA
    Model['VSYS']=VSys
    Model['VROT']=VRot
    Model['SURFDENS']=SD
    
    return Model

def LoadParamFile(FileName,FinalModel):
    print("Load Param File")
    
    f=open(FileName,"r")
    Lines=f.readlines()
    print(len(Lines))

    nModels=len(Lines)
    Models=np.array([None]*nModels)
    for i in range(nModels):
        Models[i]={'R':FinalModel['R'],'R_SD':FinalModel['R']}
        Models[i]=ParamLineToModel(Lines[i],Models[i])
        #print(FinalModel['VROT'])
        print(Models[i]['VROT'])
    f.close()
    return Models

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
        #print("hmmm",i, FinalModel[i]['R'])
    
    
   # print(FinalModel[0]['SURFDENS'])
    FinalModel[0]['SD_Obs']=FinalModel[0]['SURFDENS']
    FinalModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    FirstModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    IniModel[0]['SURFDENS']*=np.cos(FinalModel[0]['INCLINATION']*np.pi/180.)
    #print(FinalModel[0]['SURFDENS'])
    
    EnsambleModels=LoadParamFile(FileDict['EnsambleFile'],FinalModel[1])
    IniEnsambleModels=LoadParamFile(FileDict['IniEnsambleFile'],FinalModel[1])
    FullEnsambleModels=[IniEnsambleModels,EnsambleModels]
    
    ModelsDict={'IniModel':IniModel,'FirstModel':FirstModel,'FinalModel':FinalModel,'nModel':FileDict['nModels'],'Labels':FileDict['Labels']}
    
    
    ModelParamComp(ModelsDict,FullEnsambleModels)
    
    SDConv1=3.93346533/36.
    SDConv2=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    SD=15.
    SDTemp=SD/SDConv1*SDConv2
    print(SD,SDTemp)
    
Main()

