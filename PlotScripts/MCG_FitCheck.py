#!/usr/bin/env python3
import numpy as np

from decimal import Decimal


import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from matplotlib.patches import Ellipse

def SetTargs():
    #MCGModelFile="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT14/MCGSuite/MakeGalaxy_MGC_TestOutput/ba_5.5.mass_8.inc_45.0.pa_0.0.veldisp_8.0/MCG_TiltedRing.in"
    #BaroloPath="BaroloTrial/BaseBaroloTrial/"
    BaroloPath="BaroloTrial/Asymmetric_BaroloTrial_FixedINC_FixedPA/"
    MCGModelFile="/Users/nate/Dropbox/MockCubeGenerator_FourierMoments/MCGSuite/StronglyAsymmetric_Noise_1.6/ba_6.0.mass_9.5.inc_45.0.pa_0.0.veldisp_8.0.version_0/MCG_TiltedRing_0.in"
    
    MCGCubeFile="/Users/nate/Dropbox/MockCubeGenerator_FourierMoments/MCGSuite/StronglyAsymmetric_Noise_1.6/ba_6.0.mass_9.5.inc_45.0.pa_0.0.veldisp_8.0.version_0/ba_6.0.mass_9.5.inc_45.0.pa_0.0.veldisp_8.0.version_0_ConvolvedSourceCube.fits"
    
    
    
    BaroloModel=BaroloPath+"rings_final2.txt"
    BaroloModelDens=BaroloPath+"densprof.txt"
   
    #WRKPModel="SingleGalaxyFitTests/TestGalaxy/TiltedRingParams.txt"
    WRKPModel="AsymmetricGalaxyFitTests/TestGalaxy/TiltedRingParams.txt"
    
    #MCGCubeFile="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT14/MCGSuite/MakeGalaxy_MGC_TestOutput/ba_5.5.mass_8.inc_45.0.pa_0.0.veldisp_8.0/ba_5.5.mass_8.inc_45.0.pa_0.0.veldisp_8.0.fits"
    
    Targs={'MCGModel':MCGModelFile,'WRKPModel':WRKPModel,'MCGCube':MCGCubeFile,'BaroloModel':BaroloModel,'BaroloModelDens':BaroloModelDens}
    return Targs
    
def LoadMCGInput(MCGFile,MCGCube):
    R,RA,RA,DEC,Inc,PA,VSys,VRot,VRad,Vvert,VDisp,dvdz,Sigma,z0,zGrad=np.loadtxt(MCGFile,skiprows=10,unpack=True)
    #   To do some of the conversions we need the cube header
    MCGHDU=fits.open(MCGCube)
    MCG_WCS=wcs.WCS(MCGHDU[0].header)
    print(MCG_WCS.wcs.cdelt[2])
    MCGHDU.close()
    
    Sigma=Sigma/(MCG_WCS.wcs.cdelt[2]/1000.)
    crd=[[RA[0],DEC[0],VSys[0]]]
    pixcrd=MCG_WCS.wcs_world2pix(crd, 1)
    print(pixcrd)
    #   Get the
    
    X=R/R*pixcrd[0][0]
    Y=R/R*pixcrd[0][0]
    
    TRModel={'R':R,'R_SD':R,'XCENT':X,'YCENT':Y,'INCLINATION':Inc,'POSITIONANGLE':PA,'VSYS':VSys,'VROT':VRot,'VRAD':VRad,'SURFDENS':Sigma}
    return TRModel
    
def LoadWRKPOutput(WRKPFile):
    R,RW,X,Y,Inc,PA,VSys,VRot,VRad,Vvert,VDisp,dvdz,Sigma,z0,zGrad=np.loadtxt(WRKPFile,skiprows=12,unpack=True)
    
    #print("WRKPFile ",WRKPFile)
    #Sigma=Sigma/(30.**2.)
    Sigma=Sigma/(4.)
    TRModel={'R':R,'R_SD':R,'XCENT':X,'YCENT':Y,'INCLINATION':Inc,'POSITIONANGLE':PA,'VSYS':VSys,'VROT':VRot,'VRAD':VRad,'SURFDENS':Sigma}
    return TRModel
    
def LoadBaroloOutput(BaroloFile,BaroloDensFile):
    print("Loading Barolo Model")
    Params=np.loadtxt(BaroloFile,skiprows=1)
    RR,SD,SDErr=np.loadtxt(BaroloDensFile,skiprows=15,usecols=(0,9,10),unpack='True')
    R=Params[:,1]
    X=Params[:,9]
    Y=Params[:,10]
    Inc=Params[:,4]
    PA=Params[:,5]
    VSys=Params[:,11]
    VRot=Params[:,2]
    VRad=Params[:,12]
    Sigma=SD
    RR=RR
    
    
    TRModel={'R':R,'R_SD':RR,'XCENT':X,'YCENT':Y,'INCLINATION':Inc,'POSITIONANGLE':PA,'VSYS':VSys,'VROT':VRot,'VRAD':VRad,'SURFDENS':Sigma}
    return TRModel
    
def MakeCompPlot(TRDict):
    #   Name the output plot
    OutName="MCG_Barolo_Asymmetric_Comp.png"
    #   Set up the initial canvas
    Fw=10
    Fh=10
    fig=plt.figure(figsize=(Fw,Fh))
    #   Set up the plot placement values
    base=0.1
    left=0.1
    width=0.7
    height=width/2.
    
    #   Make the rotation curve plot
    RCbox=[left,base+height,width,height]
    Key="VROT"
    keywordPlot(fig,RCbox,TRDict,Key)
 
     #   Make the surface density plot
    SDbox=[left,base,width,height]
    Key="SURFDENS"
    keywordPlot(fig,SDbox,TRDict,Key)


    #   Add text comparisons
    #AddTextValues(fig,TRDict)

    plt.savefig(OutName, format='png',bbox_inches='tight')
    plt.close()

def keywordPlot(fig,box,TRDict,key):
    ax=fig.add_axes(box)
    for i in range(len(TRDict)):
        Col,DictSelect=SelectCol_and_Dict(i)
        if key=='SURFDENS':
            R=TRDict[DictSelect]['R_SD']
        else:
            R=TRDict[DictSelect]['R']
        ax.plot(R,TRDict[DictSelect][key],color=Col)
        
def SelectCol_and_Dict(step):
    if step==0:
        col='k'
        dictSelect="MCG"
    elif step==1:
        col='red'
        dictSelect="Barolo"
    elif step==2:
        col='blue'
        dictSelect="WRKP"
        
    elif step==3:
        col='green'
        dictSelect="WRKP2"
        
    elif step==4:
        col='magenta'
        dictSelect="WRKP3"
    elif step==5:
        col='cyan'
        dictSelect="WRKP4"
    return col,dictSelect
    
    
def AddTextValues(fig,TRDict):
       
    yTextLoc=0.6
    xTextLoc=[0.95,1.2]
    LabelStr="Input Values"
    fig.text(xTextLoc[0],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=25,color='k')
    
    LabelStr="WRKP Values"
    fig.text(xTextLoc[1],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=25,color='k')
    
    xTextLoc=[0.85,1.05,1.25,1.45]
    yTextLoc-=0.05
    key='XCENT'
    LabelStr="X"
    LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict)
    
    yTextLoc-=0.05
    key='YCENT'
    LabelStr="Y"
    LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict)
    
    yTextLoc-=0.05
    key='INCLINATION'
    LabelStr="Inc"
    LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict)
    
    yTextLoc-=0.05
    key='POSITIONANGLE'
    LabelStr="PA"
    LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict)
    
    yTextLoc-=0.05
    key='VSYS'
    LabelStr="VSys"
    LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict)
    
def LabelSet(fig,LabelStr,key,xTextLoc,yTextLoc,TRDict):
    fig.text(xTextLoc[0],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=20,color='k')
 
    LabelStr=str(round(TRDict['MCG'][key][0],2) )
    fig.text(xTextLoc[1],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=20,color='k')
    
    LabelStr=str(round(TRDict['Barolo'][key][0],2) )
    fig.text(xTextLoc[3],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=20,color='k')
    
    #LabelStr=str(round(TRDict['WRKP'][key][0],2) )
    #fig.text(xTextLoc[2],yTextLoc, LabelStr,ha='left',rotation=0,va='center',size=20,color='k')
   
   

   
def Main():
    print("Comparing WRKP result to MCG model")
    #   Get the target files
    CompFiles=SetTargs()
    #   Load the input model
    MCG_TRModel=LoadMCGInput(CompFiles['MCGModel'],CompFiles['MCGCube'])
    #   Load the WRKP output
    WRKP_TRModel=LoadWRKPOutput(CompFiles['WRKPModel'])
    #WRKP_TRModel2=LoadWRKPOutput(CompFiles['WRKPModel2'])
    #WRKP_TRModel3=LoadWRKPOutput(CompFiles['WRKPModel3'])
    #WRKP_TRModel4=LoadWRKPOutput(CompFiles['WRKPModel4'])
    #   Load in a barolo model
    Barolo_TRModel=LoadBaroloOutput(CompFiles['BaroloModel'],CompFiles['BaroloModelDens'])
    #   Make the comparison plot
    TRDict={'MCG':MCG_TRModel,'Barolo':Barolo_TRModel,'WRKP':WRKP_TRModel}
    
    MakeCompPlot(TRDict)
    
    






Main()
