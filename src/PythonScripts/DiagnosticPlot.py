import numpy as np
from decimal import Decimal
import argparse


import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

import AstroFncs as AF
import PlotFncs as PF
import CubeAnalysis as CA

def MakeDiagnosticPlot(AvgModel,DataCube,ModelCube,PlotName,GalaxyParams):
    print(PlotName)


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
    

    AnalysisFncs=AF.FuncDict()



    #   Open the plot
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    
    #   Make all the quick measurements from the cube -- they will be used for labels and in some subplots
    CubeMeasureDict=MakeQuickCubeMeasures(DataCube,AnalysisFncs,AvgModel,GalaxyParams)
    
    #   Write on all the labels due to the cube fit
    yTextLoc=AddObjectLabels(fig,CubeMeasureDict,AvgModel,GalaxyParams)
    
        #   Set the placement scheme
    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    
    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    PF.KeywordPlot_SingleModel(fig,placement,Key,AvgModel,CubeMeasureDict)
    
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    PF.KeywordPlot_SingleModel(fig,placement,Key,AvgModel,CubeMeasureDict)
    
    
    
    AvgPVPlots(fig,AvgModel,DataCube,ModelCube,left,base,w,h,buf)
    
    
    ScaleFactor=1.25
    wnew=w*ScaleFactor
    hnew=h*ScaleFactor
    left=left-(wnew-w)/2.
    base=base-(hnew-h)/1.5
    placement=[left-0.*(w+buf),base+1.0*(h+buf),wnew,hnew]
    Moment=0
    MomMaps=[]
    
    
    #Center=[AvgModel['XCENTER'][0],AvgModel['YCENTER'][0]]
    Center=[39.,39.]
    PixSize=np.abs(DataCube['CubeHeader']['CDELT1'])*3600.
    
    
    ax,MomMaps,MomExtent=PF.MomentPlot(fig,placement,Moment,DataCube,MomMaps,Center,PixSize)
    #ax,MomMaps,MomExtent=PF.MomentPlot(fig,placement,Moment,ModelCube,MomMaps,Center,PixSize)
    ax=PF.AddCenterToMomMap(ax,AvgModel['XCENTER'][0],AvgModel['YCENTER'][0])
    
  
    placement=[left+1.*(w+buf),base+1.*(h+buf),wnew,hnew]
    Moment=1
    ax,MomMaps,MomExtent=PF.MomentPlot(fig,placement,Moment,DataCube,MomMaps,Center,PixSize)
    #ax,MomMaps,MomExtent=PF.MomentPlot(fig,placement,Moment,ModelCube,MomMaps,Center,PixSize)

    ax=PF.AddCenterToMomMap(ax,AvgModel['XCENTER'][0],AvgModel['YCENTER'][0])
    
    
    #       Load in the average model moment map
    ModelVelMapFile=GalaxyParams['Folder']+"/"+GalaxyParams['Name']+"_Mom1_v"+GalaxyParams['Version']+".fits"
    #Mom1Header,Mom1Data,ModelMom1Map=CA.BasicCubeLoad(ModelVelMapFile)
    Mom1Header,Mom1Data=CA.BasicCubeLoad(ModelVelMapFile)
    
    #       Use the average model velocity map to draw model contours onto the observed moment map
    ax=PF.AddContoursToMomentPlot(ax,Moment,Mom1Data,AvgModel['VSYS'][0],MomMaps[1],DataCube['CubeVels'],MomExtent)
    
    #ax=PF.AddArrowToMomMap(ax,0.,0.,AvgModel['POSITIONANGLE'][0],0.5*CubeMeasureDict['Diameter']*30.)
    ax=PF.AddArrowToMomMap(ax,AvgModel['XCENTER'][0],AvgModel['YCENTER'][0],AvgModel['POSITIONANGLE'][0],0.5*CubeMeasureDict['Diameter'])
   
    #   Save the plot
    plt.savefig(PlotName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()


def AddObjectLabels(fig,CubeMeasureDict,AvgModel,GalaxyParams):
    #   First add the name of the plot
    TextSize=18
    fig.text(.6,.95, GalaxyParams['Name'] , ha='center',rotation=0,va='center',size=27)
    yTextLoc=0.70
    
    YTextStep=0.06
    MStr='{:0.2e}'.format(CubeMeasureDict['Mass'])
    LabelStr="M$_{HI}$\t\t= \t".expandtabs()+MStr+r" M$_{\odot}$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    
    #LabelStr="D$_{Hubble}$ \t=\t".expandtabs()+str(round(CubeMeasureDict['Distance'],1) )+" Mpc"
    #fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    #yTextLoc-=YTextStep
    
    
    
    #LabelStr="RHI$_{pred}$ \t\t=\t".expandtabs()+str(round(CubeMeasureDict['RHI'],1))+" '' "
    #fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    #yTextLoc-=(YTextStep+0.05)
    
    
    LabelStr="RA_kin \t\t=\t".expandtabs()+str(AvgModel['RA'][0])+" $\pm$ " +str(AvgModel['RA_ERR'][0])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep

    LabelStr="DEC_kin \t=\t".expandtabs()+str(AvgModel['DEC'][0])+" $\pm$ " +str(AvgModel['DEC_ERR'][0])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    
    LabelStr="Inc_kin \t=\t".expandtabs()+str(AvgModel['INCLINATION'][0])+" $\pm$ " +str(AvgModel['INCLINATION_ERR'][0])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    
    LabelStr="PA_kin \t\t=\t".expandtabs()+str(AvgModel['POSITIONANGLE'][0])+" $\pm$ " +str(AvgModel['POSITIONANGLE_ERR'][0])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    
    LabelStr="VSys_kin \t=\t".expandtabs()+str(AvgModel['VSYS'][0])+" $\pm$ " +str(AvgModel['VSYS_ERR'][0])+" km/s"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep

    return yTextLoc


def MakeQuickCubeMeasures(CubeInfo,AnalysisFncs,AvgModel,GalaxyParams):
    #   Set up an empty dictionary
    CubeMeasureDict={}
    #   Add the galaxy diameter, radius, and radius in beams...this assumes that the
    #       python script is passed the size in arcsecs
    CubeMeasureDict['Diameter']=float(GalaxyParams['Size'])/30.
    CubeMeasureDict['Radius']=CubeMeasureDict['Diameter']/2.*30.
    CubeMeasureDict['DiameterR']=CubeMeasureDict['Diameter']*30.

    #   Similarly, the script is passed the distance in Mpc
    CubeMeasureDict['Distance']=float(GalaxyParams['Distance'])
    #      And the mass in logarithm units
    CubeMeasureDict['logM']=float(GalaxyParams['logM'])
    CubeMeasureDict['Mass']=10.**CubeMeasureDict['logM']
    #   Get the predicted HI radius from the scaling relation found in AstroFncs.py
    RHITest=AnalysisFncs['CalcRHI'](CubeMeasureDict['logM'])
    #CubeMeasureDict['RHI']=AnalysisFncs['ProjectedSize_ArcSecond'](RHITest,CubeMeasureDict['Distance'])
    #   Finally set the central pixel from the model
    CubeMeasureDict['CentPix']=[AvgModel['XCENTER'][0],AvgModel['YCENTER'][0]]

    return CubeMeasureDict
    
    
    
def AvgPVPlots(fig,AvgModel,DataCube,ModelCube,left,base,w,h,buf):


    #   Make the major axis PV plot
    placement=[left+0.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=0
    ax=AvgModelPVPlot(fig,placement,AvgModel,DataCube,MajorMinorSwitch,ModelCube)
    #       For the major axis, add in the inclination corrected rotation curve
    PF.AddSingleRC_to_PVPlot(ax,AvgModel,DataCube)
    
    #   Make the minor axis PV plot
    placement=[left+1.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=1
    ax=AvgModelPVPlot(fig,placement,AvgModel,DataCube,MajorMinorSwitch,ModelCube)



def AvgModelPVPlot(fig,placement,AvgModel,DataCube,MajorMinorSwitch,ModelCube):

    "Make PV plot"
    
    #   First use the average model to find the center to be used for calculating the PV diagram
    CentPix=[AvgModel['XCENTER'][0],AvgModel['YCENTER'][0]]

    #   Now set the Position angle for the PV plot
    if MajorMinorSwitch==0:
        PAUse=AvgModel['POSITIONANGLE'][0]
        XLabel=r"Major Axis ('')"
        Center=AvgModel['XCENTER'][0]
    elif MajorMinorSwitch ==1:
        PAUse=AvgModel['POSITIONANGLE'][0]+90.
        XLabel=r"Minor Axis ('')"
        Center=AvgModel['YCENTER'][0]
    if PAUse > 360.:
        PAUse=PAUse-360.
    elif PAUse < 0.:
        PAUse=PAUse+360.
    print(CentPix,PAUse)
    #   The PV routine needs the beamsize in pixels to figure out how things should be cut
    BeamSize_Pix=DataCube['CubeHeader']['BMAJ']/np.abs(DataCube['CubeHeader']['CDELT1'])
    #print(BeamSize_Pix)
    #print(CubeInfo['CubeVels'])
    #   Construct the PV diagram for the data cube using the average geometry
    PV=CA.ConstructPVDiagram(DataCube['Data'],PAUse,CentPix,BeamSize_Pix,DataCube['CubeVels'])
    #print("PV sum check", np.nansum(DataCube['Data']),np.nansum(ModelCube['Data']))
    
    #   Now draw the datacube PV diagram
    PixSize=np.abs(DataCube['CubeHeader']['CDELT1'])*3600.
    PVDict={'PV':PV,'CubeVels':DataCube['CubeVels'],'PixSize':PixSize}
    ax=PF.BasePVPlot(fig,placement,PVDict)
    #       Add lines for vsys and the center point to the PV diagram
    #Center=np.shape(PVDict['PV'])[0]/2.
    ax=PF.AddCentLinesToPVPlot(ax,AvgModel['VSYS'][0],Center)
    
    #   Make a PV diagram of the model cube
    ModelPV=CA.ConstructPVDiagram(ModelCube['Data'],PAUse,CentPix,BeamSize_Pix,ModelCube['CubeVels'])
    PVDict={'PV':ModelPV,'CubeVels':ModelCube['CubeVels'],'PixSize':PixSize}
    #ax=AnalysisFncs['PlotFncs']['BasePVPlot'](fig,placement,PVDict)
    ax=PF.AddPVContoursToPlot(ax,PVDict)
    ax.set_xlabel(XLabel)
    return ax
    

