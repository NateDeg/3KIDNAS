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
from astropy.coordinates import Angle


from . import CubeAnalysis as CA
from . import PVPlotFncs as PVP
from . import MomentMapPlotFncs as MMP

"""
    This module contains the routines needed to make the average model diagnostic plot for a specific galaxy.  It contains the routines:
    MakeAvgModelPlot --> This makes and saves the diagnostic plot for the 'average' model.
    NameComparisonPlot --> This function names the plot
    AddObjectLabels --> This function adds a set of labels to the plot.
    AvgPVPlots --> This function makes the pair of PV diagrams needed for the comparison plot.
"""

def MakeBootstrapModelPlot(GalaxyDict):
    """
        This function makes the diagnostic plot for the average Tilted ring model.
    """
    #   Start by setting the plot parameters to some useful defaults
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
    #   Update the default parameters to the new ones.
    matplotlib.rcParams.update(BasePlotParams)
    #   We will need to the best fitting model
    Model=GalaxyDict['BestFitModel']
    #   And we'll need both the real data cube and the best fitting model cube

    
    DataCube=CA.BasicCubeAnalysis(GalaxyDict['CubeName'])
    ModelCube=CA.BasicCubeAnalysis(Model['ModelCube'])
    #   Name the plot
    PltName=GalaxyDict['WRKP_ResultsFolder']+GalaxyDict['ObjName']+"_BSModel.png"
    #   Open the plot
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))

    #   Write on all the labels due to the cube fit
    yTextLoc=AddObjectLabels(fig,Model,GalaxyDict)


    #   Set the placement parameters that are used to organize all the different panels in the diagnostic plot
    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    #   Add the keyword plots (see DiagnosticPlots/DiagnosticPlotFuncs.py)
    #       Add the rotation curve panel
    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    KeywordPlot(fig,placement,Key,Model)

    #       Add the surface density panel
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS_FACEON"
    KeywordPlot(fig,placement,Key,Model)
    
    #       Add the pair of PV diagram plots
    AvgPVPlots(fig,Model,DataCube,ModelCube,left,base,w,h,buf)
  
    #   Rescale the sizes/boxes for the moment maps
    ScaleFactor=1.25
    wnew=w*ScaleFactor
    hnew=h*ScaleFactor
    left=left-(wnew-w)/2.
    base=base-(hnew-h)/1.5
    #   Make the moment 0 map
    placement=[left-0.*(w+buf),base+1.0*(h+buf),wnew,hnew]
    Moment=0
    #       For the moment maps, it's necessary to give it the model center in pixels to center the maps
    XC=Model['XCENTER'][0]
    YC=Model['YCENTER'][0]
   
    #   Call the moment map function to make the basic map
    ax,MomData=MMP.MomentPlot(fig,placement,Moment,DataCube,XC,YC,GalaxyDict)
    #   Add a marker to the center of the map
    ax=MMP.AddCenterToMomMap(ax)
    placement=[left+1.*(w+buf),base+1.*(h+buf),wnew,hnew]
   
    #   Make the moment 1 map
    Moment=1
    ax,MomData=MMP.MomentPlot(fig,placement,Moment,DataCube,XC,YC,GalaxyDict)
    #   Add the center to moment 1 map
    ax=MMP.AddCenterToMomMap(ax)
    #   Add the model velocity contours to the observed moment 1 map
    ax=MMP.AddVelContoursToMomentPlot(ax,ModelCube,Model['VSYS'][0],MomData,XC,YC,Model)
    #   Add an arrow to the PV diagram indicating the position angle of the model
    ax=MMP.AddArrowToMomMap(ax,Model['POSITIONANGLE'][0],0.5*Model['R'][-1])
    
    #   Save the plot
    plt.savefig(PltName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()


def AddObjectLabels(fig,Model,GalaxyDict):
    """
        This function adds a number of labels to the plot
    """
    #   First add the name of the plot
    TextSize=18
    fig.text(.6,.95, GalaxyDict['ObjName'] , ha='center',rotation=0,va='center',size=27)
    #   Set a top location
    yTextLoc=-0.25
    #   Set size of the vertical steps for the labels
    YTextStep=0.06
    #   Set the horizontal location of the labels
    XTextLoc=0.35
    
    #   Convert the RA into hours:minutes:seconds
    RA_S=RA_DEC_Str(Model['RA'][0],Model['RA_ERR'][0],0)
    #   Add the RA label to the plot
    LabelStr="RA_model \t=\t".expandtabs()+RA_S
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep

    #   Convert the DEC into degrees:arcminutes:arcseconds
    DEC_S=RA_DEC_Str(Model['DEC'][0],Model['DEC_ERR'][0],1)
    #   Add the DEC string to the plot
    LabelStr="DEC_model \t=\t".expandtabs()+DEC_S
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the model inclination and error to the plot
    LabelStr="Inc_model \t=\t".expandtabs()+str(Model['INCLINATION'][0])+" $\pm$ " +str(Model['INCLINATION_ERR'][0])+" $^\circ$"
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the global position-angle to the plot along with it's uncertainty.
    LabelStr="PA_model_g \t=\t".expandtabs()+str(Model['PA_GLOBAL'])+" $\pm$ " +str(Model['PA_GLOBAL_ERR'])+" $^\circ$"
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the systemic velocity to the plot
    LabelStr="VSys_model \t=\t".expandtabs()+str(Model['VSYS'][0])+" $\pm$ " +str(Model['VSYS_ERR'][0])+" km/s"
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    
    #   Add the velocity dispersion to the plot
    LabelStr="VDisp_model \t=\t".expandtabs()+str(Model['VDISP'][0])+" $\pm$ " +str(Model['VDISP_ERR'][0])+" km/s"
    fig.text(XTextLoc,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    return yTextLoc

    
def AvgPVPlots(fig,Model,DataCube,ModelCube,left,base,w,h,buf):
    """
        This function draws the major and minor axis PV panels onto the diagnostic plot.  It uses routines found in DiagnosticPlots/PVPlotFuncs.py
    """
    #   For the PV plots we need to load in the full resolution MCG model cubelet
   
  
    #   Make the major axis PV plot
    placement=[left+0.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=0
    ax=PVP.AvgModelPVPlot(fig,placement,Model,DataCube,MajorMinorSwitch,ModelCube)
    
    #       For the major axis, add in the inclination corrected rotation curve
    PVP.AddSingleRC_to_PVPlot(ax,Model,DataCube)
    #   Make the minor axis PV plot
    placement=[left+1.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=1
    ax=PVP.AvgModelPVPlot(fig,placement,Model,DataCube,MajorMinorSwitch,ModelCube)
    


def RA_DEC_Str(ValU,Err,DType):
    
    #   Convert the string into an angle
    Val_Test=Angle(ValU,u.deg)
    if DType==0:
        #   For RA, use hours, minutes, seconds
        Val_Fmt=Val_Test.hms
        LabelS=["h","m","S"]
    elif DType==1:
        #   For DEC, use degrees, arcminutes, arcseconds
        Val_Fmt=Val_Test.signed_dms
        LabelS=[r"$^{\circ}$","'","''"]
        
    #   Format the string
    ValStr=""
    for i in range(len(Val_Fmt)):
        Val=Val_Fmt[i]
        if Val < 0:
            ValS="-"
        if i < 2:
            Val=round(Val,0)
            ValS=str(Val).split('.')[0]
        else:
            Val=round(Val,1)
            ValS=str(Val)
        ValStr+=ValS
        if i==0:
            ValStr+=LabelS[0]
        elif i == 1:
            ValStr+=LabelS[1]
        elif i ==2:
            ValStr+=LabelS[2]
    #   Now deal with the uncertainty
    Err_A=Angle(Err,u.deg)
    #RA_Err_Str=RA_Err.to_string(unit=u.degree)
    Val_Err_Str=Err_A.dms
    Val_Err_Str=str(round(Val_Err_Str[2],1))
    Val_Err_Str=Val_Err_Str+"''"
    #       Do the same for the uncertainty in RA and turn it into arcseconds
    FullStr=ValStr+" $\pm$ " +Val_Err_Str
    return FullStr



def KeywordPlot(fig,placement,Key,Model):
    """
        This function makes a plot showing a specific model's values for a given parameter indicated by the keyword
    """
    #   Draw the axis
    ax=fig.add_axes(placement)
    #   Set line and marker sizes
    LW=1
    MW=10
    #   Set Ymax to zero for setting axis limits
    YMax=0.
    #   Now loop through all fits
        #   Get the line color -- this sets it to black
    linecol='k'
        #   Set the X and Y values
    X,Y,YErr=SetXY(Key,Model)
        #   Adjust the max value of Y
    if np.max(Y) > YMax:
        YMax=np.max(Y)

    nElem=len(Y)
    L1=ax.plot(X[:nElem],Y,marker='.',ls='-',color=linecol,lw=LW,markersize=MW)
    ax.errorbar(X[:nElem],Y,yerr=YErr,marker='.',ls='-',color=linecol,lw=LW,markersize=MW)

    ax.set_ylim(bottom=0.)
        #   Add minor tick marks to the panel.
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
    if Key=='VROT':
        ylabel =r"V (km/s)"
    elif Key=="SURFDENS_FACEON":
        ylabel = r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)"
    ax.set_xlabel(r"R ('')")
    ax.set_ylabel(ylabel)

    #KeywordPlotFmt(ax,Key,YMax,CatVal,DiagnosticSwitch)
    return ax


def SetXY(Key,FitParams):
    """
        This function sets the X and Y values for plotting
    """
    #   Start by seting to the bass X and Y
    X=FitParams['R']
    Y=FitParams[Key]
    ErrKey=Key+"_ERR"
    YErr=FitParams[ErrKey]
    #   If the key is the surface density, adjust to using the 'R_SD' for X
    if Key == 'SURFDENS':
        X=FitParams['R_SD']
    return X,Y,YErr
