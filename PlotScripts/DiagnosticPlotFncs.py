import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from matplotlib.patches import Ellipse


def PlotFnc():
    PltFncDict={'ColorSelection':ColorSelection,'MomentPlot':MomentPlot,'AddCenterToMomMap':AddCenterToMomMap,'AddArrowToMomMap':AddArrowToMomMap,'AddContoursToMomentPlot':AddContoursToMomentPlot,'BasePVPlot':BasePVPlot,'AddRCs_to_PVPlot':AddRCs_to_PVPlot,'AddPVContoursToPlot':AddPVContoursToPlot,'AddSingleRC_to_PVPlot':AddSingleRC_to_PVPlot,'AddCentLinesToPVPlot':AddCentLinesToPVPlot,'KeywordPlot_SingleModel':KeywordPlot_SingleModel}
    return PltFncDict



      
 
    
def KeywordPlot(fig,placement,Key,ModelDict):
    #   Draw the axis
    ax=fig.add_axes(placement)

    #   Set line and marker sizes
    LW=1
    MW=10


    #   Now loop through all fits
        #   Get the line color
    YMax=0.
    for i in range(ModelDict['nModel']):
        #linecol=ColorSelection(i)
        linecol='black'
        AL=1.
        if i >=1:
            linecol='blue'
            AL=0.3
        #   Get the line label
        linelabel=ModelDict['Labels'][i]
        #print(linelabel)
        #   Set the X and Y values
        for j in range(3):
            ls=LineStyleSelection(j)
            if j ==0:
                ModelUse=ModelDict['IniModel'][i]
            elif j ==1 :
                ModelUse=ModelDict['IniModel'][i]
            elif j ==2:
                ModelUse=ModelDict['FinalModel'][i]
            X,Y=SetXY(Key,ModelUse)
        #   Plot the X-Y points
            ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW,label=linelabel,alpha=AL)
        
            YMax=np.max([np.nanmax(Y),YMax])
            
    #   Finally format the plot
    DiagnosticSwitch=1
    KeywordPlotFmt(ax,Key,YMax,DiagnosticSwitch)
    return ax
    
    
def KeywordPlot_Ensamble(fig,placement,Key,ModelDict,EnsambleModels):
    #   Draw the axis
    ax=fig.add_axes(placement)

    #   Set line and marker sizes
    LW=1
    MW=10
    YMax=0.
    for j in range(len(EnsambleModels)):
        for i in range(len(EnsambleModels[j])):
            ModelUse=EnsambleModels[j][i]
            X,Y=SetXY(Key,ModelUse)
            
            ls='-'
            if j==0:
                linecol='blue'
            elif j==1:
                linecol='green'
            ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW,alpha=0.2)
            YMax=np.max([np.nanmax(Y),YMax])
        #print(i,X,Y)
    #   Now loop through all fits
        #   Get the line color
    for i in range(ModelDict['nModel']):
        linecol=ColorSelection(i)
        #   Get the line label
        linelabel=ModelDict['Labels'][i]
        #print(linelabel)
        #   Set the X and Y values
        for j in range(3):
            ls=LineStyleSelection(j)
            if j ==0:
                ModelUse=ModelDict['IniModel'][i]
            elif j ==1 :
                ModelUse=ModelDict['IniModel'][i]
            elif j ==2:
                ModelUse=ModelDict['FinalModel'][i]
            X,Y=SetXY(Key,ModelUse)
        #   Plot the X-Y points
            ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW,label=linelabel )
        
            YMax=np.max([np.nanmax(Y),YMax])
            
    #   Finally format the plot
    DiagnosticSwitch=1
    KeywordPlotFmt(ax,Key,YMax,DiagnosticSwitch)
    return ax
    
    
def KeywordPlot_CodeComp(fig,placement,Key,ModelDict):
    #   Draw the axis
    ax=fig.add_axes(placement)

    #   Set line and marker sizes
    LW=1
    MW=10
    YMax=0.
    
    ModelUse=ModelDict['TrueModel']
    linecol='black'
    ls='-'
    X,Y=SetXY(Key,ModelUse)
    YMax=np.max([np.nanmax(Y),YMax])
    ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW )
    
    #print(np.shape(ModelDict['Models']))
    nModels=np.shape(ModelDict['Models'])[0]
    nCodes=np.shape(ModelDict['Models'])[1]
    #print(nModels,nCodes)
    
    AL=0.2
    for i in range(nCodes):
        for j in range(nModels):
            ModelUse=ModelDict['Models'][j,i]
            if i ==0:
                linecol='blue'
            elif i ==1:
                linecol='red'
            elif i ==2:
                linecol='green'
            X,Y=SetXY(Key,ModelUse)
            YMax=np.max([np.nanmax(Y),YMax])
            ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW,alpha=AL)
            if Key=='SURFDENS_FACEON':
                KeySpec='SURFDENS_FACEON2'
                X,Y=SetXY(KeySpec,ModelUse)
                YMax=np.max([np.nanmax(Y),YMax])
                ax.plot(X,Y,marker='.',ls=ls,color=linecol,lw=LW,markersize=MW,alpha=AL)
                
    #   Finally format the plot
    DiagnosticSwitch=1
    KeywordPlotFmt(ax,Key,YMax,DiagnosticSwitch)
    return ax

def SetXY(Key,FitParams):
    #   Start by seting to the bass X and Y
    #print("Fit Params",FitParams)
    X=FitParams['R']
    Y=FitParams[Key]
    #   If the key is the surface density, adjust to using the 'R_SD' for X
    if Key == 'SURFDENS' or Key == 'SURFDENS_FACEON' or Key == 'SURFDENS_FACEON2':
        X=FitParams['R_SD']
    return X,Y

    
  
def ColorSelection(LineNum):
    if LineNum == 0:
        col='black'
    elif LineNum==1:
        col='red'
    elif LineNum==2:
        col='blue'
    elif LineNum==3:
        col='magenta'
    elif LineNum==4:
        col='cyan'

    return col
    
def LineStyleSelection(ModelType):
    if ModelType == 0:
        ls=':'
    elif ModelType ==1:
        ls='-.'
    elif ModelType ==2:
        ls='-'
    return ls


def LabelSelection(Key,DiagnosticSwitch):
    if Key == "VROT":
        #   This can be used for either normal diagnostic or single model plots with different labels for VROT and SURFDENS
        if DiagnosticSwitch == 0:
            ylabel=r"$V_{\phi}$"
        elif DiagnosticSwitch== 1:
            ylabel =r"VROT_kin (km/s)"
    
    elif Key == "SURFDENS":
        if DiagnosticSwitch == 0:
            ylabel=r"$\Sigma$ "
        elif DiagnosticSwitch== 1:
            ylabel = r"SD_kin (M$_{\odot}$ pc$^{-2}$)"
            
    elif Key=="INCLINATION":
        ylabel=r"i ($^\circ$)"
        
    elif Key=="POSITIONANGLE":
        ylabel=r"PA ($^\circ$)"
        
    elif Key=="XCENTER":
        ylabel=r"X (pixels)"
        
    elif Key=="YCENTER":
        ylabel=r"Y (pixels)"
        
    elif Key=="VSYS":
        ylabel=r"VSYS (km/s)"
        
        
    else:
        ylabel=" "
    return ylabel

def GetYLim(Ymax,Key,DiagnosticSwitch):
    #   This can be used for either normal diagnostic or model average plots with different surface density units
    YLim=Ymax*1.005
    if Key == "VROT":
        YLim=np.min((1.1*YLim,400.))
    elif Key == "SURFDENS":
        if DiagnosticSwitch== 0:
            YLim=np.min((1.1*YLim,1.e-3))
        elif DiagnosticSwitch== 1:
            YLim=np.min((1.1*YLim,5.e1))
    return YLim




def CatLine(Key,ObjDict,AstroFncs,CubeMeasureDict):
    drawCatVal=False
    CatVal=[]
    if Key == "INCLINATION":
        drawCatVal=True
        CatVal=np.arccos(ObjDict['CatEntry']['ell_min']/ObjDict['CatEntry']['ell_maj'])*180./np.pi
    elif Key == "POSITIONANGLE":
        drawCatVal=True
        CatVal=ObjDict['CatEntry']['kin_pa']+180.
        if CatVal > 360.:
            CatVal=CatVal-360.
    elif Key == "VSYS":
       drawCatVal=True
       Freq=ObjDict['CatEntry']['freq']
       RestFreq=1.42040575179E+09
       CatVal=AstroFncs['RedShiftConv'](Freq,RestFreq)/1000.
       
    elif Key == "XCENTER":
        drawCatVal=True
        CatVal=CubeMeasureDict['CentPix'][0]
    elif Key == "YCENTER":
        drawCatVal=True
        CatVal=CubeMeasureDict['CentPix'][1]


    return drawCatVal,CatVal


def KeywordPlotFmt(ax,Key,YMax,DiagnosticSwitch):
    #   Label the plot
    ax.set_xlabel(r"R ('')")
    ylabel=LabelSelection(Key,DiagnosticSwitch)
    ax.set_ylabel(ylabel)
    #   Set the Ymax
    YMax=GetYLim(YMax,Key,DiagnosticSwitch)
    ax.set_ylim(top=YMax)
    if Key == "VROT" or Key == "SURFDENS":
        ax.set_ylim(bottom=0.)
    
    if Key == "SURFDENS" or Key == 'SURFDENS_FACEON':
        ApproxNoise=1.6
        ApproxNoise/=1000.
        BeamArea=np.pi*(5.)**2.
        PixelArea=6.**2.
        VelWidth=3.93346533
        
        SDConv1=3.93346533/36.
        SDConv2=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
        SDConv3=BeamArea
        ApproxNoise=ApproxNoise/SDConv3*VelWidth/SDConv2
        print("Approximate SD Noise", ApproxNoise)
        ax.axhline(y=ApproxNoise,ls='--',color='k')
        ax.set_xlim(0,200)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
