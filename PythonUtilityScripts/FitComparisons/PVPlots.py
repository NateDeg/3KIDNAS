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


from . import ParameterPlots as PP

def MakePVPlots(CatID,WallCat,StrDict,ParamDicts,ObsDicts):
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
    
    
    print("Making PV Plots",np.shape(ParamDicts))
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    
    OutName=StrDict['BaseFileName']+"_"+NameSuffix+"_PVPlots.png"
    
    #       Diameter and radius calculation
    Diameter=WallCat.ell_maj[CatID]/(5.)
    Radius=Diameter/2.*30.

    WallCoords=[[WallCat.ra[CatID],WallCat.dec[CatID],0]]
    WallPix=ObsDicts['CubeWCS'].wcs_world2pix(WallCoords,0)
    WallCat.WallPix=WallPix
                
                

    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    #       Label the plot
    fig.text(.93,.95, Name , ha='center',rotation=0,va='center',size=27)

    nFit=np.shape(ParamDicts)[0]

    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15

    FitStep=0
    FitAdvance=0

    j=1
    k=2
    placement=[left+j*(w+buf),base+k*(h+buf),w,h]
    DataPVPlot(fig,placement,ObsDicts)


    j=0
    k=1
    for i in range(nFit):
        linelabel,FitStep,FitAdvance=PP.LineLabelSelection(i,FitStep,FitAdvance,StrDict)
        linecol=PP.ColorSelection(i)
        if ParamDicts[i]['FITAchieved']=='True':
            if FitAdvance==0:
                placement=[left+j*(w+buf),base+k*(h+buf),w,h]
                PVPanelPlot(fig,placement,WallCat,CatID,StrDict,ParamDicts[i],ObsDicts,linelabel,linecol)
                j=j+1
                if j == 3:
                    k=k-1
                    j=0


    plt.savefig(OutName, format='png',bbox_inches='tight')

    plt.close()

def DataPVPlot(fig,placement,ObsDicts):
    ax=fig.add_axes(placement)
    ax.set_ylabel(r"V")
    ax.set_xlabel(r"Major Axis")
    X=np.linspace(1,np.shape(ObsDicts['PV'])[0],np.shape(ObsDicts['PV'])[0])
    V=ObsDicts['CubeVels']/1000.
    VV,XX=np.meshgrid(V,X)
    ax.pcolormesh(XX,VV,ObsDicts['PV'],cmap='Greys')
    
    return ax,X,V


def PVPanelPlot(fig,placement,WallCat,CatID,StrDict,ParamDict,ObsDicts,linelabel,linecol):
    
    LW=3
    MW=10
    ax=fig.add_axes(placement)

    ax.set_ylabel(r"V")
    ax.set_xlabel(r"Major Axis")

    X=np.linspace(1,np.shape(ObsDicts['PV'])[0],np.shape(ObsDicts['PV'])[0])
    V=ObsDicts['CubeVels']/1000.
    VV,XX=np.meshgrid(V,X)
    ax.pcolormesh(XX,VV,ObsDicts['PV'],cmap='Greys')

    ax.text(0.5, 1.05, linelabel, color=linecol,size=18,ha='center',transform=ax.transAxes)

    MaxPV=np.max(ObsDicts['PV'])
    CLevels=np.array([0.05,0.3])*MaxPV
    lTypes=(':','--','-')
    
    #print(MaxPV,np.max(ParamDict['CUBE']['PV']))
    
    V2=ParamDict['CUBE']['CubeVels']
    VV2,XX2=np.meshgrid(V2,X)
    
    ax.contour(XX2,VV2,ParamDict['CUBE']['PV'],levels=CLevels,colors=linecol,linewidths=LW,linestyles=lTypes)
#ax.pcolormesh(XX2,VV2,ParamDict['CUBE']['PV'],cmap='jet')

