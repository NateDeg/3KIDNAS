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



def MakeAllPlots(CatID,WallCat,StrDict,ParamDicts,ObsDicts,BootstrapDicts):
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
    
    
    print("Making Plots",np.shape(ParamDicts))
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    
    OutName=StrDict['BaseFileName']+"_"+NameSuffix+"_Bootstrap_FitComparisons.png"
    
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
    fig.text(.6,.95, Name , ha='center',rotation=0,va='center',size=27)

    
    FitStep=0
    FitAdvance=0
    nFit=np.shape(ParamDicts)[0]
    yTextLoc=0.85

#LabelStr="Integrated HI\t =\t".expandtabs()+str(round(ObsDicts['IntegratedHI'],3) )+ r" Jy"
#fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    yTextLoc-=0.1


    LabelStr="Ell_Maj \t=\t".expandtabs()+str(round(Diameter,3))+" beams"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    yTextLoc-=0.1

    MStr='%2E'%Decimal(ObsDicts['Mass'])
    LabelStr="Mass\t\t= \t".expandtabs()+MStr+r" M$_{\odot}$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    yTextLoc-=0.1
    
    LabelStr="Distance \t=\t".expandtabs()+str(round(ObsDicts['Distance'],3) )+" Mpc"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    yTextLoc-=0.15

    

    fig.text(1.45,yTextLoc, r"Goodness of Fit" , ha='center',rotation=0,va='center',size=22)
    yTextLoc-=0.1
    FitAdvance=0

    for i in range(nFit-1):
        
        Label,FitStep,FitAdvance=LineLabelSelection(i,FitStep,FitAdvance,StrDict)
        LabelColor=ColorSelection(i)
        #print(i, Label, LabelColor,ParamDicts[i]['FITAchieved'])
        if ParamDicts[i]['FITAchieved']== 'True':
            chi2=ParamDicts[i]['CHI2']
            LabelStr=Label+"\t = \t".expandtabs()+str(round(chi2,3) )
            if FitAdvance ==0:
                fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=18,color=LabelColor)
                yTextLoc-=0.05

    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    
    placement=[left-0.8*(w+buf),base+2*(h+buf),w*0.8,h]
    PV_RC_Plot(fig,placement,WallCat,CatID,StrDict,ParamDicts,ObsDicts)
    
    placement=[left-0.8*(w+buf),base+1*(h+buf),w,h]
    Moment=0
    MomentPlot(fig,Moment,placement,WallCat,CatID,StrDict)
    
    placement=[left-0.8*(w+buf),base+0.*(h+buf),w,h]
    Moment=1
    MomentPlot(fig,Moment,placement,WallCat,CatID,StrDict)
    
    placement=[left-0.8*(w+buf),base-1.*(h+buf),w,h]
    Moment=2
    MomentPlot(fig,Moment,placement,WallCat,CatID,StrDict)
    
    

    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)

    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left,base+1*(h+buf),w,h]
    Key="INCLINATION"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left+(w+buf),base+1*(h+buf),w,h]
    Key="POSITIONANGLE"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left,base+0*(h+buf),w,h]
    Key="VSYS"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left+(w+buf),base+0*(h+buf),w,h]
    Key="VDISPERSION"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left,base-1*(h+buf),w,h]
    Key="XCENTER"
    KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    
    placement=[left+(w+buf),base-1*(h+buf),w,h]
    Key="YCENTER"
    ax1=KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels,bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    
    

    plt.savefig(OutName, format='png',bbox_inches='tight')

    plt.close()

def KeywordPlot(fig,Key,placement,ParamDicts,StrDict,WallCat,Radius,CatID,BootstrapDicts):
    LW=1
    MW=10
    nFit=np.shape(ParamDicts)[0]
    
    #print(Key)

    ax=fig.add_axes(placement)
    FitStep=0
    FitAdvance=0
    drawCatalogueLine,LineVal=CatLine(Key,WallCat,CatID)
    if drawCatalogueLine == 'True':
        ax.axhline(y=LineVal,color='k',ls='--',lw=LW)
    
    for i in range(nFit):
        linecol=ColorSelection(i)
        linelabel,FitStep,FitAdvance=LineLabelSelection(i,FitStep,FitAdvance,StrDict)
        if ParamDicts[i]['FITAchieved'] == 'True':
            #print(i, linelabel,linecol,FitStep,FitAdvance)
            X=ParamDicts[i]['R']
            Y=ParamDicts[i][Key]
            if Key == "SURFDENS":
                X=ParamDicts[i]['R_SD']
            ax.plot(X,Y,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
    i=nFit
    linecol=ColorSelection(i)
    linelabel,FitStep,FitAdvance=LineLabelSelection(i,FitStep,FitAdvance,StrDict)
    X=BootstrapDicts['R']
    if Key == "SURFDENS":
        X=BootstrapDicts['R_SD']
    Y=BootstrapDicts[Key][0]
    YErr=BootstrapDicts[Key][1]
    ax.errorbar(X,Y,yerr=YErr,color=linecol,lw=LW,markersize=MW,label=linelabel)

    ax.set_xlabel(r"R ('')")
    ylabel=LabelSelection(Key)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.axvline(x=Radius, color='k',ls='--',lw=LW)
    return ax

def ColorSelection(LineNum):
    if LineNum == 0:
        col='red'
    elif LineNum == 1:
        col='blue'
    elif LineNum == 2:
        col='green'
    elif LineNum == 3:
        col='magenta'
    elif LineNum == 4:
        col='purple'
    elif LineNum == 5:
        col='orange'
    elif LineNum == 6:
        col='brown'
    elif LineNum == 7:
        col='deeppink'
    elif LineNum == 8:
        col='olive'
    elif LineNum == 9:
        col='royalblue'
    elif LineNum == 10:
        col='indigo'
    elif LineNum == 11:
        col='salmon'
    elif LineNum == 12:
        col='teal'
    elif LineNum == 13:
        col='forestgreen'
    
    
    return col

def LineLabelSelection(ParamStep,FitStep,FitAdvance,StrDict):
    if StrDict['FitPaths'][FitStep][1]==0:
        linelabel=StrDict['FitPaths'][FitStep][3]
        FitAdvance =0
    elif StrDict['FitPaths'][FitStep][1]==1:
        if FitAdvance == 0:
            linelabel=StrDict['FitPaths'][FitStep][2]
            FitAdvance=1
        elif FitAdvance ==1:
            linelabel=StrDict['FitPaths'][FitStep][3]
            FitAdvance=0
    if FitAdvance == 0:
        FitStep+=1

    return linelabel,FitStep,FitAdvance

def LabelSelection(Key):
    if Key == "VROT":
        ylabel=r"$V_{\phi}$"
    elif Key == "SURFDENS":
        ylabel=r"$\Sigma$ "
    elif Key == "INCLINATION":
        ylabel=r"$I$ "
    elif Key == "POSITIONANGLE":
        ylabel=r"$PA$ "
    elif Key == "VSYS":
        ylabel=r"$V_{sys}$ "
    elif Key == "VDISPERSION":
        ylabel=r"$\sigma_{r}$ "
    elif Key == "XCENTER":
        ylabel=r"$X$ "
    elif Key == "YCENTER":
        ylabel=r"$Y$ "
    else:
        ylabel=" "
    return ylabel

def CatLine(Key,WallCat,CatID):
    drawCatVal='False'
    CatVal=[]
    if Key == "INCLINATION":
        drawCatVal='True'
        CatVal=np.arccos(WallCat.ell_min[CatID]/WallCat.ell_maj[CatID])*180./np.pi
        print(WallCat.ell_min[CatID],WallCat.ell_maj[CatID],CatVal)
    elif Key == "POSITIONANGLE":
        drawCatVal='True'
        CatVal=WallCat.kin_pa[CatID]+180.
        if CatVal > 360.:
            CatVal=CatVal-360.
    elif Key == "VSYS":
       drawCatVal='True'
       CatVal=WallCat.cz[CatID]
    elif Key == "XCENTER":
        drawCatVal='True'
        CatVal=WallCat.WallPix[0,0]
    elif Key == "YCENTER":
        drawCatVal='True'
        CatVal=WallCat.WallPix[0,1]


    return drawCatVal,CatVal

def MomentPlot(fig,Moment,placement,WallCat,CatID,StrDict):
    ax=fig.add_axes(placement)
    
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    FileName=StrDict['DataFolder']+StrDict['BaseFileName']+"_source_products/"+StrDict['BaseFileName']+"_"+NameSuffix+"_moment"+str(Moment)+".fits"
    MomMap=fits.open(FileName)
    size=MomMap[0].shape
    
    PlotName="Mom"+str(Moment)
    ax.text(0.5, 1.05,PlotName,fontsize=20,horizontalalignment='center',transform=ax.transAxes)
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.set_xlim(0,size[1]-1)
    ax.set_ylim(0,size[0]-1)
    minval=np.min(MomMap[0].data)
    maxval=1.05*np.max(MomMap[0].data)
    if Moment == 0:
        Cmap='viridis'
    elif Moment == 1:
        Cmap='jet_r'
        centfreq=WallCat.freq[CatID]
        w50=WallCat.w50[CatID]
        minval=centfreq-w50
        maxval=centfreq+w50
    elif Moment ==2 :
        Cmap='plasma'
        minval=np.nanmin(MomMap[0].data)
        maxval=np.nanmax(MomMap[0].data)
    ax.imshow(MomMap[0].data,cmap=Cmap,vmin=minval,vmax=maxval)
    if Moment ==0:
        Ell=Ellipse([7,7], 5., 5., 0.,edgecolor='cyan',facecolor='none',lw=2)
        ax.add_patch(Ell)

    MomMap.close()


def PV_RC_Plot(fig,placement,WallCat,CatID,StrDict,ParamDicts,ObsDicts):
    ax=fig.add_axes(placement)
    LW=1
    MW=10


    X=np.linspace(1,np.shape(ObsDicts['PV'])[0],np.shape(ObsDicts['PV'])[0])
    V=ObsDicts['CubeVels']/1000.
    VV,XX=np.meshgrid(V,X)
    ax.pcolormesh(XX,VV,ObsDicts['PV'],cmap='Greys')
    ax.set_ylabel(r"V")
    ax.set_xlabel(r"Major Axis")
    pixSize=abs(ObsDicts['CubeHeader']['CDELT1'])*3600.
    FitStep=0
    FitAdvance=0
    
    
    
    nFit=np.shape(ParamDicts)[0]
    for i in range(nFit):
        if ParamDicts[i]['FITAchieved'] == 'True':
            linecol=ColorSelection(i)
            linelabel,FitStep,FitAdvance=LineLabelSelection(i,FitStep,FitAdvance,StrDict)
            
            
            X1=ParamDicts[i]['XCENTER'] +   ParamDicts[i]['R']/pixSize
            V1=ParamDicts[i]['VSYS']-(ParamDicts[i]['VROT']*np.sin(ParamDicts[i]['INCLINATION']*np.pi/180.))
            ax.plot(X1,V1,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)

            X2=ParamDicts[i]['XCENTER'] -   ParamDicts[i]['R']/pixSize
            V2=ParamDicts[i]['VSYS']+(ParamDicts[i]['VROT']*np.sin(ParamDicts[i]['INCLINATION']*np.pi/180.))
            ax.plot(X2,V2,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)



