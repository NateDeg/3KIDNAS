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




def MomentPlot(fig,placement,Moment,DataCube,MomMaps,Center,PixSize):
    ax=fig.add_axes(placement)

    #FileName=ObjDict['ObjFileBaseName']+"_mom"+str(Moment)+".fits.gz"
    #MomMap=fits.open(FileName)
    #size=MomMap[0].shape
    size=np.shape(DataCube['Data'])[1:3]
    print(size)
    
    MomMaps=CalcMomMapFromCube(DataCube,Moment,MomMaps)
    
    #   Set up the plot extent


    LowLim=np.zeros(2)
    HighLim=np.zeros(2)
    for i in range(2):
        CUse=Center[i]
        LowLim[i]=-CUse*PixSize
        HighLim[i]=(size[i]-CUse)*PixSize
        print("Moment map limits", LowLim[i],HighLim[i])
        
    extent=[LowLim[0],HighLim[0],LowLim[1],HighLim[1]]
    
    PlotName="Mom"+str(Moment)
    ax.text(0.5, 1.05,PlotName,fontsize=15,horizontalalignment='center',transform=ax.transAxes)
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.set_xlim(0,size[1]-1)
    ax.set_ylim(0,size[0]-1)
    #ax.set_xlim(LowLim[0],HighLim[0])
    #ax.set_ylim(LowLim[1],HighLim[1])
    minval=np.min(MomMaps)
    maxval=1.05*np.max(MomMaps)
    if Moment == 0:
        Cmap='viridis'
    elif Moment == 1:
        Cmap='jet_r'
        if DataCube['CubeVels'][0] < DataCube['CubeVels'][-1]:
            minval=DataCube['CubeVels'][0]
            maxval=DataCube['CubeVels'][-1]
        else:
            minval=DataCube['CubeVels'][-1]
            maxval=DataCube['CubeVels'][0]
        
    elif Moment ==2 :
        Cmap='plasma'
        #minval=np.nanmin(MomMap[0].data)
        #maxval=np.nanmax(MomMap[0].data)
    
    #ax.imshow(MomMaps[Moment],cmap=Cmap,vmin=minval,vmax=maxval,extent=extent,origin='lower')
    ax.imshow(MomMaps[Moment],cmap=Cmap,vmin=minval,vmax=maxval,origin='lower')

    if Moment ==0:
        Corner=-((Center[0]-1)-7)*PixSize

        Ell=Ellipse([Corner,Corner], 30., 30., 0.,edgecolor='cyan',facecolor='none',lw=2)
        #ax.add_patch(Ell)

    ax.set_xlabel(r"X (pixels)")
    ax.set_ylabel(r"Y (pixels)")
    #MomMapData=MomMap[0].data
    #MomMap.close()
    return ax,MomMaps,extent

def CalcMomMapFromCube(DataCube,Moment,MomMaps):
    YSize=np.shape(DataCube['Data'])[1]
    XSize=np.shape(DataCube['Data'])[2]
    TempMomMap=np.zeros([YSize,XSize])

    print("Initial shape",Moment,np.shape(MomMaps))
    FluxLim=1.e-7
    
    for i in range(YSize-1):
        for j in range(XSize-1):
            if Moment == 0:
                TempMomMap[i,j]=np.nansum(DataCube['MaskedData'][:,i,j])
            elif Moment ==1:
                TempMomMap[i,j]=np.nansum(DataCube['MaskedData'][:,i,j]*DataCube['CubeVels'][:])
    if Moment == 1:
        TempMomMap=TempMomMap/MomMaps[0]
        for i in range(YSize-1):
            for j in range(XSize-1):
                if MomMaps[0][i,j] < FluxLim:
                    TempMomMap[i,j]=TempMomMap[i,j]/0.
    
    MomMaps.append(TempMomMap)
    
    return MomMaps
    
    
def AddCenterToMomMap(ax,X,Y):
    #ax.plot(0,0,marker='x',color='k',markersize=2)
    ax.plot(X,Y,marker='x',color='k',markersize=2)
    return ax
    
def AddArrowToMomMap(ax,X,Y,Angle,Length):
    print("Arrow Length",Length)
    OffSetAng=90.
    EndX=X+Length*np.cos((Angle+OffSetAng)*np.pi/180.)
    EndY=Y+Length*np.sin((Angle+OffSetAng)*np.pi/180.)
    


    StartX=X-Length*np.cos((Angle+OffSetAng)*np.pi/180.)
    #StartY=Y-Length*np.sin((Angle+90.)*np.pi/180.)
    
    dX=Length*np.cos((Angle+OffSetAng)*np.pi/180.)
    dY=Length*np.sin((Angle+OffSetAng)*np.pi/180.)
    
    StartX=ax.get_xlim()[0]
    StartY=Y+(dY/dX)*(StartX-X)
    EndX=ax.get_xlim()[1]
    EndY=Y+(dY/dX)*(EndX-X)
    print("Start And End", StartX,EndX)

    ax.plot([StartX, EndX], [StartY, EndY],ls=':',color='k',marker='',)
    ax.arrow(X,Y,dX,dY,ls='-.',color='k',shape='full',width=1)
    return ax
    
def AddContoursToMomentPlot(ax,Moment,MomMap,VSys,ObsMap,Vels,Extent):
    #   First contruct the moment map
    #       Get the map size
    Size=(np.shape(MomMap)[0],np.shape(MomMap)[1])
    #       Initialize the moment map and flux tot arrays
    #MomMap=np.zeros(Size)
    #FTot=np.zeros(Size)
    #    Construct the moment map
    #for i in range(np.shape(Cube['Data'])[1]):
    #    for j in range(np.shape(Cube['Data'])[2]):
    #        for k in range(np.shape(Cube['Data'])[0]):
    #            FTot[i,j]+=Cube['Data'][k,i,j]
    #            MomMap[i,j]+=(Cube['CubeVels'][k]/1000.)**Moment*Cube['Data'][k,i,j]
    #       Adjust the Moment map so that it is nan'd where the observed map is zero
    for i in range(Size[0]):
        for j in range(Size[1]):
            if np.isnan(ObsMap[i,j]):
                MomMap[i,j]=MomMap[i,j]/0.

    #   Figue out the width of the cube
    #dV=(Cube['CubeVels'][-1]-Cube['CubeVels'][0])/1000.
    #dV=(np.nanmax(MomMap)-np.nanmin(MomMap))/1000.

    dV=(Vels[-1]-Vels[0])/1000.
    if dV < 0.:
        dV=-dV

    lTypes=('-')
    LW=1
    #   Set an array of a X and Y
    #X=np.linspace(0,np.shape(Cube['Data'])[1],np.shape(Cube['Data'])[1])
    #Y=np.linspace(0,np.shape(Cube['Data'])[2],np.shape(Cube['Data'])[2])
    X=np.linspace(0,Size[0],Size[0])
    Y=np.linspace(0,Size[1],Size[1])
    #X=np.linspace(Extent[0],Extent[1],Size[0])
    #Y=np.linspace(Extent[2],Extent[3],Size[1])
    YY2,XX2=np.meshgrid(Y,X)
    delV=dV/10.
    #   Set the contour levels
    CLevels=np.array([VSys-2.5*delV,VSys-2.0*delV,VSys-1.0*delV,VSys,VSys+1.0*delV,VSys+2.0*delV,VSys+2.5*delV])
    #   Draw on the contours
    #ax.contour(YY2,XX2,MomMap,colors='magenta',linewidths=LW,linestyles=lTypes,levels=CLevels,extent=Extent)
    ax.contour(YY2,XX2,MomMap,colors='magenta',linewidths=LW,linestyles=lTypes,levels=CLevels)
    #       Add in a thicker contour with vsys
    CLevels=np.array([VSys])
    ax.contour(YY2,XX2,MomMap,colors='magenta',linewidths=LW*2.5,linestyles=lTypes,levels=CLevels)
    
    return ax
    
    
def BasePVPlot(fig,placement,CubeInfo):
    ax=fig.add_axes(placement)
    
    #print("Getting pixel extent", np.shape(CubeInfo['PV'])[0],int(np.shape(CubeInfo['PV'])[0]/2))
    PixExtent=int(np.shape(CubeInfo['PV'])[0]/2)
    ArcSecExtent=[-PixExtent*CubeInfo['PixSize'],PixExtent*CubeInfo['PixSize']]
    #print("ArcSecExtent",ArcSecExtent)
    
    #X=np.linspace(1,np.shape(CubeInfo['PV'])[0],np.shape(CubeInfo['PV'])[0])
    X=np.linspace(ArcSecExtent[0],ArcSecExtent[1],np.shape(CubeInfo['PV'])[0])
    #print("X Check",X[PixExtent])
    
    V=CubeInfo['CubeVels']/1000.
    VV,XX=np.meshgrid(V,X)
    ax.pcolormesh(XX,VV,CubeInfo['PV'],cmap='Greys')
    #print(XX)
    #print(VV)
    #print(CubeInfo['PV'])
    
    ax.set_ylabel(r"V (km/s)" )
    #ax.set_xlabel(r"Major Axis (pixels)")
    
    return ax
    
def AddCentLinesToPVPlot(ax,VSys,XCenter):
    ax.axhline(y=VSys,ls=':',color='green',linewidth=2)
    ax.axvline(x=0,ls=':',color='green',linewidth=2)
    #print("XCenter for PV Plot", XCenter)
    return ax
    
def AddPVContoursToPlot(ax,CubeInfo):
    MaxPV=np.max(CubeInfo['PV'])
    CLevels=np.array([0.05,0.5])*MaxPV
    lTypes=(':','--','-')
    LW=1
    V2=CubeInfo['CubeVels']/1000.
    #X=np.linspace(1,np.shape(CubeInfo['PV'])[0],np.shape(CubeInfo['PV'])[0])
    PixExtent=int(np.shape(CubeInfo['PV'])[0]/2)
    ArcSecExtent=[-PixExtent*CubeInfo['PixSize'],PixExtent*CubeInfo['PixSize']]
    X=np.linspace(ArcSecExtent[0],ArcSecExtent[1],np.shape(CubeInfo['PV'])[0])
    VV2,XX2=np.meshgrid(V2,X)
    ax.contour(XX2,VV2,CubeInfo['PV'],levels=CLevels,colors='magenta',linewidths=LW,linestyles=lTypes)
    return ax
    
            
def AddSingleRC_to_PVPlot(ax,Model,CubeInfo):
    LW=1
    MW=10
    pixSize=CubeInfo['CubeHeader']['CDELT2']*3600.
        #   Get the line color
    linecol='blue'
        #   Get the line label
    linelabel=" "
        #   If there is a successful fit, start the plotting steps
    #X1=(Model['XCENTER']+1) +   Model['R']/pixSize
    X1=Model['R']
    V1=Model['VSYS']-(Model['VROT']*np.sin(Model['INCLINATION']*np.pi/180.))

    
    ax.plot(X1,V1,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)

    #X2=(Model['XCENTER']+1) -   Model['R']/pixSize
    X2= -Model['R']
    V2=Model['VSYS']+(Model['VROT']*np.sin(Model['INCLINATION']*np.pi/180.))
    ax.plot(X2,V2,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
        
 
    
def KeywordPlot_SingleModel(fig,placement,Key,Model,CubeMeasureDict):
    #   Draw the axis
    ax=fig.add_axes(placement)

    #   Set line and marker sizes
    LW=1
    MW=10


    #   Now loop through all fits
        #   Get the line color
    linecol='k'
        #   Get the line label
    linelabel="WRKP Fit"
        #   Set the X and Y values
    X,Y=SetXY(Key,Model)
        #   Plot the X-Y points
    ax.plot(X,Y,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
        #   For some key words, add in error bars
    AddErrorBars(ax,LW,MW,linelabel,linecol,Key,X,Y,Model)
    

                
    #   After the loop, add vertical lines for the radius
    ax.axvline(x=CubeMeasureDict['Radius'], color='k',ls='--',lw=LW)
    #ax.axvline(x=CubeMeasureDict['RHI'], color='k',ls='-.',lw=LW)
    #   Finally format the plot
    DiagnosticSwitch=1
    YMax=np.nanmax(Y)
    KeywordPlotFmt(ax,Key,YMax,DiagnosticSwitch)
    return ax

def SetXY(Key,FitParams):
    #   Start by seting to the bass X and Y
    X=FitParams['R']
    Y=FitParams[Key]
    #   If the key is the surface density, adjust to using the 'R_SD' for X
    if Key == 'SURFDENS':
        X=FitParams['R_SD']
    return X,Y
    
def AddErrorBars(ax,LW,MW,linelabel,linecol,Key,X,Y,FitParams):
    #   Only add error bars for a few keyword plots
    if Key == "VROT" or Key == "INCLINATION" or Key == "POSITIONANGLE":
        KeyErr=Key+"_ERR"
        YErr=FitParams[KeyErr]
        YErr=np.abs(YErr)
        ax.errorbar(X,Y,yerr=YErr,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
    
  
  
  
        
def ColorSelection(LineNum):
    if LineNum == -1:
        col='black'
    return col


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
    else:
        ylabel=" "
    return ylabel

def GetYLim(Ymax,Key,DiagnosticSwitch):
    #   This can be used for either normal diagnostic or model average plots with different surface density units
    YLim=Ymax
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
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
