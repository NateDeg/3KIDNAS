import numpy as np
import pandas as pd


import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator


def BootstrapBoxPlot(GalaxyDict,GeneralDict,BootstrapParams):

    fig, PltOpts=PlotIni()
    #   We will need to the best fitting model
    
    
    BestFit=GalaxyDict['BestFitModel']


    PltKeys=['INCLINATION','POSITIONANGLE','VSYS','XCENTER','YCENTER','RHI','VHI']

    Labels={'INCLINATION':r"Inc ($^{\circ}$)",'POSITIONANGLE':r"PA ($^{\circ}$)",'VSYS':r"V$_{\rm{sys}}$ (km s$^{-1}$)",'RHI':r"R$_{\rm{HI}}$ ('')",'VHI':r"V$_{\rm{HI}}$ (km s$^{-1}$)",'XCENTER':r"X (pix)",'YCENTER':r"Y (pix)"}
    #Lims={'INCLINATION':[0.,90.],'POSITIONANGLE':[190.,230.],'VSYS':[820,832],'RHI':[1.1,2.1],'VHI':[8.,25.],'XCENTER':r"X (pix)",'YCENTER':r"Y (pix)"}
    #Ticks={'INCLINATION':[40.,50,60.],'POSITIONANGLE':[200.,220.],'VSYS':[825,830],'RHI':[1.5,2.0],'VHI':[10.,20.]}
    
    fig,PltOpts=PlotIni()
    
    nRow=len(PltKeys)
    nCol=nRow

    for i in range(nCol):
        xKey=PltKeys[i]
        for j in range(nCol-i):
            yKey=PltKeys[nCol-j-1]
            print(i,j,xKey,yKey)
            placement=[PltOpts['left']+i*PltOpts['w'],PltOpts['base']+j*PltOpts['h'],PltOpts['w'],PltOpts['h']]
            ax=fig.add_axes(placement)
            
            #ax.set_xlim(Lims[xKey])
            #ax.set_ylim(Lims[yKey])

            xBFVal=BestFit[xKey]
            yBFVal=BestFit[yKey]
            if xKey not in ['RHI','VHI']:
                xBFVal=xBFVal[0]
            if yKey not in ['RHI','VHI']:
                yBFVal=yBFVal[0]
              
            if xKey==yKey:
                
                HistPlt(ax,BootstrapParams,BestFit,xKey,yBFVal)

            else:
                for k in range(len(BootstrapParams)):
                    if xKey=='VHI' or xKey=='RHI':
                        XUse=BootstrapParams[k][xKey]
                    else:
                        XUse=BootstrapParams[k][xKey][0]
                    if yKey=='VHI' or yKey=='RHI':
                        YUse=BootstrapParams[k][yKey]
                    else:
                        YUse=BootstrapParams[k][yKey][0]
                    ax.scatter(XUse,YUse,marker='o',color='#f44336',s=100,alpha=0.5)
                ax.scatter(xBFVal,yBFVal,s=2000,marker='*',color='#008da9')

            PanelFmt(ax,xKey,yKey,i,j,Labels)

    
    PltName=GalaxyDict['WRKP_ResultsFolder']+GalaxyDict['ObjName']+"_GeoBoxPlot.png"

    #   Save the plot
    plt.savefig(PltName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()

def HistPlt(ax,BootstrapParams,BestFit,xKey,yBFVal):

    BUse=np.empty(len(BootstrapParams),dtype='object')
    for k in range(len(BootstrapParams)):
        if xKey=='RHI' or xKey=='VHI':
            BUse[k]=BootstrapParams[k][xKey]
        else:
            BUse[k]=BootstrapParams[k][xKey][0]
    ax.hist(BUse,bins=20, orientation='horizontal',color='#f44336',histtype='step',lw=5)

    ax.axhline(y=yBFVal,color='#008da9',lw=10)
    
def PanelFmt(ax,xKey,yKey,col,row,Labels):
    
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    
    if xKey==yKey:
        HistFmt(ax,xKey,row,col,Labels)
        
    else:
    
        if row==0:
            ax.set_xlabel(Labels[xKey])
        else:
          ax.xaxis.set_major_formatter(NullFormatter())

        if col==0:
            ax.set_ylabel(Labels[yKey])
        else:
          ax.yaxis.set_major_formatter(NullFormatter())
        
def HistFmt(ax,xKey,row,col,Labels,Lims,Ticks):
    if row==0:
        ax.set_xlabel('N')
    else:
       ax.xaxis.set_major_formatter(NullFormatter())


    if col==0:
        ax.set_ylabel(Labels[xKey])
    else:
       ax.yaxis.set_major_formatter(NullFormatter())



def PlotIni():
    BasePlotParams={'font.size': 15,'axes.linewidth':2
        ,'xtick.major.size':6,'xtick.minor.size':3
        ,'xtick.major.width':1,'xtick.minor.width':1
            ,'ytick.major.size':6,'ytick.minor.size':3
            ,'ytick.major.pad':10,'xtick.major.pad':10
            ,'ytick.major.width':1,'ytick.minor.width':1
            ,'xtick.labelsize':25 ,'ytick.labelsize':25
            ,'axes.labelsize': 40
            ,'legend.fontsize': 18
                }
    matplotlib.rcParams.update(BasePlotParams)
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    
    base=-0.1
    left=0.1
    w=0.45
    h=w
    buf=0.15
    cW=0.02
    hbuf=0.2
    PltOpts={'base':base,'left':left,'w':w,'h':h,'buf':buf,'cW':cW,'hbuf':hbuf}
    return fig, PltOpts
