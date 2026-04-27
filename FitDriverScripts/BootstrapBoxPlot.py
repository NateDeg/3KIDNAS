"""
BootstrapBoxPlot.py
=======

**Author:** Nathan Deg

**Description:**
This module contains functions used to produce a diagnostic plot showing the individual bootstrap fits used to construct the uncertainty.
"""

import numpy as np
import pandas as pd


import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator


def BootstrapBoxPlot(GalaxyDict,GeneralDict,BootstrapParams):
    """
    This function is the main function for building the box plot.
    Parameters
            ----------
            GalaxyDict : dictionary
                The dictionary containing the galaxy measurements
            GeneralDict : dictionary
                The dicationary containing the general definitions
            BoostrapParams : array
                An array of dictionaries containing all bootstrap fits
    """
    #   Initialize the plot
    fig, PltOpts=PlotIni()
    #   Grab the best fit model
    BestFit=GalaxyDict['BestFitModel']
    #   Set the keys for plotting scalar parameters
    PltKeys=['INCLINATION','POSITIONANGLE','VSYS','XCENTER','YCENTER','RHI','VHI']
    #   Set the corresponding labels
    Labels={'INCLINATION':r"Inc ($^{\circ}$)",'POSITIONANGLE':r"PA ($^{\circ}$)",'VSYS':r"V$_{\rm{sys}}$ (km s$^{-1}$)",'RHI':r"R$_{\rm{HI}}$ ('')",'VHI':r"V$_{\rm{HI}}$ (km s$^{-1}$)",'XCENTER':r"X (pix)",'YCENTER':r"Y (pix)"}

    #   Now set the number of rows and columns needed for the triangle box plots
    nRow=len(PltKeys)
    nCol=nRow
    #   Loop through each combination of scalar keys
    for i in range(nCol):
        #   Select the x-key
        xKey=PltKeys[i]
        #   Avoid double counting and build the triangle by going from nCol-i
        for j in range(nCol-i):
            #   Set the matched yKey
            yKey=PltKeys[nCol-j-1]
            #   Set the rectangle for the current panel
            placement=[PltOpts['left']+i*PltOpts['w'],PltOpts['base']+j*PltOpts['h'],PltOpts['w'],PltOpts['h']]
            #   Add the axis
            ax=fig.add_axes(placement)
            #   Get the best fit values
            try:
                xBFVal=BestFit[xKey]
                yBFVal=BestFit[yKey]
            except:
                xBFVal=[np.nan]
                yBFVal=[np.nan]
                
            try:
                yBFVal=BestFit[yKey]
            except:
                yBFVal=[np.nan]
            #   RHI and VHI are stored as one number, while the others are in arrays, so select the first element
            if xKey not in ['RHI','VHI']:
                xBFVal=xBFVal[0]
            if yKey not in ['RHI','VHI']:
                yBFVal=yBFVal[0]
            #   When x==y, add a histogram plot
            if xKey==yKey:
                HistPlt(ax,BootstrapParams,BestFit,xKey,yBFVal)
            #   Otherwise make a scatter plot
            else:
                # Loop through all bootstrap fits and plot them
                for k in range(len(BootstrapParams)):
                    #   Again, the fact taht VHI and RHI aren't in arrays, means the selection is a little different
                    if BootstrapParams[k]['FITAchieved']==False:
                        continue
                    if xKey=='VHI' or xKey=='RHI':
                        XUse=BootstrapParams[k][xKey]
                    else:
                        XUse=BootstrapParams[k][xKey][0]
                    if yKey=='VHI' or yKey=='RHI':
                        YUse=BootstrapParams[k][yKey]
                    else:
                        YUse=BootstrapParams[k][yKey][0]
                    #   Scatter the points
                    ax.scatter(XUse,YUse,marker='o',color='#f44336',s=100,alpha=0.5)
                #   Add the best fit value
                ax.scatter(xBFVal,yBFVal,s=2000,marker='*',color='#008da9')
            #   Format the panel
            PanelFmt(ax,xKey,yKey,i,j,Labels)
    #   Now we'll add in the rotation curve and surface density profiles
    #   Again, start by setting the keys
    PKeys=['VROT','SURFDENS_FACEON']
    XKeys={'VROT':'R','SURFDENS_FACEON':'R_SD'}
    EKeys={'VROT':'VROT_ERR','SURFDENS_FACEON':'SURFDENS_FACEON_ERR'}
    Labels={'VROT':r"V (km s$^{-1}$)",'SURFDENS_FACEON':r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)"}
    #   Loop through the RC and SD profiles
    for i in range(2):
        #   Select the appropriate X, Y, and YErr keys
        Key=PKeys[i]
        XKey=XKeys[Key]
        EKey=EKeys[Key]
        #   Set the placement of the panels
        placement=[PltOpts['left']+4.25*PltOpts['w'],PltOpts['base']+(5.25-i)*PltOpts['h']+(1-i)*0.25*PltOpts['h'],3.0*PltOpts['w'],1.5*PltOpts['h']]
        #   Add the panel
        ax=fig.add_axes(placement)
        #   Label and format the panels
        ax.set_ylabel(Labels[Key])
        ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
        
        if i==1:
            ax.set_xlabel(r"R ('')")
        else:
            ax.xaxis.set_major_formatter(NullFormatter())
        #   Plot all bootstrap profiles
        for j in range(len(BootstrapParams)):
            if BootstrapParams[j]['FITAchieved']:
                ax.plot(BootstrapParams[j][XKey],BootstrapParams[j][Key],color='#f44336',lw=3,alpha=0.1)
        #   Add the best fitting model
        ax.errorbar(BestFit[XKey],BestFit[Key],yerr=np.abs(BestFit[EKey]),color='#008da9',lw=5,ls='-',alpha=0.9)
        
    #   Set the plot name
    PltName=GalaxyDict['WRKP_ResultsFolder']+GalaxyDict['ObjName']+"_GeoBoxPlot.png"
    #   Save the plot
    plt.savefig(PltName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()

def HistPlt(ax,BootstrapParams,BestFit,xKey,yBFVal):
    """
    This function make the 1D histogram plots.
    Parameters
            ----------
            ax : matplotlib axes
                The panel in the figure to be plotted
            BoostrapParams : array
                An array of dictionaries containing all bootstrap fits
            BestFit : dictionary
                The best fit model dictionary
            xKey : string 
                The specific parameter to plot
            yBFVal : float
                The best fit value
    """
    #   The first step is to make an array of the bootstrap values
    #       Initialize the empty array
    BUse=np.empty(len(BootstrapParams),dtype='object')
    #   Loop through all bootstraps
    for k in range(len(BootstrapParams)):
        if BootstrapParams[k]['FITAchieved']==False:
            BUse[k]=np.nan
            continue
        #   Once again, there is the issue about VHI and RHI shapes
        if xKey=='RHI' or xKey=='VHI':
            BUse[k]=BootstrapParams[k][xKey]
        else:
            BUse[k]=BootstrapParams[k][xKey][0]
    #   With the array constructed, make the histogram.  Generally these are horizontal rather than vertical.
    try:
        ax.hist(BUse,bins=20, orientation='horizontal',color='#f44336',histtype='step',lw=5)
    except:
        pass
    #   Add a horizontal line at the best fit value
    ax.axhline(y=yBFVal,color='#008da9',lw=10)
    
def PanelFmt(ax,xKey,yKey,col,row,Labels):
    """
    This function formats a panel
    Parameters
            ----------
            ax : matplotlib axes
                The panel in the figure to be plotted
            xKey : string 
                The x-parameter
            yKey : string 
                The y-parameter
            row : integer
                The current row of the panel
            col : integer 
                The current column of the panel
            Labels : dictionary
                A dictionary of panel labels
    """
    #   Start by setting the minor axis ticks
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    #   if the x and y keys are the same, this is a histogram which has a slightly different formating
    if xKey==yKey:
        HistFmt(ax,xKey,row,col,Labels)
    else:
        #   Only add x labels on the bottom row
        if row==0:
            ax.set_xlabel(Labels[xKey])
        else:
            #   Otherwise remove the major ticks
            ax.xaxis.set_major_formatter(NullFormatter())
            #   Only add y labels on the leftmost column
        if col==0:
            ax.set_ylabel(Labels[yKey])
        else:
            #   Otherwise remove the major tick marks
          ax.yaxis.set_major_formatter(NullFormatter())
        
def HistFmt(ax,xKey,row,col,Labels):
    """
    This function formats a histogram panel
    Parameters
            ----------
            ax : matplotlib axes
                The panel in the figure to be plotted
            xKey : string 
                The x-parameter
            row : integer
                The current row of the panel
            col : integer 
                The current column of the panel
            Labels : dictionary
                A dictionary of panel labels
    """
    #   Along the x-axis set the label as N
    if row==0:
        ax.set_xlabel('N')
    else:
        #   Otherwise remove the ticks
       ax.xaxis.set_major_formatter(NullFormatter())
    #   If in the leftmost column, add a label
    if col==0:
        ax.set_ylabel(Labels[xKey])
    else:
        #   Otherwise, remove the ticks
       ax.yaxis.set_major_formatter(NullFormatter())



def PlotIni():
    """
    This function initializes a plot and contains general use plotting parameters/definitions
        ------
        Returns
            fig : matplotlib figure instance
                The actual canvas for the plot
            
            PltOpts: dictionary
                A useful set of definitions for plotting purposes
    """
    #   Start by updating the base plot parameters
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
    #   Initialize the figure
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    #   Set the basic panel definitions
    base=-0.1
    left=0.1
    w=0.45
    h=w
    buf=0.15
    cW=0.02
    hbuf=0.2
    #   Store them in a dictionary and return both the figure and the plotting options.
    PltOpts={'base':base,'left':left,'w':w,'h':h,'buf':buf,'cW':cW,'hbuf':hbuf}
    return fig, PltOpts
