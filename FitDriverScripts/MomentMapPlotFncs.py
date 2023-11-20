import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import copy as copy
from astropy.visualization.wcsaxes import SphericalCircle


import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable


import third_party.CosmosCanvas
import third_party.CosmosCanvas.velmap as vmap
from third_party.colourspace import maps


"""
    This module contains a number of plotting functions that are useful for both diagnostic and average model plots.  This module contains routines focused on making moment maps.  It contains the routines:

    MomentPlot --> This function makes a basic moment map panel for a plot
    AddCenterToMomMap --> This function adds a point to a moment map panel indicating the center
    AddArrowToMomMap --> This function adds an arrow to a moment map panel indicating the position angle
    AddVelContoursToMomentPlot --> This function adds a set of velocity contours to a moment map based on some other moment map.
    MakeMomMap --> This function makes a set of moment maps for some cube.
    QuickMaskCube --> This does a quick masking of a data cube by a galaxy's SoFiA mask cube.
    
    DetermineExtent --> This function determines the extent for a moment map.
    FormatMomMap --> This function formats and cleans up the moment map panels.
"""

def MakeAllMomentMapPlots(fig,Model,DataCube,ModelCube,PltOpts,GalaxyDict,GeneralDict):
    print("Making Moment Map Plots")
    #   Start by making moment maps for the model and data cube
    for i in range(2):
        if i==0:
            CUse=DataCube
            MSwitch=1
        elif i==1:
            CUse=ModelCube
            MSwitch=0
        CUse=QuickMaskCube(CUse,GalaxyDict)
        CUse['Mom0'],CUse['Mom1']=MakeAllMomMaps(CUse,MSwitch)
        
    #   We'll also need a grid of X and Y for the making the maps
    XX,YY=MakeGrid(DataCube)

    #   Now we'll need to set up the dimensions for the moment maps
    Size=np.shape(DataCube['Mom0'])
    PltOpts['hU']=Size[0]/Size[1]*PltOpts['w']
    PltOpts['base']=PltOpts['base']-(PltOpts['hU']+PltOpts['hbuf'])
    PltOpts['lU']=PltOpts['left']
    #   Now we can make the Mom0 panel
    MakeMomPanel(fig,Model,DataCube,ModelCube,PltOpts,0,XX,YY)
    
    #   We can also make the Moment 1 panel
    PltOpts['lU']=PltOpts['left']+PltOpts['buf']+PltOpts['w']
    MakeMomPanel(fig,Model,DataCube,ModelCube,PltOpts,1,XX,YY)

def MakeMomPanel(fig,Model,DataCube,ModelCube,PltOpts,Moment,XX,YY):
    #   Here we'll start with making the canvas
    placement=[PltOpts['lU'],PltOpts['base'],PltOpts['w'],PltOpts['hU']]
    #           Use the astropy WCS projection for the panel
    ax=fig.add_axes(placement,projection=DataCube['CubeWCS'],slices=('x', 'y', 0))
    #       We'll also make the axis for the colorbar
    placement=[PltOpts['lU'],PltOpts['base']+PltOpts['hU'],PltOpts['w'],PltOpts['cW']]
    cax=fig.add_axes(placement)
    #   Now we'll want to format the plot
    ax=WCSPanelFmt(ax)
    #   Set up the colormap
    minV,maxV,CMap=SetCMap(Model,DataCube,ModelCube,Moment)
    #   With these we can make the moment map
    if Moment==0:
        MKey='Mom0'
        cLabel=r"$\Sigma$ (M$_{\odot}$ pc$^{-2})$"
    elif Moment==1:
        MKey='Mom1'
        cLabel=r"V$_{\rm{los}}$ (km s$^{-1})$"
        ax.set_facecolor('k')

    im=ax.pcolormesh(XX,YY,DataCube[MKey],cmap=CMap,shading='auto',vmin=minV,vmax=maxV)
    cbar=fig.colorbar(im,cax=cax,orientation='horizontal')
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_ticks_position('top')
    cax.set_xlabel(cLabel)
    
    
    
    #   Add a coordinate grid
    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    #   Add an ellipse showing the beam size to the moment 0 map
    if Moment==0:
        ESize=DataCube['CubeHeader']['BMaj']/np.abs(DataCube['CubeHeader']['CDElT1'])
        
        LocX=7
        LocY=7
        #   Get the ellipse
        Ell=Ellipse([LocX,LocY], ESize, ESize, 0.,edgecolor='cyan',facecolor='none',lw=2)
        #   Add the ellipse to the map.
        ax.add_patch(Ell)
    elif Moment==1:
        AddArrowToMomMap(ax,Model,DataCube)
        AddVelContoursToMomentPlot(ax,ModelCube,Model,XX,YY)
    #   For both figures, add the center point
    AddCenterToMomMap(ax,Model)
    
    
    
def WCSPanelFmt(ax):
    #   Get the longitude and latitude
    lon = ax.coords[0]
    lat = ax.coords[1]
    #   Set the axis label
    lon.set_axislabel("RA")
    lat.set_axislabel("DEC")
    #   And set minor ticks
    lon.set_ticks(number=4)
    lat.set_ticks(number=4)
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(5)
    lon.set_minor_frequency(5)
    #   make sure that we don't have overlapping ticks
    lon.set_ticklabel(exclude_overlapping=True)
    lat.set_ticklabel(exclude_overlapping=True)
    
    return ax
    
def SetCMap(Model,DataCube,ModelCube,Moment):
    #   We want different colormaps for the different moments
    if Moment==0:
        CMap='gray'
        #   For the moment zero map, the limits should be based on the data cube moment map
        minV=0.
        maxV=1.05*np.nanmax(DataCube['Mom0'])
    elif Moment==1:
        VS=Model['VSYS'][0]
        dV=np.nanmax(Model['VROT'])*np.sin(Model['INCLINATION'][0]*np.pi/180.)
        minV=VS-dV
        maxV=VS+dV
        CMap=vmap.create_cmap_doubleVelocity(minV,maxV,div=VS,Cval_max=35)
        
        
    return minV,maxV,CMap


def QuickMaskCube(ObsCube,ObjDict):
    """
        This function masks the observed data cube by the galaxy's SoFiA mask file and stores it in the ObsCube dictionary.
    """
    #   First get rid of any nan's
    ObsCube['Data']=np.nan_to_num(ObsCube['Data'])
    #   Now open up the SoFiA masked cube --> note that this means the 'observed cube' must be at the native 4 km/s resolution.
    MaskCube=fits.open(ObjDict['MaskName'])
    #   Create a 'MaskedData' array for use in making moment maps
    ObsCube['MaskedData']=ObsCube['Data']*MaskCube[0].data
    #   Close the mask file
    MaskCube.close()
    #   Return the observed cube dictionary.
    return ObsCube

def MomentPlot(fig,placement,Moment,ObsCube,XC,YC,GalaxyDict,GeneralDict,AvgModel=False):
    """
        This function makes a moment map panel on a larger canvas.
    """
    #   First set up the moment map panel on the canvas using the placement rectangle.
    ax=fig.add_axes(placement,projection=ObsCube['CubeWCS'],slices=('x', 'y', 0))
    #   Next mask the 'observed cube' by the SoFiA mask -- the observed cube must be at 4 km/s resolution.
    ObsCube=QuickMaskCube(ObsCube,GalaxyDict)
    
    #   Make the appropriate moment map (indicated by the Moment integer).
    MomMapData=MakeMomMap(ObsCube,Moment,1)
        #   Copy the moment map into a useable map where the units can be adjusted
    MomMapDataUse=copy.copy(MomMapData)
    #   When doing the moment 0 map, adjust the units
    if Moment ==0:
        MomMapDataUse=Moment0MapUnitConvert(MomMapDataUse,ObsCube)
    #   Set the base minimum and maximum values colormap on the map
    minval=np.nanmin(MomMapDataUse)
    maxval=1.05*np.nanmax(MomMapDataUse)
    #   Adjust the moment 1 colormap limits
    if Moment == 1:
        #   If we don't have a TR model to use when making the map, use the cube size
        minval=ObsCube['CubeVels'][-1]/1000.
        maxval=ObsCube['CubeVels'][0]/1000.
        #   If we do have a model, use the inclination corrected velocity to get the limits on the colormap
        if AvgModel != False:
            VS=AvgModel['VSYS'][0]
            dV=np.nanmax(AvgModel['VROT'])*np.sin(AvgModel['INCLINATION'][0]*np.pi/180.)
            minval=VS-1.1*dV
            maxval=VS+1.1*dV
    #   Set the colored maps
    if Moment == 0:
        Cmap='gray'
    elif Moment == 1:
        Cmap='jet'
        ax.set_facecolor('black')
        VS=(maxval+minval)/2.
        Cmap=vmap.create_cmap_doubleVelocity(minval,maxval,div=VS,Cval_max=35)


        
    elif Moment ==2 :
        Cmap='plasma'
    #   Format the moment map panel
    ExtentDict=FormatMomPlot(ax,Moment,XC,YC,ObsCube)
    #   Make the basic moment map using imshow to keep the proper aspect ratio
    im1=ax.imshow(MomMapDataUse,cmap=Cmap,vmin=minval,vmax=maxval,origin='lower',extent=ExtentDict['Extent'])
    #   For the moment 0 map, add in an ellipse to represent 1 beam
    if Moment ==0:
        #   First set the size of the beam from pixels to arcseconds
        ESize=5.*abs(ExtentDict['dX'])
        #   Set the center of the ellipse to pixel (7,7)
        Loc=7.
        #   Adjust the location to be coordinates in arcseconds using the 'ExtentDict' generated in DetermineExtent via FormatMomPlot
        LocX=(Loc-ExtentDict['RefPixX'])*ExtentDict['dX']+ExtentDict['RefValX']
        LocY=(Loc-ExtentDict['RefPixY'])*ExtentDict['dY']+ExtentDict['RefValY']
        #   Get the ellipse
        Ell=Ellipse([LocX,LocY], ESize, ESize, 0.,edgecolor='cyan',facecolor='none',lw=2)
        #   Add the ellipse to the map.
        ax.add_patch(Ell)
    #  Add a colorbar to the axis
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar=fig.colorbar(im1, cax=cax, orientation='vertical')
    if Moment ==0:
        cbar.set_label(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    elif Moment == 1:
        cbar.set_label(r"$V_{\rm{los}}$ (km s$^{-1}$)")
    
    return ax,MomMapData

def DetermineExtent(XC,YC,ObsCube):
    """
        This function sets the extent for moment maps in arcseconds
    """
    #   Set up the X and Y definitions needed
    #       First get the reference pixel
    RefPixX=XC
    #       Set the center to 0
    RefValX=0.
    #       Get the size of each pixel in arcseconds -- this assumes CDELT is in arcseconds
    dX=ObsCube['CubeHeader']['CDELT1']*3600.
    #       Get the number of pixels in the x direction
    nX=ObsCube['CubeHeader']['NAXIS1']
    #   Repeat the steps for the Y direction
    RefPixY=YC
    RefValY=0.
    dY=ObsCube['CubeHeader']['CDELT2']*3600.
    nY=ObsCube['CubeHeader']['NAXIS2']
    
    #  Get the lower and upper limits in arcseconds for both X and Y
    LowX=(-RefPixX)*dX+RefValX
    HighX=((nX-1)-RefPixX)*dX+RefValX
    LowY=(-RefPixY)*dY+RefValY
    HighY=((nY-1)-RefPixY)*dY+RefValY
    #   Set the extent array to be (X_low,X_high,Y_low,Y_high)
    Extent=[LowX,HighX,LowY,HighY]
    #   Store the extents and calculated values into a dictionary in case the extra info is needed elsewhere.
    ExtentDict={'Extent':Extent,'nX':nX,'RefPixX':RefPixX,'RefValX':RefValX,'dX':dX,'nY':nY,'RefPixY':RefPixY,'RefValY':RefValY,'dY':dY}
    #   Return the extent dictionary
    return ExtentDict

    
def AddCenterToMomMap(ax,Model):
    """
        This function adds a black X mark to the center point of a moment map.  Since the maps are formatted in delta RA and delta DEC from a center point, this is always 0,0
    """
    XC=Model['XCENTER'][0]
    YC=Model['YCENTER'][0]
    PCol='magenta'
    ax.scatter(XC,YC,marker='x',color=PCol,linewidths=5.0,s=150)
    return ax
    
def AddArrowToMomMap(ax,Model,Cube):
    """
        This function adds an arrow to a moment map panel to represent the position angle
    """
    #   First set the center coordinates of the arrow as the center of the galaxy
    CentX=Model['XCENTER'][0]
    CentY=Model['YCENTER'][0]
    Angle=Model['POSITIONANGLE'][0]
    L=1.1*Model['R'][-1]/(np.abs(Cube['CubeHeader']['CDELT1'])*3600.)
    #   Adjust the angle due to flipping the X-axis in RA
    AngleUse=360.-Angle
    #   Use the length and the angle to get the length in x and y for the arrow
    dX=-L*np.cos((AngleUse+90.)*np.pi/180.)
    dY=L*np.sin((AngleUse+90.)*np.pi/180.)
    #   Set the arrow width
    awidth=1.
    #   Add the arrow to the panel
    ax.arrow(CentX,CentY,dX,dY,ls=':',color='white',shape='full',width=awidth,head_width=3*awidth)
    return ax
    
def MakeAllMomMaps(Cube,MaskSwitch):
    
    Mom0=MakeMomMap(Cube,0,MaskSwitch)
    Mom0=Moment0MapUnitConvert(Mom0,Cube)
    Mom1=MakeMomMap(Cube,1,MaskSwitch)
    return Mom0,Mom1
    
def MakeMomMap(Cube,Moment,MaskSwitch):
    """
        This function constructs a moment 0, moment 1, or moment 2 array from a cube.
    """
    #       Get the map size
    Size=(np.shape(Cube['Data'])[1],np.shape(Cube['Data'])[2])
    #       Initialize the moment map and flux tot arrays
    MomMap=np.zeros(Size)
    FTot=np.zeros(Size)
    #   Decide whether to use the full data cube or the masked data cube.
    if MaskSwitch == 0:
        DUse=Cube['Data']
    else:
        DUse=Cube['MaskedData']
    #   The moment 0, 1, and 2 maps all need to start with a total flux map
    for i in range(np.shape(Cube['Data'])[1]):
        for j in range(np.shape(Cube['Data'])[2]):
            FTot[i,j]=np.nansum(DUse[:,i,j])
    #       When using unmasked data, the low flux cells can cause problems, so artificially set them to zero
    for i in range(np.shape(Cube['Data'])[1]):
        for j in range(np.shape(Cube['Data'])[2]):
            if FTot[i,j] <1.e-7:
                FTot[i,j]=0.
    #   Now construct the moment maps
    if Moment == 0:
        #   The moment 0 map is simply the total flux map.
        MomMap=FTot
    elif Moment==1 or Moment==2:
        #   For the moment 1 or moment 2 maps, we need to get the flux weighted average velocity in each pixel
        #       Loop through all the pixels
        for i in range(np.shape(Cube['Data'])[1]):
            for j in range(np.shape(Cube['Data'])[2]):
                #Loop through the channels
                for k in range(np.shape(Cube['Data'])[0]):
                    #   Get the flux-weighted sum for each channel
                    MomMap[i,j]+=(Cube['CubeVels'][k]/1000.)*DUse[k,i,j]
            #   Set the pixels with no flux to NaN's
            if FTot[i,j] ==0.:
                FTot[i,j]=FTot[i,j]/0.
        #   Normalize the map by the fluxes to get the moment 1 map
        MomMap=MomMap/FTot
    #   For the moment 2 map, go back through the pixels one more time
    if Moment ==2:
        #   Copy the moment 1 map to a temporary array
        VAvgMap=copy.copy(MomMap)
        #   Reset the moment map array to zeros
        MomMap=np.zeros(Size)
        #   Loop through all pixels and channels
        for i in range(np.shape(Cube['Data'])[1]):
            for j in range(np.shape(Cube['Data'])[2]):
                for k in range(np.shape(Cube['Data'])[0]):
                    #   Sum up the flux weighted squared difference
                    MomMap[i,j]+=(Cube['CubeVels'][k]/1000.-VAvgMap[i,j])**Moment*DUse[k,i,j]
        #   Now get the moment 2 map by taking the root of the mean square differences
        MomMap=np.sqrt(MomMap/FTot)
    #   Return either the moment 0, moment 1, or moment 2 map
    return MomMap
    
def Moment0MapUnitConvert(MomMap,Cube):
    """
        This function adjusts the units on the moment 0 map to Msol/pc^2 from
        the initial Jy/beam from the summation
    """
    #   First get the channel size and the beam size (in pixel units)
    ChanWidth=np.abs(Cube['CubeHeader']['CDELT3']/1000.)
    BeamSize=Cube['CubeHeader']['BMAJ']/abs(Cube['CubeHeader']['CDELT1'])
    BeamArea3=2.*np.pi*(BeamSize/2.355)**2.
    #   Convert to Jy km/s pixel^-1
    MomMap=MomMap*ChanWidth/BeamArea3
    #   Get the pixel area in arcsec^2
    PixArea=(abs(Cube['CubeHeader']['CDELT1'])*3600.)**2.
    #   Convert to Jy km/s arcsec^-2
    MomMap/=PixArea
    #   Convert to Msol/pc^2
    SDConv=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    MomMap/=SDConv
    return MomMap

    
def AddVelContoursToMomentPlot(ax,Cube,Model,XX,YY):
    """
        This function adds velocity contours to a moment 1 map panel from some other cube
    """
    #   Figure out the velocity width of the cube
    #       Get the outermost velocity
    VOut=Model['VROT'][-1]
    #       Use the model to get the inclination corrected outermoste velocity
    VSinI=VOut*np.sin(Model['INCLINATION'][0]*np.pi/180.)
    #       Set the model profile width
    VWidth=2.*VSinI
    #   Set dV to be twice this width
    dV=1.0*VWidth
    #       Set the contour line type and linewidth
    lTypes=('--')
    LW=1.5

    #   Set the velocity steps for the contours
    delV=dV/5.
    #   Set the contour levels
    CLevels=np.zeros(11)
    j=0
    for i in range(-5,6):
        CLevels[j]=Model['VSYS'][0]+i*dV/7.
        j+=1
    #   Draw on the contours
    MomMap=Cube['Mom1']
    CCol='#DEB887'
    ax.contour(XX,YY,MomMap,colors=CCol,linewidths=LW,linestyles=lTypes,levels=CLevels)
    #       Add in a thicker contour with vsys
    lTypes=('-')
    CLevels=np.array([Model['VSYS'][0]])
    ax.contour(XX,YY,MomMap,colors=CCol,linewidths=2.5*LW,linestyles=lTypes,levels=CLevels)
    #   Return the panel
    return ax
    


def MakeGrid(Cube):
    
    Shape=np.shape(Cube['Mom0'])
    X=np.zeros(Shape[1])
    Y=np.zeros(Shape[0])
    
    for i in range(Shape[1]):
        X[i]=i
    for i in range(Shape[0]):
        Y[i]=i
    
    XX,YY=np.meshgrid(X,Y)
    return XX,YY
