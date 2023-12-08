import sys as sys
import os as os
import copy as copy
import numpy as np


import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

from . import CubeAnalysis as CA
from . import PVPlotFncs as PVP

from datetime import date

def WriteBootstrappedFitOutputFile_Text(GalaxyDict):
    print(GalaxyDict.keys())
    print(GalaxyDict['WRKP_ResultsFolder'])

    #   name the output file
    BootstrapOutputFile=GalaxyDict['WRKP_ResultsFolder']+GalaxyDict['ObjName']+"_BSModel.txt"
    
    
    Model=GalaxyDict['BestFitModel']
        #   First get the date that the analysis is being run
    now = date.today()
    Date= str(now.year) + "-" + str(now.month) + "-" +  str(now.day)
    #   Add the date to the dictionary
    Model['Date']=Date
    
    
    
    #   Write out a the file header
    f=open(BootstrapOutputFile,"w")
    HeaderStr="Object:\t "+GalaxyDict['ObjName']+"\n"
    HeaderStr+="Date: \t" + Date+"\n"
    HeaderStr+="KINVER: \t 3KIDNASv1\n"
    #HeaderStr+="Version: \t" + str(FittingOptions['ModelVersion'])+"\n"
    #HeaderStr+="SBID: \t" + str(AvgModel['SBID'])+"\n"
    #HeaderStr+="SRCTR: \t" + str(AvgModel['SRCTR'])+"\n"
    #HeaderStr+="SRCVER: \t" + str(AvgModel['SRCVER'])+"\n"
    #HeaderStr+="KINTR: \t" + str(AvgModel['KINTR'])+"\n"
    #HeaderStr+="KINVER: \t" + str(AvgModel['KINVER'])+"\n"
    #HeaderStr+="Flag: \t"+str(AvgModel['FitFlags'])+"\n\n"
    #   Next write out the Bootstrap info
    HeaderStr+="\nnBootstraps\t"+str(GalaxyDict['nBootstraps'])+"\n"
    HeaderStr+="nBootstrap_Fits\t"+str(Model['nFits'])+"\n"
    
    #   And now write out the S/N information
    HeaderStr+="\nCube Noise\n"
    HeaderStr+="RMS (mJy/beam) \t"+str(Model['RMS'][0])+"\n"
    HeaderStr+="SN_Int \t"+str(Model['SN_Integrated'][0])+"\n"
    HeaderStr+="SN_Peak \t"+str(Model['SN_Peak'][0])+"\n"
    HeaderStr+="SN_Avg \t"+str(Model['SN_Avg'][0])+"\n"
    HeaderStr+="SN_Median \t"+str(Model['SN_Median'][0])+"\n"
    
    #   Next write out the geometric parameters
    #       Start with a brief header
    HeaderStr+="\nGeometry Parameters\n"
    HeaderStr+="Param Name\t Value \t Error \n"
    f.write(HeaderStr)
    #   And now add each parameter to the file
    GeoStr="X_model (pixels) \t" +str(Model['XCENTER'][0])+"\t\t"+str(Model['XCENTER_ERR'][0])+"\n"
    GeoStr+="Y_model (pixels) \t" +str(Model['YCENTER'][0])+"\t\t"+str(Model['YCENTER_ERR'][0])+"\n"
    GeoStr+="RA_model (degrees) \t" +str(Model['RA'][0])+"\t\t"+str(Model['RA_ERR'][0])+"\n"
    GeoStr+="DEC_model (degrees) \t" +str(Model['DEC'][0])+"\t\t"+str(Model['DEC_ERR'][0])+"\n"
    GeoStr+="Inc_model (degrees) \t" +str(Model['INCLINATION'][0])+"\t\t"+str(Model['INCLINATION_ERR'][0])+"\n"
    GeoStr+="PA_model (degrees) \t" +str(Model['POSITIONANGLE'][0])+"\t\t"+str(Model['POSITIONANGLE_ERR'][0])+"\n"
    GeoStr+="PA_model_g (degrees) \t" +str(Model['PA_GLOBAL'])+"\t\t"+str(Model['PA_GLOBAL_ERR'])+"\n"
    GeoStr+="VSys_model (km/s) \t" +str(Model['VSYS'][0])+"\t\t"+str(Model['VSYS_ERR'][0])+"\n"
    GeoStr+="VDisp_model (km/s) \t" +str(Model['VDISP'][0])+"\t\t"+str(Model['VDISP_ERR'][0])+"\n\n"
    f.write(GeoStr)
    #   Then do the rotation curve
  
    nR=len(Model['VROT'])
    #       Add a brief header first
    HeaderStr="Rotation Curve\n"
    HeaderStr+="nR=\t"+str(nR)+"\n"
    HeaderStr+="Rad \t\t VROT_model \t e_VRot_model  \n"
    HeaderStr+="('') \t\t (km/s) \t\t    (km/s) "
    f.write(HeaderStr)
    #       Now start the profile string
    ProfileStr="\n"
    #       And loop through the radial profile adding the rotation curve
    for i in range(nR):
        ProfileStr+=str(Model['R'][i])+"\t\t"+str(Model['VROT'][i])+"\t\t"+ str(Model['VROT_ERR'][i])+"\n"
    f.write(ProfileStr)
    
    #   Finally do the same thing for the surface density
    #       Again write a small header
    HeaderStr="\nSurface Density Profile\n"
    HeaderStr+="nR=\t"+str(nR)+"\n"
    HeaderStr+="Rad  \t\t SD_FO_model \t e_SD_FO_inc_model \n"
    HeaderStr+="('')\t\t (Msol/pc^2) \t (Msol/pc^2)  "
    f.write(HeaderStr)
    #   Then go through the full surface density profile
    ProfileStr="\n"
    for i in range(nR):
        ProfileStr+=str(Model['R_SD'][i])+"\t\t"+str(Model['SURFDENS_FACEON'][i])+"\t\t"+str(Model['SURFDENS_FACEON_ERR'][i])+"\n"
    f.write(ProfileStr)
    
    
    
    #       Add another set for the projected SD profile
    ExtendSDProfile=GalaxyDict['ExtendedSDProfile']
    nR_SD=len(ExtendSDProfile['R_SD'])
    HeaderStr="\nProjected Surface Density Profile from Mom0 map\n"
    HeaderStr+="nR=\t"+str(nR_SD)+"\n"
    HeaderStr+="Rad  \t\t SD_projected_model \t e_SD_projected_model \t SD_FO_projected_model \t e_SD_FO_inc_projected_model \n"
    HeaderStr+="('')\t\t (Msol/pc^2) \t (Msol/pc^2)  \t (Msol/pc^2) \t (Msol/pc^2) "
    f.write(HeaderStr)
    #   Then go through the full surface density profile
    ProfileStr="\n"
    for i in range(nR_SD):
        ProfileStr+=str(round(ExtendSDProfile['R_SD'][i],2))+"\t\t"+str(round(ExtendSDProfile['SURFDENS'][i],2))+"\t\t"+str(round(ExtendSDProfile['SURFDENS_ERR'][i],2))+"\t\t"+str(round(ExtendSDProfile['SURFDENS_FACEON'][i],2))+"\t\t"+str(round(ExtendSDProfile['SURFDENS_FACEON_ERR'][i],2))+"\n"
    f.write(ProfileStr)

    #       Finish off with the Scaling Relation results
    ScalingDict=GalaxyDict['ScalingDict']
    ScalingStr="\nRHI Extraction Method (0=3D profile, 1=2D map, 2=2D map End point)\n"
    ScalingStr+=str(ScalingDict['SDMethod'])+"\n"
    ScalingStr+="RHI and limits (arcsec) \n"
    ScalingStr+=str(round(ScalingDict['RHI_CorrArr'][1],2))+"\t\t"+str(round(ScalingDict['RHI_CorrArr'][0],2))+"\t\t"+str(round(ScalingDict['RHI_CorrArr'][2],2))+"\n"
    ScalingStr+="RHI and limits (kpc) \n"
    ScalingStr+=str(round(ScalingDict['RHI_kpc'][1],2))+"\t\t"+str(round(ScalingDict['RHI_kpc'][0],2))+"\t\t"+str(round(ScalingDict['RHI_kpc'][2],2))+"\n"
    ScalingStr+="\nVHI Extraction flag (0=interpolation, 1=extrapolation\n"
    ScalingStr+=str(ScalingDict['VHIFlag'])+"\n"
    ScalingStr+="VHI and limits (km/s) \n"
    ScalingStr+=str(round(ScalingDict['VHIArr'][1],2))+"\t\t"+str(round(ScalingDict['VHIArr'][0],2))+"\t\t"+str(round(ScalingDict['VHIArr'][2],2))+"\n"
    
    f.write(ScalingStr)
    
    #   Finally close the kinematic model file.
    f.close()
    
def SavePVDiagrams(GalaxyDict,GeneralDict):
    #   Set a list of file suffixes
    FileSuffixes=["_PVMajor_Data.fits","_PVMinor_Data.fits","_PVMajor_Model.fits","_PVMinor_Model.fits"]
    #   Set the base name of the files
    BaseName=GalaxyDict['WRKP_ResultsFolder']+GalaxyDict['ObjName']
    
    #   Set the model
    Model=GalaxyDict['BestFitModel']
    #   Load in the data cube and the model cube
    DataCube=CA.BasicCubeAnalysis(GalaxyDict['CubeName'])
    ModelCube=CA.BasicCubeAnalysis(Model['ModelCube'])

    #   From the model get the center point and base PA
    CenterX=Model['XCENTER'][0]
    CenterY=Model['YCENTER'][0]
    PA=Model['POSITIONANGLE'][0]
    #   Now put everything into a set of lists
    CUseArr=[DataCube,DataCube,ModelCube,ModelCube]
    CentArr=[CenterX,CenterY,CenterX,CenterY]
    PAArr=[PA,PA+90.,PA,PA+90.]
    
    
    #   The PV routine needs the beamsize in pixels to figure out how things should be cut
    BeamSize_Pix=DataCube['CubeHeader']['BMAJ']/np.abs(DataCube['CubeHeader']['CDELT1'])
    #   It also needs the central pixel
    CentPix=[Model['XCENTER'][0],Model['YCENTER'][0]]
    
    #   Now loop through the cubes and make each PV file
    for i in range(4):
        FName=BaseName+FileSuffixes[i]
        CUse=CUseArr[i]
        CentUse=CentArr[i]
        PAUse=PAArr[i]-180.
        #   Correct the PA to be in 0-360
        if PAUse > 360.:
            PAUse=PAUse-360.
        elif PAUse < 0.:
            PAUse=PAUse+360.

        
        PV=CA.ConstructModelBasedPVDiagram(CUse['Data'],PAUse,Model,BeamSize_Pix,DataCube)
        PVDict={'PV':PV,'CubeVels':DataCube['CubeVels']}
        #   Now generate a fits object with the data
        hdu = fits.PrimaryHDU(PVDict['PV'].T)
        #   Go through and set a lot of the keywords
        
        hdu.header.set('CRPIX1',np.shape(PV)[0]/2)
        hdu.header.set('CRVAl1',0.)
        hdu.header.set('CDELT1',np.abs(DataCube['CubeHeader']['CDELT1']))
        hdu.header.set('CTYPE1','Offset')
        
        ObjKeys2=['CUNIT1','CRPIX2','CRVAL2','CDELT2','CTYPE2','CUNIT2','BUNIT','BMAJ','BMIN','BPA','OBJECT']
        MatchKeys=['CUNIT1','CRPIX3','CRVAL3','CDELT3','CTYPE3','CUNIT3','BUNIT','BMAJ','BMIN','BPA','OBJECT']
        for ii in range(len(ObjKeys2)):
            hdu.header.set(ObjKeys2[ii],DataCube['CubeHeader'][MatchKeys[ii]])

        hdu.header.set('CTYPE1','3KIDNAS')
        hdul = fits.HDUList([hdu])
        hdul.writeto(FName,overwrite=True)

        
