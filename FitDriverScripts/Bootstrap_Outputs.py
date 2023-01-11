import sys as sys
import os as os
import copy as copy
import numpy as np



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
    HeaderStr+="KINVER: \t WRKPv1\n"
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

    #   Finally close the kinematic model file.
    f.close()
    
