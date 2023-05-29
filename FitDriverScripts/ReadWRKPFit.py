import sys as sys
import os as os
import copy as copy
import numpy as np

def ReadWRKPOutputFile(GeneralDict,GalaxyDict):
    print("Reading in WRKP results file")
    #   Get the name of the WRKP Output File and model cube file
    WRKP_Results_File,WRKP_CubeFile=GetOutputFileName(GalaxyDict)
    #   Now load in the best fitting model file
    GalaxyModel=LoadBestFitModelFile(WRKP_Results_File)
    #   Also store the model cube file name in the GalaxyModel dictionary
    GalaxyModel['ModelCube']=WRKP_CubeFile
    return GalaxyModel

def GetOutputFileName(GalaxyDict):
    TargFolder=GalaxyDict['TargFolderU']+"/"+GalaxyDict['ObjNameU']+"/"
    FName=TargFolder+GalaxyDict['ObjNameU']+"_AvgModel_v1.txt"
    CName=TargFolder+GalaxyDict['ObjNameU']+"_AverageModel_v1.fits"
    
    return FName,CName


def LoadBestFitModelFile(FileName):
    print("Loading Avg Model",FileName)
    if os.path.isfile(FileName) == False:
        FitAchieved=False
        AvgDict={'FITAchieved':FitAchieved}
        return AvgDict
    else:
        FitAchieved=True
        #   Read the file
        Fit=open(FileName,"r")
        Lines=Fit.readlines()
        #   Load in the SN estimates
        
        #   Set the lines that start the noise, geometry, and profile values
        NoiseLineIndx=5
        ProfileStartLine=25
        GeoLineStart=13
        #   Figure out the number of profile values
        nR=len(Lines)-ProfileStartLine
        #   Load in the geometric parameters
        XCent,XErr=GeoLineAssign(Lines[GeoLineStart],nR)
        YCent,YErr=GeoLineAssign(Lines[GeoLineStart+1],nR)
        RA,RAErr=GeoLineAssign(Lines[GeoLineStart+2],nR)
        DEC,DECErr=GeoLineAssign(Lines[GeoLineStart+3],nR)
        Inc,IncErr=GeoLineAssign(Lines[GeoLineStart+4],nR)
        PA,PAErr=GeoLineAssign(Lines[GeoLineStart+5],nR)
        VSys,VSysErr=GeoLineAssign(Lines[GeoLineStart+6],nR)
        VDisp,VDispErr=GeoLineAssign(Lines[GeoLineStart+7],nR)
        #print("Geo load test", XCent[0],YCent[0],RA[0],DEC[0],Inc[0],PA[0],VSys[0])

        #   Initialize the radial profiles
        R=np.zeros(nR)
        VRot=np.zeros(nR)
        VErr=np.zeros(nR)
        VErrInc=np.zeros(nR)
        SD=np.zeros(nR)
        SDErr=np.zeros(nR)
        #   Load in the radial profile values
        for i in range(nR):
            j=i+ProfileStartLine
            R[i],VRot[i],VErr[i],VErrInc[i],SD[i],SDErr[i]=ProfileLineAssign(Lines[j])
            #print(R[i],VRot[i])
        #   Close the file
        Fit.close()
        #   Assign values to a standard TR dictionary
        AvgDict={'R':R,'R_SD':R,'XCENTER':XCent, 'XCENTER_ERR':XErr\
            ,'YCENTER':YCent,'YCENTER_ERR':YErr\
            ,'INCLINATION':Inc,'INCLINATION_ERR':IncErr\
            ,'POSITIONANGLE':PA,'POSITIONANGLE_ERR':PAErr\
            ,'VSYS':VSys,'VSYS_ERR':VSysErr\
            ,'VROT':VRot,'VROT_ERR':VErr\
            ,'SURFDENS':SD,'SURFDENS_ERR':SDErr\
            ,'RA':RA,'RA_ERR':RAErr
            ,'DEC':DEC,'DEC_ERR':DECErr
            ,'VDISP':VDisp,'VDISP_ERR':VDispErr
            ,'FITAchieved':FitAchieved,'CHI2':-1}
        AvgDict['SURFDENS_FACEON']=AvgDict['SURFDENS']
        AvgDict['SURFDENS_FACEON_ERR']=AvgDict['SURFDENS_ERR']
        #       Add in the noise values to the dictionary
        NoiseKeys=["RMS","SN_Integrated","SN_Peak","SN_Avg","SN_Median"]
        for key in NoiseKeys:
            AvgDict[key]=NoiseLineAssign(Lines[NoiseLineIndx])
            NoiseLineIndx+=1
        return AvgDict

def GeoLineAssign(Line,nR):
    S=Line.split()
    Val=np.full(nR,float(S[2].strip()))
    Err=np.full(nR,float(S[3].strip()))
    return Val,Err
    
def NoiseLineAssign(Line):
    S=Line.split()
    Val=np.full(1,float(S[-1].strip()))
    return Val
    
def ProfileLineAssign(Line):
    #S=re.split('\t\t\t|\t\t|\t',Line)
    S=Line.split()
    R=float(S[0].strip())
    V=float(S[1].strip())
    VErr=float(S[2].strip())
    VErrInc=float(S[3].strip())
    SD=float(S[4].strip())
    SDErr=float(S[5].strip())
    return R,V,VErr,VErrInc,SD,SDErr


def LoadBootstrappedFit(FileName):
    print("Loading Avg Model",FileName)
    if os.path.isfile(FileName) == False:
        FitAchieved=False
        AvgDict={'FITAchieved':FitAchieved}
        return AvgDict
    else:
        FitAchieved=True
        #   Read the file
        Fit=open(FileName,"r")
        Lines=Fit.readlines()
        #   Load in the SN estimates
        
        #   Set the lines that start the noise, geometry, and profile values
        NoiseLineIndx=8
        ProfileStartLine=30
        GeoLineStart=16
        #   Figure out the number of profile values
        nR=nRLineSet(Lines[27])
        #   Load in the geometric parameters
        XCent,XErr=GeoLineAssign(Lines[GeoLineStart],nR)
        YCent,YErr=GeoLineAssign(Lines[GeoLineStart+1],nR)
        RA,RAErr=GeoLineAssign(Lines[GeoLineStart+2],nR)
        DEC,DECErr=GeoLineAssign(Lines[GeoLineStart+3],nR)
        Inc,IncErr=GeoLineAssign(Lines[GeoLineStart+4],nR)
        PA,PAErr=GeoLineAssign(Lines[GeoLineStart+5],nR)
        PA_G,PA_GErr=GeoLineAssign(Lines[GeoLineStart+6],nR)
        VSys,VSysErr=GeoLineAssign(Lines[GeoLineStart+7],nR)
        VDisp,VDispErr=GeoLineAssign(Lines[GeoLineStart+8],nR)
        #print("Geo load test", XCent[0],YCent[0],RA[0],DEC[0],Inc[0],PA[0],VSys[0])

        #   Initialize the radial profiles
        R=np.zeros(nR)
        VRot=np.zeros(nR)
        VErr=np.zeros(nR)
        VErrInc=np.zeros(nR)
        SD=np.zeros(nR)
        SDErr=np.zeros(nR)
        #   Load in the radial profile values
        for i in range(nR):
            j=i+ProfileStartLine
            R[i],VRot[i],VErr[i]=VProfileLineAssign(Lines[j])
        
        nR_SD=nRLineSet(Lines[j+3])
        ProfileStartLine=j+6
        R_SD=np.zeros(nR_SD)
        SD=np.zeros(nR_SD)
        SDErr=np.zeros(nR_SD)
        SD_FO=np.zeros(nR_SD)
        SD_FO_Err=np.zeros(nR_SD)
        for i in range(nR_SD):
            j=i+ProfileStartLine
            R_SD[i],SD_FO[i],SD_FO_Err[i]=SDProfileLineAssign(Lines[j])
            
            
        #   Close the file
        Fit.close()
        #   Assign values to a standard TR dictionary
        AvgDict={'R':R,'R_SD':R,'XCENTER':XCent, 'XCENTER_ERR':XErr\
            ,'YCENTER':YCent,'YCENTER_ERR':YErr\
            ,'INCLINATION':Inc,'INCLINATION_ERR':IncErr\
            ,'POSITIONANGLE':PA,'POSITIONANGLE_ERR':PAErr\
            ,'VSYS':VSys,'VSYS_ERR':VSysErr\
            ,'VROT':VRot,'VROT_ERR':VErr\
            ,'SURFDENS':SD_FO,'SURFDENS_ERR':SD_FO_Err\
            ,'RA':RA,'RA_ERR':RAErr
            ,'DEC':DEC,'DEC_ERR':DECErr
            ,'VDISP':VDisp,'VDISP_ERR':VDispErr
            ,'FITAchieved':FitAchieved,'CHI2':-1}
        AvgDict['SURFDENS_FACEON']=AvgDict['SURFDENS']
        AvgDict['SURFDENS_FACEON_ERR']=AvgDict['SURFDENS_ERR']
        #       Add in the noise values to the dictionary
        NoiseKeys=["RMS","SN_Integrated","SN_Peak","SN_Avg","SN_Median"]
        for key in NoiseKeys:
            AvgDict[key]=NoiseLineAssign(Lines[NoiseLineIndx])
            NoiseLineIndx+=1
        return AvgDict

def nRLineSet(Line):
    lSplit=Line.split('=')
    nR=int(lSplit[1].strip())
    return nR

    
def VProfileLineAssign(Line):
    #S=re.split('\t\t\t|\t\t|\t',Line)
    S=Line.split()
    R=float(S[0].strip())
    V=float(S[1].strip())
    VErr=float(S[2].strip())
    #SD=float(S[4].strip())
    #SDErr=float(S[5].strip())
    return R,V,VErr
    
def SDProfileLineAssign(Line):
    #S=re.split('\t\t\t|\t\t|\t',Line)
    S=Line.split()
    R=float(S[0].strip())
    SD_FO=float(S[1].strip())
    SD_FO_Err=float(S[2].strip())
    return R,SD_FO,SD_FO_Err
