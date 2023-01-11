#!/usr/bin/env python3
import numpy as np
from . import  BBaroloResults as BR
from . import FATResults as FR
from . import ParameterPlots as PP
from . import CubeInformation as CI
from . import PVPlots as PVP
from . import BootstrapParameterPlots as BPP

def BootstrapErrors(CatID,StrDict,WallCat):
    Name=WallCat.name[CatID]
    NameSuffix=Name.split(' ')[1]
    print("Comparing Fits for galaxy ", CatID,Name)
    print(StrDict['FitPaths'])


    CubeInfo=CI.CubeAnalysis(CatID,StrDict,WallCat)


    nFits=np.shape(StrDict['FitPaths'])[0]


    TRParams=[]
        #for i in range(7,8,1):
    for i in range(1):
        if StrDict['FitPaths'][i][1] == 0:
            TRParams.append(BR.GetBBaroloFit(CatID,WallCat,StrDict['FitPaths'][i][0]
                                             ,StrDict['FitPaths'][i][2],StrDict['FitPaths'][i][4],CubeInfo))
        elif StrDict['FitPaths'][i][1] == 1:
            P1,P2=FR.GetFATFit(CatID,WallCat,StrDict['FitPaths'][i][0]\
                               ,StrDict['BaseFileName'],CubeInfo,StrDict['FitPaths'][i][4])
            TRParams.append(P1)
            TRParams.append(P2)

    nEstimates=StrDict['FitPaths'][1][5]
    print("number of estimates", nEstimates)

    BootstrapParams=[]
    for i in range(100):
        j=i+1
        FitFolder=StrDict['FitPaths'][1][0]+str(j)+"/"+NameSuffix+"/"
        print(FitFolder)
        ParamFile=[FitFolder+"ringlog2.txt"]
        DensFile=FitFolder+"densprof.txt"
        MainParams,FitAchieved=BR.LoadBBaroloParamFile(ParamFile)
        RR,SD=BR.LoadBBaroloSurfDens(DensFile)
        BBaroloDict=BR.BBaroloParamsToDict(MainParams,SD,CubeInfo)
        BBaroloDict['R_SD']=RR
        BootstrapParams.append(BBaroloDict)

    BootStrapQuantities=BootStrapAvg(BootstrapParams)


    BPP.MakeAllPlots(CatID,WallCat,StrDict,TRParams,CubeInfo,BootStrapQuantities)
#    PVP.MakePVPlots(CatID,WallCat,StrDict,TRParams,CubeInfo)


def BootStrapAvg(BootstrapParams):
    print(len(BootstrapParams))
    nEstimates=len(BootstrapParams)
    BootstrapSample={'R':BootstrapParams[0]['R'],'R_SD':BootstrapParams[0]['R_SD']}

    Key='VROT'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='SURFDENS'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='INCLINATION'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='POSITIONANGLE'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='VSYS'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='VDISPERSION'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='XCENTER'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)
    Key='YCENTER'
    BootstrapSample[Key]=AvgQuantity(BootstrapParams,nEstimates,Key)

    return BootstrapSample


def AvgQuantity(BootstrapParams,nEstimates,Key):
    Avg=np.zeros(len(BootstrapParams[0]['R']))
    for i in range(nEstimates):
        Avg+=BootstrapParams[i][Key]
    Avg=Avg/nEstimates
    Err=np.zeros(len(BootstrapParams[0]['R']))
    for i in range(nEstimates):
        Err+=(BootstrapParams[i][Key]-Avg)**2.
    Err=np.sqrt(Err/nEstimates)
    return Avg,Err



