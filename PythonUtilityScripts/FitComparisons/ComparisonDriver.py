#!/usr/bin/env python3
import numpy as np
from . import  BBaroloResults as BR
from . import FATResults as FR
from . import ParameterPlots as PP
from . import CubeInformation as CI
from . import PVPlots as PVP

def CompareFits(CatID,StrDict,WallCat):
    Name=WallCat.name[CatID]
    print("Comparing Fits for galaxy ", CatID,Name)
    print(StrDict['FitPaths'])


    CubeInfo=CI.CubeAnalysis(CatID,StrDict,WallCat)


    nFits=np.shape(StrDict['FitPaths'])[0]


    TRParams=[]
        #for i in range(7,8,1):
    for i in range(nFits):
        if StrDict['FitPaths'][i][1] == 0:
            TRParams.append(BR.GetBBaroloFit(CatID,WallCat,StrDict['FitPaths'][i][0]
                                             ,StrDict['FitPaths'][i][2],StrDict['FitPaths'][i][4],CubeInfo))
        elif StrDict['FitPaths'][i][1] == 1:
            P1,P2=FR.GetFATFit(CatID,WallCat,StrDict['FitPaths'][i][0]\
                               ,StrDict['BaseFileName'],CubeInfo,StrDict['FitPaths'][i][4])
            TRParams.append(P1)
            TRParams.append(P2)

    PP.MakeAllPlots(CatID,WallCat,StrDict,TRParams,CubeInfo)
#    PVP.MakePVPlots(CatID,WallCat,StrDict,TRParams,CubeInfo)
