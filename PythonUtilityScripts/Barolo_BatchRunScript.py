import numpy as np
import os




#DataFolder="J103522-244506_BootstrapSamples"
#BaseName="J103522-244506"
#BaseOutFolder="BBaroloBatchFits_J103522-244506"
#BBaroloParamFile="WALLABY_PS_Hya_DR1_103_BBarolo_Chandra_SMask_NoIni.par"


DataFolder="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT14/MCGSuite/J103523-281855Models"
BaseName="ba_3.501.mass_9.47666.inc_51.86.pa_313.2.veldisp_8.0.version_"
#BaseOutFolder="BBaroloBatchFits_MCG_J103523-244506_NoIni"
#BBaroloParamFile="MCG_J103523-281855_Chandra.par"
BaseOutFolder="BBaroloBatchFits_MCG_J103523-244506_Chandra_ConstantMask_NoAngLims_VSysFree"
BBaroloParamFile="MCG_J103523-281855_Chandra.par"


#Load in the basic parameter file

BParFile=open(BBaroloParamFile,"r")

BParams=BParFile.readlines()



BParFile.close()

TempParamFile="TempParams.par"

for i in range(0,100,1):
    TPF=open(TempParamFile,"w")
    #DataName=DataFolder+"/"+BaseName+"_"+str(i)+".fits"
    DataName=DataFolder+"/"+BaseName+str(i)+"/"+BaseName+str(i)+".fits"
    OutName=BaseOutFolder+"/"+str(i)+"/"
    print(DataName)
    print(OutName)
    
    BParams[0]="FITSFILE   " +DataName+"\n"
    BParams[1]="OUTFOLDER   "+OutName+"\n"

    BParOut=" "
    for j in range(len(BParams)):
        BParOut+=str(BParams[j])

    TPF.write(BParOut)
    TPF.close()

    os.system("./Programs/BBarolo -p "+TempParamFile)
