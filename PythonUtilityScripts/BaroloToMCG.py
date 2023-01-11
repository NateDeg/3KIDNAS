import numpy as np



def ReadBarolo(Path):
    file=Path+"ringlog2.txt"
    BBaroloParams=np.loadtxt(file,skiprows=1,usecols=(1,2,3,4,5,7,8,9,10,11,12))
    file=Path+"densprof.txt"
    BBaroloParams[:,6]=np.loadtxt(file,skiprows=14,usecols=7)
    return BBaroloParams

def WriteMCGTiltedRingInput(ModelParams,InputFileName):
    #   Column Format
    HeaderStr="# The format of the input file — 1=row, 2= column \n 2 \n"
    #   Get the number of ring
    nRings=np.shape(ModelParams)[0]
    #   Add the number of strings to the input file
    HeaderStr+="#  Number of Rings in the model\n "+str(nRings) +"\n"
    
    #   The cloud mode and base cloud density
    HeaderStr+="# The cloud mode you are using\n 0 \n"
    HeaderStr+="# The base cloud surface density\n 10. \n"
    #   The tilted ring unit switches
    HeaderStr+="# Central Position Switch (0=degrees, 1=arcsec), Inc/PA Unit switch (0=degrees, 1=arcsec), Velocity Units (0=m/s, 1=km/s), Brightness Units (0=Jy km/s arcsec^-2)\n 0 0 1 0 \n"
    HeaderStr+="# The parameters in each radial bin\n# Rmid    Rwidth    Xcent    Ycent    Inc    PA    VSys    VRot    VRad    Vvert    VDisp    dvdz    Sigma        z0    zGradStart\n"
    
    RWidth=ModelParams[1,0]-ModelParams[0,0]
    print (RWidth)

    EntryStr=""
    for i in range(nRings):
        EntryStr+=str(ModelParams[i,0])+"\t"+str(RWidth) + "\t" \
            +"0. \t 0. \t" + str(ModelParams[i,3])+"\t" + str(ModelParams[i,4])\
            +"\t" + str(ModelParams[i,9]) + "\t" + str(ModelParams[i,1]) \
            +"\t" + str(ModelParams[i,10]) + "\t 0. \t" +str(ModelParams[i,2])\
            +"\t 0. \t" + str(ModelParams[i,6]) +"\t 1. \t 5." \
            +"\n"

    file=open(InputFileName,'w')
    file.write(HeaderStr)
    file.write(EntryStr)
    file.close()

def WriteMCGDataCubeInput(ModelParams,InputFileName)


TestFile="WALLABY_PS_Hya_DR1__Outputs_010_Vel_None/"
TestParams=ReadBarolo(TestFile)

InputFileName="TempTR.in"
WriteMCGTiltedRingInput(TestParams, InputFileName)

InputFileName="TempDC.in"
WriteMCGDataCubeInput(TestParams,InputFileName)
