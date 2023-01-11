import sys as sys
import os as os
import copy as copy

def LoadDefaultWRKPFiles(GeneralDict):
    print("Loading in default WRKP control files")
    
    #   Set the default file names using the general dictionary
    Infile=GeneralDict['WRKP_GeneralMainIn']
    OptionsFile=GeneralDict['WRKP_GeneralOptionsIn']
    
    #   Open up the default main input file and read in all the lines and store them
    MainIn=open(Infile,"r")
    Main_In_Lines=MainIn.readlines()
    MainIn.close()
    #   Do the same for the options file
    OptionsIn=open(OptionsFile,"r")
    WRKP_Options_Lines=OptionsIn.readlines()
    OptionsIn.close()
    #   Store all the default line information into the general dictionary and return it
    GeneralDict['MainWRKPInputLines']=Main_In_Lines
    GeneralDict['MainWRKPOptionsLines']=WRKP_Options_Lines
    return GeneralDict
    
def NameWRKPRunFiles(GalaxyDict,BootstrapSwitch):
    if BootstrapSwitch==0:
        GalaxyDict['MainInputFile']=GalaxyDict['ObjNameU']+"_WRKP_Input.txt"
        GalaxyDict['FittingOptionsFile']=GalaxyDict['ObjNameU']+"_WRKP_Options.txt"
    
    return GalaxyDict

def WriteWRKPMainFile(WorkingMainLines,GalaxyDict):

    #   Set the name of the cube to be fit
    WorkingMainLines[1]=GalaxyDict['CubeNameU']+"\n"
    #   Set the name of the fitting options file
    WorkingMainLines[3]=GalaxyDict['FittingOptionsFile']+"\n"
    #   Set the name of the object
    WorkingMainLines[15]=GalaxyDict['ObjNameU']+"\n"
    #   Set the name of the folder that will contain the outputs
    WorkingMainLines[13]=GalaxyDict['TargFolderU']+"\n"
    
    #   Set up the portion of the file to deal with SoFiA shape estimates
    WorkingMainLines[8]=str(1)+"\n"
    if len(WorkingMainLines) <= 18:
        WorkingMainLines.insert(10,GalaxyDict['SoFiAShapeFile']+"\n")
    else:
        WorkingMainLines[10]=GalaxyDict['SoFiAShapeFile']+"\n"

    #   Now that we have the main input file contents written, we can write out the file
    mFile=open(GalaxyDict['MainInputFile'],'w')
    for x in WorkingMainLines:
        mFile.write(x)
    mFile.close()
    
def WriteWRKPOptionsFile(WorkingOptionsLines,GalaxyDict,BootstrapSwitch):
    
    #   Set the name of the cube to be used for initial estimates beyond the inclination and position angle
    WorkingOptionsLines[9]=GalaxyDict['CubeNameU']+"\n"
    #   And give it the name of the mask file
    WorkingOptionsLines[10]=GalaxyDict['MaskNameU']+"\n"
    #   If we are doing a bootstrap, we want to set the number of rings to the original number
    if BootstrapSwitch ==1:
        nRTarg=len(GalaxyDict['BestFitModel']['R'])
        #print("Bootstrap WRKP options", BootstrapSwitch)
        WorkingOptionsLines[28]=str(nRTarg)+"\n"
        #print(nRTarg,WorkingOptionsLines[28])
    
    #   Now write the contents to the temporary options file
    oFile=open(GalaxyDict['FittingOptionsFile'],'w')
    for x in WorkingOptionsLines:
        oFile.write(x)
    oFile.close()



def RunWRKP(GeneralDict,GalaxyDict,BSSwitch):
    print("Running WRKP")
    
    #   First copy the general Main lines to a local set
    WorkingMainLines=copy.copy(GeneralDict['MainWRKPInputLines'])
    #   Next name the temporary WRKP Files for the run
    GalaxyDict=NameWRKPRunFiles(GalaxyDict,0)
    #   Now write up the main input file
    WriteWRKPMainFile(WorkingMainLines,GalaxyDict)
    #   Once this is done, we can get rid of the WorkingMainLines
    del WorkingMainLines
    #   Now we'll do a similar process for the fitting options:
    #       First copy general options lines to a local variable
    WorkingOptionsLines=copy.copy(GeneralDict['MainWRKPOptionsLines'])
    #       And next we'll write the specific fitting options for this run
    WriteWRKPOptionsFile(WorkingOptionsLines,GalaxyDict,BSSwitch)
    #   As before, we'll clean by getting rid of the working set of options lines
    del WorkingOptionsLines
    
    #   Now we can run WRKP
    #       Start with making the runtime command
    WRKPCmd=GeneralDict['FitterExecPath']+ " "+GalaxyDict['MainInputFile']
    #       And now run it
    os.system(WRKPCmd)
    #   And clean up the input files
    ClnCmd="rm "+GalaxyDict['MainInputFile']+" " + GalaxyDict['FittingOptionsFile'] + " " +GalaxyDict['SoFiAShapeFile']
    os.system(ClnCmd)
    #   Finally return the galaxy dictionary
    return GalaxyDict
