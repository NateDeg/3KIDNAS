import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import argparse
import LoadModel as LM
import CubeAnalysis as CA

def Main():
    XIndx=13
    YIndx=14

    parser = argparse.ArgumentParser(description='Fix the RA & DEC')
    parser.add_argument('BestFitFolder', metavar='Best Fit Folder',help='some string')
    parser.add_argument('ModelBase', metavar='Model Cube',help='some string')
    parser.add_argument('Version', metavar='Version Number',help='some string')
    parser.add_argument('DataCube', metavar='Data Cube',help='some string')
    
    args = parser.parse_args()
    
    ModelFile=args.BestFitFolder+"/"+args.ModelBase+"_AvgModel_v"+args.Version+".txt"

    
    #   Load in the X and Y coordinates
    ModelContents=LoadModelLines(ModelFile)
    X,Y=GetX_Y(ModelContents,XIndx,YIndx)
    #   Load in the DataCube header
    CubeHeader=GetHeader(args.DataCube)
   
    
    RA,DEC=RA_DEC_Calc(X,Y,CubeHeader)
    
    WriteRA_DEC(RA,DEC,ModelFile,XIndx,YIndx,ModelContents)
    
def LoadModelLines(FileName):
    f=open(FileName,"r")
    Lines=f.readlines()
    f.close()
    return Lines
    
    
def GetX_Y(Lines,XIndx,YIndx):
    X=ParseLine_ForVal(Lines[XIndx])
    Y=ParseLine_ForVal(Lines[YIndx])
    return X,Y
    
def ParseLine_ForVal(Line):
    Entries=Line.split()
    Val=float(Entries[2])
    return Val
    
def GetHeader(CubeFile):
    Cube=fits.open(CubeFile)
    Header=Cube[0].header
    return Header

def RA_DEC_Calc(X,Y,CubeHeader):

    print("Correcting RA & DEC")
    CubeWCS = wcs.WCS(CubeHeader)
    #   Set the location of the center in pixels as well as the average model RA and DEC
    Cent=CubeWCS.pixel_to_world(X,Y,0)
    RACent=Cent[0].ra.deg
    DECCent=Cent[0].dec.deg
    
    RACent=round(RACent,7)
    DECCent=round(DECCent,7)
    
    return RACent,DECCent

def WriteRA_DEC(RA,DEC,ModelFile,XIndx,YIndx,Lines):


    RALine="RA_kin (degrees)\t"+str(RA)+"\t0.00\n"
    DECLine="DEC_kin (degrees)\t"+str(DEC)+"\t0.00\n"
    Lines[XIndx+2]=RALine
    Lines[YIndx+2]=DECLine
    
    f=open(ModelFile,"w")
    for line in Lines:
        f.write(line)
    f.close()



Main()
