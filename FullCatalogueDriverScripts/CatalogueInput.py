import numpy as np
import pandas as pd
import os.path
from os import path

from astropy.io import fits
from astropy.table import Table

"""
    This module contains routines for reading in CSV, Excel, and FITS catalogue files.  It has the routines:
        LoadCSV --> This function loads in a CSV file and puts the catalogue into a dictionary
        LoadXLS --> This function loads in an excel file and puts it into a dictionary
"""

def LoadCatalogue(CatFileName):
    #   First check that the catalogue file exists
    if os.path.isfile(CatFileName)==False:
        print("The given catalogue file is missing:", CatFileName)
        print("Exiting program")
        exit()
    #   Figure out the file type
    FType=CheckCatFileType(CatFileName)
    #   Load the catalogue using the correct file type
    if FType==0:
        Cat=LoadFITSCatalogue(CatFileName)
        TargKeys=['name','team_release']
        Cat=ConvertFITSCatStrings(Cat,TargKeys)
    elif FType==1:
        Cat=LoadCSV(CatFileName,',')
    return Cat

def CheckCatFileType(CatFileName):
    NameString=CatFileName.split(".")
    print(NameString)
    PossibleTypes=['fits','csv']
    CorrectType=False
    for FType in range(len(PossibleTypes)):
        if NameString[-1]==PossibleTypes[FType]:
            CorrectType=True
            break
            
    if CorrectType==False:
        print("Catalogue file is not an acceptable file type")
        print("Acceptable types are:", PossibleTypes)
        print("Exiting program")
        exit()
    return FType
    


def LoadCSV(File,sep):
    """
        This function loads in a CSV file and puts the catalogue into a dictionary
    """
    #   Read the CSV file
    Contents=pd.read_csv(File,sep=sep)
    #   Set the file contents to a dictionary
    Contents = Contents.rename(columns=lambda x: x.strip())
    return Contents

def LoadFITSCatalogue(File):
    #   Read in the data using the astropy table structure
    data= Table.read(File)
    #   Convert the astropy table to a pandas dataframe
    data = data.to_pandas()
    #   Fix the column names by getting rid of whitespace
    data = data.rename(columns=lambda x: x.strip())
    return data

def ConvertFITSCatStrings(Cat,TargKeys):
    for x in TargKeys:
        for i in range(len(Cat[x])):
            Cat[x][i]=Cat[x][i].decode('utf-8')
    return Cat
