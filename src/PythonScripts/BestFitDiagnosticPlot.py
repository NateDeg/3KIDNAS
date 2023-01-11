import numpy as np
from decimal import Decimal
import argparse

import LoadModel as LM
import CubeAnalysis as CA
import DiagnosticPlot as DP

def Main():
    #   Get the set of arguments needed to produce the model cube
    print("Make Diagnostic plot of best fit")
    parser = argparse.ArgumentParser(description='Build a diagnostic Plot')
    parser.add_argument('BestFitFolder', metavar='Best Fit Folder',help='some string')
    parser.add_argument('ModelBase', metavar='Model Cube',help='some string')
    parser.add_argument('Version', metavar='Version Number',help='some string')
    parser.add_argument('DataCube', metavar='Data Cube',help='some string')
    parser.add_argument('MaskCube', metavar='Mask Cube',help='some string')
    parser.add_argument('Distance', metavar='Galaxy Distance',help='some string')
    parser.add_argument('Size', metavar='Galaxy Size',help='some string')
    parser.add_argument('Mass', metavar='logarithmic Galaxy Mass',help='some string')
    

    args = parser.parse_args()
    print(args.BestFitFolder)
    
    ModelFile=args.BestFitFolder+"/"+args.ModelBase+"_AvgModel_v"+args.Version+".txt"
    ModelCubeFile=args.BestFitFolder+"/"+args.ModelBase+"_AverageModel_v"+args.Version+".fits"
    print(ModelFile)
    
    PlotName=args.BestFitFolder+"/"+args.ModelBase+"_AverageModel_v"+args.Version+".png"
    
    GalaxyParams={'Distance':args.Distance, 'Size':args.Size, 'logM':args.Mass,'Name':args.ModelBase,'Version':args.Version,'Folder':args.BestFitFolder}
    
    
    #   Load in the best fitting model
    BestFitModel=LM.LoadBestFitModelFile(ModelFile)
    
    #   Load in the DataCube
    DataCube=CA.BasicCubeAnalysis(args.DataCube)
    #   Load in the MaskCube
    MaskCube=CA.BasicCubeAnalysis(args.MaskCube)
    #   Make sure the mask cube is either 1 or 0
    QuickMaskValCheck(MaskCube)
    #   Now set up the masked data array for moment map calculations
    DataCube['MaskedData']=DataCube['Data']*MaskCube['Data']
    
    #   Load in the Model Cube
    print("Model Cube file", ModelCubeFile)
    ModelCube=CA.BasicCubeAnalysis(ModelCubeFile)
    
    #   Make the diagnostic plot
    DP.MakeDiagnosticPlot(BestFitModel,DataCube,ModelCube,PlotName,GalaxyParams)
   
   
def QuickMaskValCheck(MaskCube):
    FluxLim=1.e-10
    
    #   Set everything below the threshold to 0
    Indics = MaskCube['Data'] < FluxLim
    MaskCube['Data'][Indics] = 0
    #   Set everything above the threshold to 1
    Indics = MaskCube['Data'] > FluxLim
    MaskCube['Data'][Indics] = 1

    #print("Mask sum test",np.sum(MaskCube['Data']))
    return MaskCube
   
Main()

