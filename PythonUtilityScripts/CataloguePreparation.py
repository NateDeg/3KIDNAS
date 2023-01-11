import numpy as np
import argparse

import LoadCatalogue as LC
import CatalogueClass as CC
import SpecialCataloguePrep as SCP

def CataloguePrepInitialize():
    print("Catalogue Prepartion Initialization")

    parser = argparse.ArgumentParser(description='PrepCatalogueOptions')
    #    parser.add_argument('MainFolder', metavar='Name',help='The name of the catalogue')
    #parser.add_argument('BaseName', metavar='Name',help='The name of the catalogue')
    parser.add_argument('CatFile', metavar='Name',help='The name of the catalogue')
    parser.add_argument('OutCatName', metavar='Name',help='The name of the catalogue')
    args = parser.parse_args()
    
    #    MainFolder=args.MainFolder
    CatFile=args.CatFile
    OutCatName=args.OutCatName
    
    #CatFile=MainFolder+"/"+BaseName+"catalog.xls"
    #print(CatFile)

    return CatFile,OutCatName

def CatPrepMain():
    #   Start the main portion
    print("Catalogue Prep Scipt")
    #   Prep the catologue file
    CatFile,OutCatName=CataloguePrepInitialize()
    #   Read the WALLABY SoFiA Catalogue
    WallabyCat=LC.LoadWallabyCat(CatFile)
    #   Select the galaxies for modelling
    ModellingSwitch=[113,74,29]
    #WallabyCat=SCP.LoadSpecCatalogue()
    #ModellingSwitch=[0]

    #   Open up the output catalogue
    file=open(OutCatName,"w")
    #   Write up a catalogue header
    HeaderStr=str(len(ModellingSwitch))+"\n"
    file.write(HeaderStr)

    HeaderStr="ObjectID \t NumberString \t RA \t DEC \t z \t rms \t w20 \t w50 \t"\
        +"ell_maj \t ell_min \t ell_pa \t kin_pa \t freq"+\
        "\n"
    file.write(HeaderStr)
    #Fill in all the catalogue values
    for i in range(len(ModellingSwitch)):
        j=ModellingSwitch[i]
        #    j=33
        print(j,len(WallabyCat.name))
        NumStr=WallabyCat.name[j].split(' ')[1]
        #NumStr=WallabyCat.name[j]

        EntryStr=str(int(WallabyCat.id[j]))+"\t"+NumStr+"\t"+str(WallabyCat.ra[j]) \
                + "\t" + str(WallabyCat.dec[j]) +"\t" +str(WallabyCat.cz[j]) \
                + "\t" + str(WallabyCat.rms[j]) +"\t" + str(WallabyCat.w20[j])\
                + "\t" + str(WallabyCat.w50[j]) +"\t" + str(WallabyCat.ell_maj[j])\
                + "\t" + str(WallabyCat.ell_min[j]) +"\t" + str(WallabyCat.ell_pa[j])\
                + "\t" + str(WallabyCat.kin_pa[j]) +"\t" + str(WallabyCat.freq[j])\
                +"\n"

        file.write(EntryStr)
    file.close()
    

CatPrepMain()


