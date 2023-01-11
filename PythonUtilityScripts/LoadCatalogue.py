import numpy as np
import pandas as pd

import CatalogueClass as CC

def LoadWallabyCat(fname):
    print("Loading Catalogue from", fname)
    xls = pd.read_excel(fname,sheet_name="Catalog",header=16)
    Catalogue=CC.Catalogue(len(xls['name'])-2)
    for i in range(len(xls.columns)):
        vars(Catalogue)[xls.columns[i]]=xls[xls.columns[i]][2:2+Catalogue.nGalaxies].values
    return Catalogue
