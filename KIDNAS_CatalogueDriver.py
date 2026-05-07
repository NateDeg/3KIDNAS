"""
KIDNAS_CatalogueDriver.py
=======

**Author:** Nathan Deg

**Description:**
This module is the main catalogue driver for the 3KIDNAS.  It will take a SoFiA catalogue + data products and fit them all.  It will then take examine the fits and produce a catalogue of more robust/reliable models.
"""

#!/usr/bin/env python3
from FitDriverScripts import *
from FullCatalogueDriverScripts import *
import multiprocessing as mp
from multiprocessing import freeze_support


def main():
    """
        This is the main function.  It simply calls the driver in FullCatalogueDriverScripts/CatalogueDrivery.py
    """
    CD.CatalogueDriverMain()

#   This bit is necessary for multiprocessing
if __name__=="__main__":
    freeze_support()
    mp.set_start_method('spawn')
    main()

