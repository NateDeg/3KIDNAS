"""
KIDNAS_GalaxyFitDriver.py
=======

**Author:** Nathan Deg

**Description:**
This module is the main driver for fitting a single galaxy with 3KIDNAS.  It takes in a cube + mask and some initial estimates and fits the cube.  It also determines the uncertainties using the flipping bootstrap and will attempt to extract scaling parameters.
"""

#!/usr/bin/env python3
from FitDriverScripts import *
import multiprocessing as mp
from multiprocessing import freeze_support


def main():
    """
        This is the main function.  It simply calls the driver in FitDriverScripts/FullSingleGalaxyFit.py
    """
    FSGF.GalaxyFit()

#   This bit is necessary for multiprocessing
if __name__=="__main__":
    freeze_support()
    main()
