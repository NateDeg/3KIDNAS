#!/usr/bin/env python3
from FitDriverScripts import *
from FullCatalogueDriverScripts import *
import sys as sys
import os as os
import multiprocessing as mp
from multiprocessing import freeze_support


def Main():
    #GalFitModules={}
    CD.CatalogueDriverMain()
    



if __name__=="__main__":
    freeze_support()
    Main()

