#!/usr/bin/env python3
from FitDriverScripts import *
import sys as sys
import os as os
import multiprocessing as mp
from multiprocessing import freeze_support


def Main():
    FSGF.GalaxyFit()

if __name__=="__main__":
    freeze_support()
    Main()
