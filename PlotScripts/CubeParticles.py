#!/usr/bin/env python3
import numpy as np

from decimal import Decimal


import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator
from mpl_toolkits.mplot3d import Axes3D

   

   
def Main():
 
    File="Asym_DifferenceCube_Particles.txt"
    X,Y,VRot,Flux=np.loadtxt(File,skiprows=0,unpack=True)
    
    
    left=0.15
    base=0.2
    w=0.7
    h=w
    wBuf=0.02

    fig=plt.figure(figsize=(5,5))
    ax=fig.add_axes([left,base,w,h],projection='3d')

    ax.scatter(X,Y,VRot,alpha=0.5,marker='.',c=Flux,cmap='plasma')

    #plt.show()
    plt.savefig("CubeParts.png", format='png', bbox_inches='tight')
    plt.close()


Main()
