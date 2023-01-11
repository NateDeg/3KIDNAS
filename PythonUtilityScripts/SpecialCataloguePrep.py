import numpy as np
import argparse
import CatalogueClass as CC


def LoadSpecCatalogue():
    nObj=0
    Cat=CC.Catalogue(1)
    Cat.id=[0]
    Cat.ra=[0.]
    Cat.dec=[0.]
    Cat.z=[-1.]
    Cat.rms=[1.6]
    Cat.w20=[100.]
    Cat.w50=[100.]
    Cat.ell_maj=[105.0]
    Cat.ell_min=[74.24]
    Cat.ell_pa=[0.]
    Cat.kin_pa=[0.]
    Cat.freq=[-1]
    Cat.cz=[0.]
    Cat.name=["ba_5.5.mass_8.inc_45.0.pa_0.0.veldisp_8.0"]
    

    return Cat
