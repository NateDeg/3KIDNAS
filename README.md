# 3KIDNAS
 The full pipeline for the WALLABY Kinematics analysis


---
Requirements

There are a number of standard python packages used in 3KIDNAS -- numpy, scipy, astropy, matplotlib, multiprocessing.  
Additionally, the colourspace and CosmosCanvas packages are stored in third_party/

The code currently is fixed to use python3.9 (this will be adjusted in the future).



---
Quick Installation Guide

The first step is to compile the key 3rd party software -- the easiest approach is to go to terminal then

1) cd third_party/fftw-3.3.8/
    make clean
    ./configure
    make
2) cd third_party/cfitsio/
    make clean
    ./configure
    make
3) cd third_party/SoFiA-2-master_2_5_1/
    make clean
    ./compile.sh
    
Now it should be possible to compile the main code

1) cd src/
    make clean
    make
    
---
Running and testing the code

There is a set of test data and sample inputs currently available at:


To run on a single galaxy go to the folder with the GalaxyInputFile.py and run:
python $(PATH)/WRKP_GalaxyFitDriver.py $(GalaxyInputFile)

To run on a catalogue:
python $(PATH)/WRKP_CatalogueDriver.py $(CatalogueInputFile)

---
Known Bugs
1) There is a known issue with linux vs Mac architechtures.  To adjust the code for linux:
    a) Change line 161 in src/Inputs/DataCubeInput.f from "integer\*8 nInts(3)" to "integer nInts(3)"
    b) Change line 118 in src/Inputs/DataCubeOutputs.f from "integer\*8 naxesT(3),naxisT" to "integer naxesT(3),naxisT"
    
2) There is a potential incompatibility with the third_party/colourspace package and the current version of matplotlib.


