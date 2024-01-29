Basic Usage
=================================

Using 3KIDNAS
-------------

3KIDNAS is meant to be run in pipeline mode, but it is possible to run the various components on their own.

Pipeline Mode
^^^^^^^^^^^^^^^^^

Running 3KIDNAS in pipeline mode has been designed to be fairly straightforward.  It requires a user provided runtime input file that will point to the data; a SoFiA-2 catalogue and a set of cubelets and associated masks.  With those present, then the pipeline can be run using

.. code-block:: python ~/$PathToKIDNAS/WRKP_CatalogueDriver.py $PipelineInputFile

Single Galaxy Mode
^^^^^^^^^^^^^^^^^

It is also possible to fit a single galaxy rather than a full suite of them.  It also requires an input file that specifies some runtime options.  With that file in place a fit to a single galaxy can be run using the command

.. code-block:: python ~/$PathToKIDNAS/WRKP_GalaxyDriver.py $SingleGalaxyInputFile

BootStrapSampler
^^^^^^^^^^^^^^^^^

The bootstrap resampler code is designed to generate a single bootstrap resampled cubelet.  It can used with a data cubelet, a model cubelet, and the model geometry.  With the input text file set up, the code can be run via

.. code-block:: ~/$PathToKIDNAS/Programs/BootStrapSampler $BootstrapSamplerInputFile

SingleGalaxyFitter
^^^^^^^^^^^^^^^^^

The program at the heart of KIDNAS is SingleGalaxyFitter.  This program finds the best fitting model for some cubelet.  As with all the other programs, there are a set of input files that configure the run.  Once set, the program is run via

.. code-block:: ~/$PathToKIDNAS/Programs/SingleGalaxyFitter $SingleGalaxyFitterInFile





