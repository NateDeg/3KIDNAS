3KIDNAS Outputs
=================================


There are a variety of data products produced by 3KIDNAS.  When run in pipeline mode, there are a set of products produced for each galaxy regardless of whether the fit is acceptable.  These products are stored in individual folders for each galaxy.  

After fitting all galaxies, 3KIDNAS checks whether the fits are acceptable.  All accepted fits are stored in a subfolder.  Within this subfolder, individual folders are made for each galaxy.  Additionally a catalogue file is generated contains the model parameters for all successful fits. And, for ease in browsing, copies of all diagnostic plots are placed in a separate folder (both for all fits and for only the accepted fits).  Essentially the structure is:

| MainOutputFolder
| ├── AcceptedGalaxyModels
| │   ├── CatalogueFile
| │   ├── DiagnosticPlotFolderForAcceptedFits
| │   ├── AcceptedGalaxy1
| │   └── AcceptedGalaxy2
| │   └── ...
| ├── DiagnosticPlotFolderForAllFits
| ├── README.md
| ├── AttemptedGalaxyFit1
| └── AttemptedGalaxy2
| └── ...    


The AcceptedGalaxyModels folder contains all acceptable models.  Within that folder, the catalogue file is the single most important 3KIDNAS output.


Accepted Model Catalogue
--------------
The accepted model catalogue is the most important output from 3KIDNAS as it is the one that users will likely interact with the most.  It contains the model parameters for each individual galaxy fit.  While this information is available in the individual output folders, it is much more accessible in the output csv catalogue.  The various columns in the catalogue are

.. csv-table:: The 3KIDNAS Catalogue File Columns
   :file: KIDNAS_CatalogueColumnName.csv
   :widths: 30, 70
   :header-rows: 1


Diagnostic Plots
--------------
3KIDNAS produces a diagnostic plot for every galaxy that it models.  An example of such a plot is:

.. image:: WALLABY_J130906+051433_BSModel.png

The upper left panel shows the rotation curve of the model, while the upper right panel shows the surface density profile.  The upper leftmost map shows the moment0 map of the data.  The cyan ellipse is the beam, while the yellow dashed ellipse shows the outermost contour of the model moment0 map.  The upper rightmost map shows the model moment1 map.  The magenta x shows the model centre, while the white arrow points in the direction of the position angle.  The dashed yellow contours are from the model moment1 map, and the solid yellow contour is at the systemic velocity of the system.

The next set of four panels show the major axis PV diagram (left, upper) and minor axis PV diagram (right, upper), and corresponding residual (bottom left and right respectively).  The greyscale in the upper panels are the PV diagrams of the observed maps, while the lower panels show the PV diagrams of the difference cube.  The red contours show the model PV diagrams, and the dotted blue points in the major axis PV diagram show the projected rotation curve.

Finally, the diagram lists the various geometric parameters.  

Accepted Model Outputs
--------------
There are a number of output files generated for any accepted kinematically modelled galaxy.  All the output files follow a naming convention of:
Galaxy Name + Team Release Flag + Output File Type
where the Team Release Flag is the kinTR set in the 3KIDNAS runtime file and the output file type is the specific file.

The first output file is the model output file, which has the output name *_AvgMod.txt.  

