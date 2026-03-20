cccccc
c       This module contains the global variable definitions for the
c           bootstrap generator program
ccccccccccc

      module BootstrapGlobals
c               Use the various object definition routines
      use DataCubeMod
      use BeamMod
      use ParticleMod

      implicit none

c           This particular variable stores the center
c               used for bootstrap resampling routines.  They should be
c               in units of pixels and radians
      Type BootstrapCenter
        real CentX
        real CentY
        real CentV
        real PA
        real Inc
      end Type


      Type(DataCube) ObservedCube, ModelCube,DifferenceCube
      Type(DataCube) BootstrapCube, MaskCube

      Type(BootstrapCenter) BS_Cent


      Type(Beam2D) ObservedBeam

      integer idum, nBootstrap

      real RMS

      character(600) ObservedCubeFile,ModelCubeFile
      character(600) MaskCubeFile
      character(600) OutputFolder
      character(600) BaseOutName
      character(600) DifferenceCubeName
      character(600) SampleFile

      real SpatialBlackSize,VelBlockSize

      end module BootstrapGlobals
ccccc
