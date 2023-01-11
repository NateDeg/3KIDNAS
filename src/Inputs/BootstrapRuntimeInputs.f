cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting general
c       inputs for the pipeline
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module BootstrapInputMod

      use BootstrapGlobals
      use DataCubeInputMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating bootstrap samples
      subroutine BootstrapIn()
      implicit none
      character(500) TopLevelInputFile
      integer timearray(3)
 

      print*, "Getting the bootstrap generator inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Bootstrap Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the name of the observed cube
      read(10,*)
      read(10,'(A)') ObservedCubeFile
c           Get the name of the model cube file
      read(10,*)
      read(10,'(A)') ModelCubeFile
c           Get the name of the cube mask file
c      read(10,*)
c      read(10,'(A)') MaskCubeFile
c           Get the base name for the output files
      read(10,*)
      read(10,'(A)') BaseOutName

c       Get the size of the blocks in terms of beams and channels
      read(10,*)
c      read(10,*) SpatialBlackSize,VelBlockSize
      read(10,*) VelBlockSize
c       Get the center of the cube for geometry based resampling
      read(10,*)
      read(10,*) BS_Cent%CentX,BS_Cent%CentY,BS_Cent%CentV
     &              ,BS_Cent%PA
      BS_Cent%PA=BS_Cent%PA*Pi/180.

c           Get the base name for the outputs
c      read(10,*)
c      read(10,'(A)') BaseOutName
c           Get the base name of the objects
c      read(10,*)
c      read(10,'(A)') OutputFolder
c       Get the number of bootstrap samples to make
c      read(10,*)
c      read(10,*) nBootstrap
c       Get the noise in mJy
c      read(10,*)
c      read(10,*) RMS
c       Get the random seed
c      read(10,*)
c      read(10,*) idum


      close(10)

c           If the seed is positive, set it by the time array
      if(idum .ge. 0) then
        print*, "using time to generate random seed"
        call itime(timeArray)
        idum=abs(timeArray(1)*idum)
        idum=-int(idum*timeArray(2))-timeArray(3)
      endif


      return
      end subroutine
cccccccc

ccccccc
      subroutine LoadCubesForBootstrap()
      implicit none
      integer MaskSwitch
      Type(Beam2D) TempBeam

      MaskSwitch=0

      print*, "Loading in data cube"
      call ReadFullDataCube(ObservedCube,ObservedBeam
     &                      ,ObservedCubeFile,MaskSwitch)


c      print*, "Loading in model cube"
      call ReadFullDataCube(ModelCube,ObservedBeam
     &                      ,ModelCubeFile,MaskSwitch)
c      call Allocate_Beam2D(ObservedBeam,ObservedCube%DH%nPixels)
c      call DCBrightnessConversion(ObservedCube,ObservedBeam)
c      call DCBrightnessConversion(ModelCube,ObservedBeam)

c      call ReadFullDataCube(MaskCube,TempBeam
c     &                      ,MaskCubeFile,MaskSwitch)
      return
      end subroutine
cccccccc


      end module BootstrapInputMod
