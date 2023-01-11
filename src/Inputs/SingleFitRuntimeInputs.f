cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting general
c       inputs for fitting a single galaxy
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SingleFitInputMod

      use PipelineGlobals
      use FittingOptionsMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating a mock cube
      subroutine SingleFitIn(MaskSwitch)
      implicit none
      character(500) TopLevelInputFile
      integer MaskSwitch


      print*, "Getting the cube inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Pipeline Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the name of the folder containing everything
      read(10,*)
      read(10,'(A)') DataCubeFile
c           Get the name of the fitting options file
      read(10,*)
      read(10,'(A)') FittingOptionsFile
c           Get the mask switch
      read(10,*)
      read(10,*) MaskSwitch

      read(10,*)
      if(MaskSwitch .eq. 1) then
        read(10,'(A)') MaskName
      endif
c       See whether there is a catalogue file to load in
      read(10,*)
      read(10,*) PFlags%CatFlag
      read(10,*)
      if(PFlags%CatFlag .eq. 1) then
        read(10,'(A)') CatalogueFile
      endif



c       Get a value for the random seed
      read(10,*)
      read(10,*) idum
c       Read in the name of the main output folder to store different fits in
      read(10,*)
      read(10,'(A)') MainOutputFolder

c       Read in the name of the object itself to store different fits in
      read(10,*)
      read(10,'(A)') GalaxyDict%GalaxyName
c       Read in the cube noise (this is temporary until I implement a noise calculation
      read(10,*)
      read(10,*) GalaxyDict%RMS
    

      close(10)


      call FittingOptionsIn()

      GalaxyDict%DataCubeFile=DataCubeFile

      return
      end subroutine
cccccccc

      end module
