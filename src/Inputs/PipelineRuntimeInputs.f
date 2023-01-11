cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting general
c       inputs for the pipeline
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module PipelineInputMod

      use PipelineGlobals
      use SofiaInputMod
      use FittingOptionsMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating a mock cube
      subroutine PipelineIn(PID,Primary)
      implicit none
      character(500) TopLevelInputFile
      integer,INTENT(IN) :: PID,Primary


      print*, "Getting the cube inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Pipeline Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the name of the folder containing everything
      read(10,*)
      read(10,'(A)') MainFolderPath
c           Get the name of the SoFiA catalogue
      read(10,*)
      read(10,'(A)') CubeletFolder

c           Get the name of the SoFiA catalogue
      read(10,*)
      read(10,'(A)') CatalogueFile
c           Get the base name of the objects
      read(10,*)
      read(10,'(A)') ObjectBaseName
c       Get the format of the input SoFiA catalogue (1=ascii, 2= use python script)
      read(10,*)
      read(10,*) SCatalogue%SofiaSwitch
c           If using a python script, get the name of the script from the user.
      if(SCatalogue%SofiaSwitch .eq. 2) then
        read(10,*)
        read(10,'(A)') SCatalogue%PrepScript
      endif

c           Get the name of the fitting options file
      read(10,*)
      read(10,'(A)') FittingOptionsFile

      close(10)

c       Set the catalogue file name to the full name
      CatalogueFile=trim(MainFolderPath)//"/"
     &          //trim(CatalogueFile)
c      print*, trim(CatalogueFile)

c       Name the main output folder  using the object base name
      call MakeOutputFolder(PID,Primary)
c       Also try to make the output folder


c      call GeneralSoFiAIn()

      call FittingOptionsIn()

      return
      end subroutine
cccccccc


cccccc
c
      subroutine MakeOutputFolder(PID,Primary)
      implicit none
      integer,INTENT(IN) :: PID,Primary
      character(500) MakeScript
c       Name the main output folder  using the object base name
      MainOutputFolder=trim(ObjectBaseName)//"Outputs/"
c       If the processor is the primary one, try to make the directory
      if(PID .eq. Primary) then
        MakeScript="mkdir "//trim(MainOutputFolder)
        call system(trim(MakeScript))
      endif

      return
      end subroutine

cccccc

      end module
