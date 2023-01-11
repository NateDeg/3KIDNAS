cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading a
c         SoFiA catalogue file
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SofiaInputMod

      use PipelineGlobals

      implicit none

      contains

ccccccc
c           This routine selects the specific format of the SoFiA catalogue file
      subroutine GeneralSoFiAIn()
      implicit none

      print*, "Read in SoFiA catalogue"

      if(SCatalogue%SofiaSwitch .eq. 1) then
        call ReadSoFiATextFile()
      elseif(SCatalogue%SofiaSwitch .eq. 2) then
        call PythonCataloguePrep()
      endif

      return
      end subroutine
cccccccc

cccccccc
c           This routine reads in a SoFiA text catalogue file
      subroutine ReadSoFiATextFile()
      implicit none
      logical FileCheck
      integer i, nCols
      real,ALLOCATABLE :: CatItems(:)

      print*, "Reading SoFiA text catalogue"
c           Check that the catalogue file exists
      INQUIRE(file=trim(CatalogueFile),EXIST=FileCheck)
      if(FileCheck .eqv. .False.) then
        print*, "Missing catalogue file "
     &          , trim(CatalogueFile)
        stop
      endif
c       Now figure out how many items are in the catalogue
      open(10,file=trim(CatalogueFile),status='old')
      i=0
      do
        read(10,*,END=100)
        i=i+1
      enddo
100   continue
c           Allocate the catalogue array
      SCatalogue%nObjects=i-4
      print*, "number of objects",SCatalogue%nObjects
      call SoFiACatalogueAllocation(SCatalogue)
c           Go back to the beginning of the file
      rewind(10)
c       Skip past the first four lines
      do i=1,4
        read(10,*)
      enddo
c       Finally read the catalogue
c           The catalogue will have more columns than we need, so first set up
c           an array to store those items temporarily
      nCols=33
      ALLOCATE(CatItems(nCols))
      do i=0,SCatalogue%nObjects-1
        read(10,*) CatItems(1:nCols)
c           Get the relevant catalogue parameters
        SCatalogue%Objects(i)%ObID=CatItems(1)
        SCatalogue%Objects(i)%RA=CatItems(6)
        SCatalogue%Objects(i)%DEC=CatItems(7)
        SCatalogue%Objects(i)%ZCent=CatItems(8)
        SCatalogue%Objects(i)%XMin=CatItems(9)
        SCatalogue%Objects(i)%XMax=CatItems(10)
        SCatalogue%Objects(i)%YMin=CatItems(11)
        SCatalogue%Objects(i)%YMax=CatItems(12)
        SCatalogue%Objects(i)%ZMin=CatItems(13)
        SCatalogue%Objects(i)%ZMax=CatItems(14)
        SCatalogue%Objects(i)%EllipseMaj=CatItems(18)
        SCatalogue%Objects(i)%EllipseMin=CatItems(19)
        SCatalogue%Objects(i)%EllipsePA=CatItems(20)
        SCatalogue%Objects(i)%RMS=CatItems(29)
        SCatalogue%Objects(i)%nChannels=int(CatItems(16))
        SCatalogue%Objects(i)%n_LOS=int(CatItems(17))
        print*,SCatalogue%Objects(i)%EllipseMaj
     &          ,SCatalogue%Objects(i)%EllipseMin
     &          ,SCatalogue%Objects(i)%EllipsePA
     &          ,SCatalogue%Objects(i)%RMS

      enddo
c       And close the file
      close(10)

      return
      end subroutine
ccccccccccc



cccccc
c       This routine runs a python script to convert the input catalogue to
c           the format needed for the pipeline
      subroutine PythonCataloguePrep()
      implicit none
      character(1000) script
      print*, "Prepping catalogue using ", trim(SCatalogue%PrepScript)


      OutCatName="KinematicCatalogue.txt"

      script="python "//trim(SCatalogue%PrepScript)//" "
     &              //trim(CatalogueFile)// " "
c     &              //trim(MainFolderPath)// " "
c     &              //trim(BaseName) //" "
     &              //trim(OutCatName)
      call system(trim(script))

      call PythonCatalogueRead()

      return
      end subroutine
cccccc


cccccc
c       This routine reads in a text catalogue prepped by the python script
c
      subroutine PythonCatalogueRead()
      implicit none
      integer i
      integer ID
      character(50) NumString
      real Items(11)

      Items=0.
      open(10,file=trim(OutCatName),status='old')
      read(10,*) SCatalogue%nObjects

      call SoFiACatalogueAllocation(SCatalogue)
      read(10,*)

      do i=0, SCatalogue%nObjects-1
        read(10,*) ID, NumString,Items

        SCatalogue%Objects(i)%ObID=ID
        SCatalogue%Objects(i)%ObjName=trim(NumString)
        SCatalogue%Objects(i)%RA=Items(1)
        SCatalogue%Objects(i)%DEC=Items(2)
        SCatalogue%Objects(i)%ZCent=Items(3)
        SCatalogue%Objects(i)%RMS=Items(4)
        SCatalogue%Objects(i)%w20=Items(5)
        SCatalogue%Objects(i)%w50=Items(6)
        SCatalogue%Objects(i)%EllipseMaj=Items(7)
        SCatalogue%Objects(i)%EllipseMin=Items(8)
        SCatalogue%Objects(i)%EllipsePA=Items(9)
        SCatalogue%Objects(i)%kinPA=Items(10)
        SCatalogue%Objects(i)%CentralFreq=Items(11)
c        print*, "Catalogue item", i, ID,trim(NumString)

      enddo

      close(10)

      return
      end subroutine
ccccccccc


      end module
