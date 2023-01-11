cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting the
c       datacube objects (cube, mask, and possibly maps)
c       needed for analyzing a particular galaxy
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FullGalaxyInputMod
      use PipelineGlobals
      use DataCubeInputMod
      use MaskCubeMod


      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating a mock cube
      subroutine FullGalaxyIn(SoFiAObjectID)
      implicit none
      integer,INTENT(IN) :: SoFiAObjectID
      integer j
      Type(Beam2D) TempBeam

      character(10) CatalogueName
      integer MaskSwitch

c      print*, "Load in full galaxy inputs"
c     &          ,SCatLocal%Objects(SoFiAObjectID)%ObID

      call GetFileNames(SoFiAObjectID)      !This file
c       Build local velocity cubes
      call BuildVelCubes(SoFiAObjectID)     !This file

c       Load the datacube
      MaskSwitch=0
      print*, "DC File ", trim(DataCubeFile)
c      call ReadFullDataCube(ObservedDC,ObservedBeam,trim(DataCubeFile))
      call ReadFullDataCube(ObservedDC,ObservedBeam
     &                  ,trim(TempCubeFile),MaskSwitch)
      print*, "Gal Input Beam", SoFiAObjectID,ObservedBeam%BeamMajorAxis
c           The actual data cube must be in units of Jy/pixel
      call DCBrightnessConversion(ObservedDC,ObservedBeam)


c       Allocate the Beam now for simplicity
      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
      print*, "Beam Allocated"

c       Read in the mask
c      call ReadFullDataCube(DataCubeMask,TempBeam,trim(MaskName))
      MaskSwitch=1
      call ReadFullDataCube(DataCubeMask,TempBeam
     &          ,trim(TempMaskFile),MaskSwitch)



      return
      end subroutine
cccccccc

ccccccccc
c           This routine names the various data files that will be used
      subroutine GetFileNames(SoFiAObjectID)
      implicit none
      integer,INTENT(IN) :: SoFiAObjectID
c      character(1000) BaseFileName

      integer j

      j=SoFiAObjectID
c      print*, j
      BaseName=
     &              trim(MainFolderPath)//"/"
     &              //trim(CubeletFolder)//"/"
     &              //trim(ObjectBaseName)
     &              //trim(SCatLocal%Objects(j)%ObjName)
c      print*, "hmmm ",trim(ObjectBaseName)
c      print*, "hmmm2 ",trim(SCatLocal%Objects(j)%ObjName)
c      print*, "Base FName",j,trim(BaseName)

      DataCubeFile=trim(BaseName)//"_cube.fits"
      MaskName=trim(BaseName)//"_mask.fits"
      MapNames(0)=trim(BaseName)//"_moment0.fits"
      MapNames(1)=trim(BaseName)//"_moment1.fits"
      MapNames(2)=trim(BaseName)//"_moment2.fits"

c      print*, "Datacube file", trim(DataCubeFile)


      return
      end subroutine
ccccccc


cccccc
c
      subroutine BuildVelCubes(SoFiAObjectID)
      implicit none
      integer,INTENT(IN) :: SoFiAObjectID

      real RA,DEC

      character(1000) ConversionScript
      character(20) RA_Str, Dec_Str
      integer j
      j=SoFiAObjectID


      RA=SCatLocal%Objects(j)%RA
      DEC=SCatLocal%Objects(j)%DEC
      write(RA_Str,'(F10.5)') RA
      write(DEC_Str,'(F10.5)') DEC
      ConversionScript="python PythonUtilityScripts/CubePreparation.py "
     &          //trim(BaseName)//"_"
     &          //" "//trim(RA_Str)//" "//trim(Dec_Str)

c      print*, "Sanity checks", RA, DEC
c      print*, "Str checks ", RA_Str, " ", DEC_Str
c      print*, "Conversion Script ", trim(ConversionScript)
      call system(ConversionScript)

      TempCubeFile="TempDataCube.fits"
      TempMaskFile="TempMaskCube.fits"
      TempMapNames(0)="TempMom0.fits"
      TempMapNames(1)="TempMom1.fits"
      TempMapNames(2)="TempMom2.fits"


      return
      end subroutine
cccccccc

      end module
