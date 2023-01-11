ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This code is meant to test routines for fitting a single galaxy
c       using the routines that are built into the WALLABY Kinematic Pipeline
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program SingleGalaxyTester


      use PipelineGlobals
      use SingleFitInputMod
      use DataCubeInputMod
      use PreGalaxyAnalysisMod
      use PostGalaxyAnalysisMod
      use DataCubeOutputsMod

      use FitOutputMod

      use MakeMomentMapsMod

      use CommonConsts

      use ParameterVectorToTiltedRingMod
      implicit none

  
c      integer nProcessors,PID
c      integer Primary
c      parameter(Primary = 0)

      integer MaskSwitch,DataSwitch
      integer i
      character(100) TestMapName

      integer j,k

      real SDConv

      Type(Beam2D) TempBeam

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "Single Galaxy Fit Tests"

      Version=1
c       Get the runtime inputs
      call SingleFitIn(MaskSwitch)
c       Load the cube
      DataSwitch=0
      call ReadFullDataCube(ObservedDC,ObservedBeam
     &          ,DataCubeFile,DataSwitch)
      print*, "Total DC brightness", sum(ObservedDC%Flux)
c       Allocate the beam
      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
c           The actual data cube must be in units of Jy/pixel
      call DCBrightnessConversion(ObservedDC,ObservedBeam)
      print*, "Total DC brightness after unit conversions"
     &              , sum(ObservedDC%Flux)


c       If using a mask, load it in
      if(MaskSwitch .eq. 1) then
c        print*, "Reading Mask File"
        call ReadFullDataCube(DataCubeMask,TempBeam
     &          ,MaskName,MaskSwitch)
c        print*, "Done mask input"
        GalaxyDict%MaskFile=MaskName
c               Set the cube flux to either 1 or zero
        do i=0,DataCubeMask%DH%nPixels(0)-1
            do j=0,DataCubeMask%DH%nPixels(1)-1
                do k=0,DataCubeMask%DH%nChannels-1
                    if(DataCubeMask%Flux(i,j,k) .ge. 1.e-7) then
                        DataCubeMask%Flux(i,j,k)=1.
                    else
                        DataCubeMask%Flux(i,j,k)=0.
                    endif
                enddo
            enddo
        enddo
      elseif(MaskSwitch .eq. 3) then
c           Otherwise set the mask equal to the observed cube
        DataCubeMask=ObservedDC
        DataCubeMask%Flux=1.
        GalaxyDict%MaskFile=DataCubeFile
      endif

      GalaxyFit=>GalaxyFit_Simple
      OutputFit=>OutputBestFit_Simple
c      GalaxyFit=>BBaroloFit

c       Set up a catalogue object to match the general pipeline structure
      call MakeCatalogue()

c       Call the pre-galaxy analysis for the object
      call PreGalaxyAnalysis(SCatLocal%Objects(0))


c      print*, "initial vector", PVIni%Param(0:PVIni%nParams-1)
c      Call Galaxy fit for the test cube
      call GalaxyFit(SCatLocal%Objects(0))
c      print*, "Still initial vector", PVIni%Param(0:PVIni%nParams-1)

c      print*, "Second sanity flux check",sum(ObservedDC%Flux)
c      print*, "Sanity Check part 2", sum(ModelDC%Flux)

      call PostGalaxyAnalysis()


    


c       Output the single galaxy test fit
      call system("mkdir "//trim(MainOutputFolder))
      call OutputFit(SCatLocal%Objects(0))

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccc
      subroutine MakeCatalogue()
      use PipelineGlobals
      use SoFiACatalogueMod
      implicit none

c      print*, "Make catalogue entry"

      SCatLocal%nObjects=1
      call SoFiACatalogueAllocation(SCatLocal)

      call CatalogueEntryDefinition(SCatLocal%Objects(0))

      SCatLocal%SourceName="MCGSuite"

      if(PFlags%CatFlag .eq. 1) then
        open(10,file=trim(CatalogueFile),status='old')
        read(10,*)
        read(10,*) SCatLocal%Objects(0)%EllipseInc
        read(10,*)
        read(10,*) SCatLocal%Objects(0)%EllipsePA
        close(10)
      endif
    

      return
      end subroutine
cccccc

ccccc
      subroutine CatalogueEntryDefinition(CatObj)
      use SoFiACatalogueMod
      use CommonConsts
      use PipelineGlobals
      implicit none
      Type(CatalogueItem), INTENT(INOUT) :: CatObj

c       Set the catalogue object ID
      CatObj%ObID=0
c       Set the catalogue object name using the input name
      CatObj%ObjName=trim(GalaxyDict%GalaxyName)
c       Set the catalogue object RMS using the input RMS
      CatObj%RMS=GalaxyDict%RMS


      return
      end subroutine
cccccc


