cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the pre-galaxy analysis flow.
c       The goal is to get good initial guesses for the various
c       model parameters.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module PreGalaxyAnalysisMod
      use PipelineGlobals
      use InitialAnalysisPreProcessingMod
      use InitialAnalysisMod


      implicit none
      contains

ccccc
      subroutine PreGalaxyAnalysis(CatItem)
      implicit none
      Type(CatalogueItem),INTENT(INOUT) :: CatItem

      real SDConv1,SDTemp
      integer i

c       Before running the profile estimation, the SD limits for the cube
c           must be calculated
      SDConv1=1./(abs((ObservedDC%DH%PixelSize(0)
     &              *ObservedDC%DH%PixelSize(1))))
      do i=1,2
        SDTemp=LinSDLims_MSol(i)
        call MSolPc2_To_JyAS2(SDTemp,SDTemp)
        SDTemp=SDTemp/SDConv1
        SDTemp=SDTemp/abs(ObservedDC%DH%ChannelSize)
        LinSDLims(i)=SDTemp
        LogSDLims(i)=log10(SDTemp)
        print*, "SD limits", i, LinSDLims_MSol(i),LinSDLims(i)
     &              ,LogSDLims(i)
      enddo
    
c       For testing purposes, lock the preanalysis switch to 1
      PreAnalysisSwitch=1

      if(PreAnalysisSwitch .eq. 1) then
c       Prep any data needed for the preanalysis
        call InitialAnalysis_Prep() !/src/PreAnalysis/InitialAnalysis_PreProcessing.f
c       Run the initial analysis to get the quantities needed for the initial parameter estimate
        call InitialAnalysis(CatItem) !/src/PreAnalysis/InitialAnalysis.f
c           Set up the TR fitting options and initial parameter vector needed for the
c               fitting algorithm
c           THe output from InitialAnalysis is a model tilted ring
      elseif(PreAnalysisSwitch .eq. 2) then
c            call TiltRing_DeAllocate(ModelTiltedRing)
            call LoadTRFile()
      endif
       
c       The model tilted ring gives an initiial parameter vector
      call SetUpFittingParameters() !/src/PreAnalysis/ModellingInitializations.f

      return
      end subroutine
cccccc


cccccccc
c
      subroutine LoadTRFile()
      use TiltedRingInputMod
      use UnitConvertMod
      implicit none

      character(100) TempName
      integer i
      real SDConv1,SDTemp


      print*, "Load pre-analysis estimate from file"

      TempName="TestTRFile.txt"
      call ReadTiltedRingModel(ModelTiltedRing,TempName)

      SDConv1=1./(abs((ObservedDC%DH%PixelSize(0)
     &              *ObservedDC%DH%PixelSize(1))))

     
      TR_FittingOptions%nRings=ModelTiltedRing%nRings
      call TiltRingFittingOptions_Allocate(TR_FittingOptions)

      do i=0, ModelTiltedRing%nRings-1
c           Convert to pixels
        ModelTiltedRing%R(i)%Rmid=ModelTiltedRing%R(i)%Rmid
     &              /abs(ObservedDC%DH%PixelSize(0))
        ModelTiltedRing%R(i)%Rwidth=ModelTiltedRing%R(i)%Rwidth
     &              /abs(ObservedDC%DH%PixelSize(0))
c           Convert Surface density
        SDTemp=ModelTiltedRing%R(i)%Sigma
c        print*, "SD_MSol_pc2",SDTemp
c        call JyAS2_To_MSolPc2(SDTemp,SDTemp)
        call MSolPc2_To_JyAS2(SDTemp,SDTemp)
c        print*, "SD_HyAs2",SDTemp,1./SDConv1
        SDTemp=SDTemp/SDConv1
        SDTemp=SDTemp/abs(ObservedDC%DH%ChannelSize)
c        print*, "SD_HyAs2",SDTemp


        ModelTiltedRing%R(i)%Sigma=SDTemp
        ModelTiltedRing%R(i)%logSigma=log10(SDTemp)
c           Convert to radians
        ModelTiltedRing%R(i)%Inclination=
     &          ModelTiltedRing%R(i)%Inclination*Pi/180.
        ModelTiltedRing%R(i)%PositionAngle=
     &          ModelTiltedRing%R(i)%PositionAngle*Pi/180.+Pi/2.


        print*, ModelTiltedRing%R(i)%Rmid,ModelTiltedRing%R(i)%Rwidth
     &          ,ModelTiltedRing%R(i)%Sigma
     &          ,ModelTiltedRing%R(i)%VDisp

      enddo

      ObservedDC%DH%Uncertainty=0.0002

      return
      end subroutine
cccccccc

      end module
