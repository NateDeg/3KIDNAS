cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the pre-galaxy analysis flow.
c       The goal is to get good initial guesses for the various
c       model parameters.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module InitialAnalysisPreProcessingMod
      use PipelineGlobals
      use GetMomentMapsMod
      use VelProfileAnalysisMod
      use MaskCubeMod
      use DataCubeInputMod
      use EstimateCubeNoiseMod


      implicit none
      contains

ccccc
      subroutine InitialAnalysis_Prep()
      implicit none


      integer TempMaskSwitch
      real RMSTemp

      integer i,j,k
      real FluxLim



      print*, "Doing preanalysis data processing"

      if(PreAnalysisSwitch .eq. 1) then
c       Get the moment maps and profile needed for getting the
c           initial estimates
          if(MomentMapSwitch .eq. 0) then
            print*, "Constructing pre-analysis maps from cube and mask"

            call ConstructProjectionsFromCube(ObservedDC,DataCubeMask)
          elseif(MomentMapSwitch .eq. 1
     &              .or. MomentMapSwitch .eq. 2) then
            print*, "Loading pre-analysis maps from files"
c        call GetMomentMaps()
c        call ReadSoFiAVelocityProfile()
          elseif(MomentMapSwitch .eq. 3) then
            print*, "Constructing map from alternate cube/mask"
c               Load in the secondary cube
            TempMaskSwitch=0
            call ReadFullDataCube(SecondaryObservedDC
     &              ,SecondaryBeam
     &              ,SecondaryCubeFile,TempMaskSwitch)
c       Allocate the beam
            call Allocate_Beam2D(SecondaryBeam
     &          ,SecondaryObservedDC%DH%nPixels)
c           The actual data cube must be in units of Jy/pixel
            call DCBrightnessConversion(SecondaryObservedDC
     &                  ,SecondaryBeam)
c               Load in the mask that we'll use
            TempMaskSwitch=1
            call ReadFullDataCube(DataCubeMask
     &              ,SecondaryBeam
     &              ,SecondaryMaskFile,TempMaskSwitch)
c               TEMPORY FOR TESTING PREANALYSIS CODE
            FluxLim=1.e-4
            do i=0,DataCubeMask%DH%nPixels(0)-1
                do j=0,DataCubeMask%DH%nPixels(1)-1
                    do k=0,DataCubeMask%DH%nChannels-1
                        if(DataCubeMask%Flux(i,j,k)
     &                      .lt. FluxLim) then
                            DataCubeMask%Flux(i,j,k)=0.
                        else
                            DataCubeMask%Flux(i,j,k)=1.
                        endif
                    enddo
                enddo
            enddo
c               In order to get a later profile, it's necessary to get the
c                   cube noise
            call EstimateNoise(SecondaryObservedDC,RMSTemp)
            SecondaryObservedDC%DH%Uncertainty=RMSTemp

c               Use this data cube + mask to get the maps
            call ConstructProjectionsFromCube(SecondaryObservedDC
     &                  ,DataCubeMask)
c           Finally deallocate the secondary mask
            call DeAllocateDataCube(SecondaryObservedDC)
            
          else
            print*, "No acceptable method of "
     &                  ,"getting moment maps selected"
            print*, "Ending program"
            stop
          endif
      endif



      return
      end subroutine
cccccc




ccccc
      subroutine ConstructProjectionsFromCube(ObsDC,MaskDC)
      use DataCubeOutputsMod
      implicit none
      Type(DataCube), INTENT(IN) :: ObsDC,MaskDC
      character(50) Mom0Name

c       First set up a masked cube
      MaskedObservedDC=ObsDC
c       First apply the mask to the cube
      call MaskCube(MaskedObservedDC,MaskDC) !/src/CompareCubes/MaskCube.f
c      print*, "Masked cube flux", sum(MaskedObservedDC%Flux)

c       Obtain the moment maps
      call ConstructMomentMaps(MaskedObservedDC,ObservedMaps)
c           TEMPORARY -- WRITE OUT THE VARIOUS MOM0 Map
c      Mom0Name="TempMom0.fits"
c      call WriteDCSliceToFITS(ObservedMaps,ObservedBeam,Mom0Name,0
c     &                  ,"TestMap")
c       Obtain a velocity profile as well
      call MakeVelProfile(MaskedObservedDC,ObservedVelocityProfile)
    
      return
      end subroutine
ccccccc

      end module
