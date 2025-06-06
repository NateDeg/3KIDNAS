cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for
c       comparing the flux contained in 2 different cubes.
c        The cubes must have the same dimensions.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FullModelComparisonMod
      use PipelineGlobals
      use ParameterVectorToTiltedRingMod
      use TiltedRingGenerationMod
      use FillDataCubeWithTiltedRingMod
      use CubeKernelConvolutionMod
      use MaskCubeMod
      use CubeCompareMod


      contains
cccccccc
c
      subroutine TiltedRingModelComparison(TestParams,chi2)
      implicit none
      real,INTENT(IN) :: TestParams(1:PVModel%nParams)
      real,INTENT(OUT) :: chi2
      real pixelarea
      integer i

      logical BadModelFlag


c      print*, "Comparing some parameter vector to the observed data"
c       First reset the model cube
      ModelDC%Flux=0.
c       The model cube should already be allocated ... but the particle arrays will
c           need to be deallocated after running this routine
c       Now put the test params array into the parameter vector notation
c      print*, shape(TestParams)
c      print*, PVModel%nParams


      PVModel%Param=TestParams
c       Set the tilted ring parameters from the vector and the fitting options
      call ParamToTiltedRing(PVModel,ModelTiltedRing
     &          ,TR_FittingOptions)


      if(PFlags%Linear_Log_SDSwitch .eq. 0) then
        ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma=
     &      ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%SigUse
      elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
        ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma=
     &   10.**(ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%SigUse)
      endif
    
c      print*, "Model SD"
c     &          ,ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma

      call BadModelCheck(BadModelFlag)
      if(BadModelFlag) then
        chi2=1.e20
        return
      endif

c      print*, "Sigma profile"
c     &  ,ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma
c      print*, "Siguse profile"
c     &      ,ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%SigUse
c       Finish making the tilted ring model
c      pixelarea=abs(ModelDC%DH%PixelSize(0)*ModelDC%DH%PixelSize(1))



c      CheckForPhysicality()
c      if Unphysical
c        chi2=+1e20


      call BuildTiltedRingModel(ModelTiltedRing,idum)
c       Create the point-source data cube
      call FillDataCubeWithTiltedRing(ModelDC,ModelTiltedRing)
c      print*, "Filled DC", sum(ModelDC%Flux)
c        Convolve the cube with the beam
c           Note that it is assumed that the real beam kernel has already been calculated
      call CubeBeamConvolution(ModelDC,ObservedBeam)
c      print*, "Post Convolve", sum(ModelDC%Flux), sum(ObservedDC%Flux)
c           Now compare the observed and model cubes
      call CubeCompare(ObservedDC,ModelDC,chi2
     &                      ,ObservedDC%DH%Uncertainty)
c      print*, "cube compare", chi2
c           To finish off, deallocate the particle arrays for each ring
c
      do i=0, ModelTiltedRing%nRings-1
        call Ring_ParticleDeAllocation(ModelTiltedRing%R(i))
      enddo

      return
      end subroutine
cccccc


ccccc
c
c       Check that the model doesn't hit any hard, non-physical limits
      subroutine BadModelCheck(BadModelFlag)
      implicit none

      logical, INTENT(INOUT) :: BadModelFlag
      integer i,j
      
      BadModelFlag=.False.

      do i=0, ModelTiltedRing%nRings-1
c       Check that the inclination is between 0 and 90
        if(ModelTiltedRing%R(i)%Inclination .lt. 0.
     &     .or. ModelTiltedRing%R(i)%Inclination.gt. Pi/2.) then
            BadModelFlag=.True.
        endif
c           Check that the rotation is greater than 0
        if(ModelTiltedRing%R(i)%VRot .lt. 0.) then
            BadModelFlag=.True.
        endif
c           Check that the surface density is above 0
        if(ModelTiltedRing%R(i)%Sigma .lt. 0.) then
            BadModelFlag=.True.
        endif
c           Check that the cube center is inside the cube
        do j=0,1
            if(ModelTiltedRing%R(i)%CentPos(j) .lt. 0.
     &      .or. ModelTiltedRing%R(i)%CentPos(j)
     &      .gt. ObservedDC%DH%nPixels(j)-1) then
                BadModelFlag=.True.

            endif
        enddo
      enddo



      return
      end subroutine

ccccc
      end module
