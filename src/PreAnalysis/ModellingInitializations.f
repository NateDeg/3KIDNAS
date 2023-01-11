cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the pre-galaxy analysis flow.
c       The goal is to get good initial guesses for the various
c       model parameters.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ModellingInitializationsMod
      use PipelineGlobals
      use ParameterVectorToTiltedRingMod
      use CommonConsts
      use InputUnitConversionsMod

      implicit none
      contains

ccccc


ccccc
c           This routine sets up the modelling parameters needed for the fitting routine
      subroutine SetUpFittingParameters()
      implicit none

c       The method of setting up the fitting parameters will depend on the
c           specific fitting algorithm.
      if (PFlags%CoreCodeSwitch .eq. 1) then
c           When using the inbuilt downhill simplex code (option 1), it is
c           necessary to convert a tilted ring model into a parameter vector,
c           and set the limits on said vector.
c               First, set up the tilted ring fitting parameters
        call SetupTRFittingOptionsFromModelTR(
     &              PFlags%Linear_Log_SDSwitch)
c                Now set the various parameter limits
        call SetupTRParamLimitsFromModelTR(
     &          PFlags%Linear_Log_SDSwitch
     &          ,2*ObservedBeam%BeamMajorAxis)
c               Finally set the initial parameter vector
        call SetInitialParamVector()
      endif
    
      return
      end subroutine

    



cccccc
c
      subroutine SetupTiltedRingStructures(nRings)
      implicit none
      integer i
      integer,INTENT(IN):: nRings
      real RWidth
      real ApproxSize
      
c       Use the number of rings to set the size flag
      ApproxSize=nRings/TR_FittingOptions%nRingsPerBeam
      if(ApproxSize .lt. SizeLims(1)) then
        GalaxyDict%Flags%SizeFlag=1
      elseif(ApproxSize .gt. SizeLims(2)) then
        GalaxyDict%Flags%SizeFlag=2
      else
        GalaxyDict%Flags%SizeFlag=0
      endif
c       Set the number of rings in both the Tilted Ring and Tilted Ring fitting options structures
      ModelTiltedRing%nRings=nRings
      TR_FittingOptions%nRings=nRings
c       Allocate the tilted ring and tilted ring fitting options structures
      call TiltRing_Allocate(ModelTiltedRing)
      call TiltRingFittingOptions_Allocate(TR_FittingOptions)
c       Calculate the width of each ring
      RWidth=ObservedBeam%BeamMajorAxis
     &              /TR_FittingOptions%nRingsPerBeam
c       Set up the radius and midpoint for the tilted ring
      do i=0,TR_FittingOptions%nRings-1
        ModelTiltedRing%R(i)%Rwidth=RWidth
      enddo

      return
      end subroutine
ccccccccc


cccccc
c
      subroutine SetupTRFittingOptionsFromModelTR(SDSwitch)
      implicit none
      integer,INTENT(IN) :: SDSwitch
      integer i


c           Simply set the radial profiles of the TR_FittingOptions to
c               the model TR objects
      do i=0,TR_FittingOptions%nRings-1
        TR_FittingOptions%RadialProfiles(i)=ModelTiltedRing%R(i)
c        print*, "Fitting options Rmid"
c     &          , TR_FittingOptions%RadialProfiles(i)%Rmid
c     &          , TR_FittingOptions%RadialProfiles(i)%VDisp
c           Depending on the SD switch, set SigUse
        if(SDSwitch .eq. 0) then
            TR_FittingOptions%RadialProfiles(i)%SigUse=
     &              TR_FittingOptions%RadialProfiles(i)%Sigma
        elseif(SDSwitch .eq. 1) then
            TR_FittingOptions%RadialProfiles(i)%SigUse=
     &              TR_FittingOptions%RadialProfiles(i)%logSigma
        endif
        print*, "Initial Paramvector"
     &      ,TR_FittingOptions%RadialProfiles(i)%Sigma
     &      ,TR_FittingOptions%RadialProfiles(i)%logSigma
      enddo

c       Now use the previously set fitting options to set up the parameter indexing
c           This indexing also gets the number of free parameters
      call LogicalTiltedRingIndexing(TR_FittingOptions) !/ObjectDefinitions/TiltedRing.f


      print*, "Number of free parameters"
     &          ,TR_FittingOptions%nFittedParamsTotal
      return
      end subroutine
ccccc

ccccc
c
      subroutine SetupTRParamLimitsFromModelTR(SDSwitch,CentSize)
      implicit none
      integer,INTENT(IN) :: SDSwitch
      real,INTENT(IN) :: CentSize
      integer i

      real AngLowLims(2),AngHighLims(2)

      real SDTemp,SDConv1


      TR_FittingOptions%CyclicSwitch(0:12)=0
c           Set the limits on the center to be within +/- 2 beams

      do i=0,1
        TR_FittingOptions%ParamLowerLims(i)
     &              =ModelTiltedRing%R(0)%CentPos(i)
     &              -CentSize
        TR_FittingOptions%ParamUpperLims(i)
     &              =ModelTiltedRing%R(0)%CentPos(i)
     &              +CentSize

        TR_FittingOptions%ParamRange(i)=5.        !   Instead of limits, try ranges
c        print*, "hmmm", i,ModelTiltedRing%R(0)%CentPos(i),CentSize
      enddo
c       Set the limits on the inclination
      TR_FittingOptions%ParamLowerLims(2)=0.
      TR_FittingOptions%ParamUpperLims(2)=Pi/2.

      TR_FittingOptions%ParamRange(2)=10.*Pi/180.
c      AngLowLims(1)=0.
c      AngHighLims(1)=Pi/2.
c      AngLowLims(2)=ModelTiltedRing%R(0)%Inclination
c     &              -10.*Pi/180.
c      AngHighLims(2)=ModelTiltedRing%R(0)%Inclination
c     &              +10.*Pi/180.
c      TR_FittingOptions%ParamLowerLims(2)=maxval(AngLowLims)
c      TR_FittingOptions%ParamUpperLims(2)=minval(AngHighLims)
c       Set the limits on the position angle
      TR_FittingOptions%ParamLowerLims(3)=0.
      TR_FittingOptions%ParamUpperLims(3)=2.*Pi

      TR_FittingOptions%ParamRange(3)=10.*Pi/180.
c      AngLowLims(1)=0.
c      AngHighLims(1)=2.*Pi
c      AngLowLims(2)=ModelTiltedRing%R(0)%PositionAngle
c     &              -90.*Pi/180.
c      AngHighLims(2)=ModelTiltedRing%R(0)%PositionAngle
c     &              +90.*Pi/180.
c      TR_FittingOptions%ParamLowerLims(3)=maxval(AngLowLims)
c      TR_FittingOptions%ParamUpperLims(3)=minval(AngHighLims)


      TR_FittingOptions%CyclicSwitch(3)=1   !   The PA is cyclic
c       Set the limits on VSys using the profile edges
c           Use the minval and maxvel as velocity profiles may be increasing or decreasing
      TR_FittingOptions%ParamLowerLims(4)
     &              =minval(ObservedDC%Channels)
      TR_FittingOptions%ParamUpperLims(4)
     &              =maxval(ObservedDC%Channels)


      TR_FittingOptions%ParamRange(4)=20.
c       Set the general limits on VRot in km/s
      TR_FittingOptions%ParamLowerLims(5)=VRotLims(1)
      TR_FittingOptions%ParamUpperLims(5)=VRotLims(2)

      TR_FittingOptions%ParamRange(5)=50.
c       Set the limits on VRad in km/s
      TR_FittingOptions%ParamLowerLims(6)=-50.
      TR_FittingOptions%ParamUpperLims(6)=50.

      TR_FittingOptions%ParamRange(6)=50.
c       Set the limits on VDisp in km/s
      TR_FittingOptions%ParamLowerLims(7)=0.
      TR_FittingOptions%ParamUpperLims(7)=20.

      TR_FittingOptions%ParamRange(7)=5.
c       Set the limits on Vvert in km/s
      TR_FittingOptions%ParamLowerLims(8)=-20.
      TR_FittingOptions%ParamUpperLims(8)=20.

      TR_FittingOptions%ParamRange(8)=5.
c       Set the limits on dvdz
      TR_FittingOptions%ParamLowerLims(9)=0.
      TR_FittingOptions%ParamUpperLims(9)=5.

      TR_FittingOptions%ParamRange(9)=1.
c       Set the limits on Sigma         -- This needs work
      if(SDSwitch .eq. 0) then
        TR_FittingOptions%ParamLowerLims(10)=LinSDLims(1)
        TR_FittingOptions%ParamUpperLims(10)=LinSDLims(2)

        TR_FittingOptions%ParamRange(10)=0.002

        SDTemp=2.
        SDConv1=1./(abs((ObservedDC%DH%PixelSize(0)
     &              *ObservedDC%DH%PixelSize(1))))
        call MSolPc2_To_JyAS2(SDTemp,SDTemp)
        SDTemp=SDTemp/SDConv1
        SDTemp=SDTemp/abs(ObservedDC%DH%ChannelSize)
        TR_FittingOptions%ParamRange(10)=SDTemp
c        print*, "New SD Initialization Range"
c     &          ,TR_FittingOptions%ParamRange(10)



      elseif(SDSwitch .eq. 1) then
        TR_FittingOptions%ParamLowerLims(10)=LogSDLims(1)
        TR_FittingOptions%ParamUpperLims(10)=LogSDLims(2)

        TR_FittingOptions%ParamRange(10)=1.
      endif
c       Set the limits on z0
      TR_FittingOptions%ParamLowerLims(11)=0.
      TR_FittingOptions%ParamUpperLims(11)=5.

      TR_FittingOptions%ParamRange(11)=1.
c       Set the limits on zGradiantStart
      TR_FittingOptions%ParamLowerLims(12)=0.5
      TR_FittingOptions%ParamUpperLims(12)=10.

      TR_FittingOptions%ParamRange(12)=1.

      return
      end subroutine
cccccc


cccccc
c
      subroutine SetInitialParamVector()
      implicit none
c       First allocate the parameter vector
      PVIni%nParams=TR_FittingOptions%nFittedParamsTotal
      call AllocateParamVector(PVIni)
      call TiltedRingOptionsToPV(PVIni,TR_FittingOptions)
      return
      end subroutine
cccccccc



      end module
