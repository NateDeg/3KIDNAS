cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the pre-galaxy analysis flow.
c       The goal is to get good initial guesses for the various
c       model parameters.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module InitialAnalysisMod
      use PipelineGlobals
      use MaskCubeMod

      use EstimateCubeNoiseMod
      use EstimateShapeMod
      use EstimateProfileMod
      use VelProfileAnalysisMod
      use ParameterVectorToTiltedRingMod

      use SoFiACatalogueMod

      use CommonConsts
      use DataCubeInputMod

      use InputUnitConversionsMod

      use ModellingInitializationsMod

      use SofiaProfileInputMod


      implicit none
      contains

ccccc
c
c       The initial analysis must generate a set of tilted ring parameters and store
c           them in the ModelTiltedRing variable.
c
      subroutine InitialAnalysis(CatItem)
      implicit none
      Type(CatalogueItem),INTENT(INOUT) :: CatItem

      real Incl,PA, RMax,VSys, Noise
      real Center(0:1),VEdges(0:1), TestPos(2)
      real,ALLOCATABLE :: EstimatedProfiles(:,:)
      real,ALLOCATABLE :: EstimatedRadialProfiles(:,:)
      integer nRings,i



      print*, "Running the initial analysis to get"
     &          ," the initial parameter estimate(s)"
c       It is always necessary to get the noise
c       Get the noise
      call GetNoise(Noise,CatItem)
c       Set the observed data cube noise value to the noise value
      ObservedDC%DH%Uncertainty=Noise

c       If the preanalysis switch == 1, use the moment maps and velocity profile
c           to get the initial model parameters
      if (PreAnalysisSwitch.eq.1) then
        call ProjectedAnalysis(CatItem)
      endif


      return
      end subroutine
cccccc

cccccccc
c           This subroutine uses the moment maps/velocity profile and/or
c               catalogue values to get the initial parameter estimates
      subroutine ProjectedAnalysis(CatItem)
      implicit none
      Type(CatalogueItem),INTENT(INOUT) :: CatItem

      real Incl,PA, RMax,VSys, Noise,VDisp
      real Center(0:1),VEdges(0:1),minV,maxV
      real,ALLOCATABLE :: EstimatedProfiles(:,:)
      real,ALLOCATABLE :: EstimatedRadialProfiles(:,:)
      integer ProfShape(2)
      integer nRings

      real SDTemp,SDConv1
      integer i

c       Find the center position of the object
      call GetObjectCenter(Center,CatItem)

c           Get the inclination and position angle estimates
      call GetGalaxyShape(Incl,PA,Center,CatItem)

c           Estimate Vsys from the observed velocity profile

      ProfShape=shape(ObservedVelocityProfile)
      call EstimateVSysAndEdges(ProfShape(2)
     &                  ,ObservedVelocityProfile
     &                  ,0.3,VSys,VEdges)   !/src/PreAnalysis/VelProfileAnalysis.f
c       Check that VSys lives between the limits
      minV=minval(ObservedVelocityProfile(0,:))
      maxV=maxval(ObservedVelocityProfile(0,:))
      if ((VSys .lt. minV) .or. (VSys .gt. maxV)) then
        VSys=(minV+maxV)/2.
      endif


c           Estimate the rotation curve and surface density profile
c               use the model maps, center, and shape
c               It is important that the observed maps have an appropriate
c               uncertainty measure based on the cube that is used to produce the
c               maps.  When using secondary cubes/masks, this is done in the preprocessing step
      if(MomentMapSwitch .eq. 0) then
        ObservedMaps%DH%Uncertainty=ObservedDC%DH%Uncertainty
      endif


      call EstimateProfiles(ObservedMaps,ObservedBeam,Center
     &                  ,Incl,PA,RMax
     &                  ,EstimatedProfiles,EstimatedRadialProfiles
     &                  ,VSys,TR_FittingOptions,nRings
     &                  ,LinSDLims,VRotLims,NoiseSigmaLim
     &                  )!/src/PreAnalysis/EstimateRadialProfiles.f
c       Get an estimate of the initial velocity dispersion
      call EstimateFlatVDisp(VDisp)
c       Store the initial parameters into a tilted ring structure
      call ConvertFlatDiskProfilesToTR(nRings,Center,Incl,PA
     &              ,VSys,VDisp,EstimatedRadialProfiles)
c       Finally clean the memory and deallocate the profiles
      DEALLOCATE(EstimatedProfiles)
      DEALLOCATE(EstimatedRadialProfiles)
      DEALLOCATE(ObservedVelocityProfile)

      return
      end subroutine
cccccc



ccccccc
c       This routine acts as an interface to the various
c           routines for getting an initial estimate
c           of the galaxy center
      subroutine GetObjectCenter(Center,CatItem)
      implicit none
      real,INTENT(INOUT)  :: Center(0:1)
      Type(CatalogueItem),INTENT(IN) :: CatItem
      integer Flag

      GalaxyDict%Flags%CenterFlag=0

      if(PFlags%CenterSource .eq. 0) then
        call Iter_EstimateCenter(ObservedMaps,Center,Flag) !/src/PreAnalysis/EstimateShape.f
        GalaxyDict%Flags%CenterFlag=Flag
c        print*, "Iter center estimate", Center
      elseif(PFlags%CenterSource .eq. 1) then
c        Center(0)=CatItem%XCent
c        Center(1)=CatItem%YCent
        Center(0)=ObservedDC%DH%RefLocation(0)
        Center(1)=ObservedDC%DH%RefLocation(1)
      elseif(PFlags%CenterSource .eq. 2) then
        Center=ObservedMaps%DH%PixelCent
      endif

c      Center(0)=77
c      Center(1)=85
      print*, "Initial Center Guess",PFlags%CenterSource
     &          ,Center

      return
      end subroutine
cccccccc


ccccccc
c       This routine acts as an interface to the various
c           routines for getting an initial estimate
c           of the galaxy shape parameters
      subroutine GetGalaxyShape(Incl,PA,Center,CatItem)
      implicit none
      real,INTENT(INOUT)  :: Incl,PA
      real,INTENT(IN)  :: Center(0:1)
      Type(CatalogueItem),INTENT(IN) :: CatItem
      real Ellip


      if(PFlags%ShapeSource .eq. 0) then
        call EstimateShape(ObservedMaps,Center
     &                      ,Incl,PA)
      elseif(PFlags%ShapeSource .eq. 1) then

        PA=CatItem%kinPA*Pi/180.+Pi/2.
        if(PA .gt. 2.*Pi) PA=PA-2.*Pi
        Ellip=CatItem%EllipseMin/CatItem%EllipseMaj
        Incl=acos(ellip)
      endif
    
c      PA=Pi/2.
c      Incl=75.5*Pi/180.
c      PA=56.5*Pi/180.
c      PA=PA+Pi/2.
      print*, "Shape huh", Incl*180./Pi,PA*180./Pi
      print*, "SoFiA shape", CatItem%EllipseInc,CatItem%EllipsePA
      Incl=CatItem%EllipseInc*Pi/180.
      PA=(CatItem%EllipsePA-90.)*Pi/180.
      if(PA .lt. 0.) then
        PA=PA+2.*Pi
      endif

      return
      end subroutine
cccccccc




ccccc
c       This routine gets the noise
c
c           THIS IS NOT COMPLETE YET AND ASSUMES THE USE OF A CATALOGUE VALUE
      subroutine GetNoise(Noise,CatItem)
      use sort
      implicit none
      Type(CatalogueItem),INTENT(INOUT) :: CatItem
      real,INTENT(INOUT) :: Noise
      real BeamArea
      real RMSTemp

      real SPeak,SInt
      integer nCells,medIndx

      real,ALLOCATABLE :: FlatFlux(:)
      integer,ALLOCATABLE :: FlatIndx(:)

c       Start by estimating the noise from the cube (*Note that this will not work with a noiseless cube)
c      Noise=0.02
c      RMSTemp=0.02
      call EstimateNoise(ObservedDC,RMSTemp)    !/src/PreAnalysis/EstimateCubeNoise.f
c       Set the Noise variable
      CatItem%RMS=RMSTemp
      Noise=RMSTemp

c       Get a number of S/N estimates
c           Start by getting the peak and total flux
      SPeak=maxval(MaskedObservedDC%Flux)
      SInt=sum(MaskedObservedDC%Flux)
      nCells=int(sum(DataCubeMask%Flux))
c           We also want the median flux, which requires some sorting of the masked array

      ALLOCATE(FlatIndx(nCells))
      ALLOCATE(FlatFlux(nCells))
c       Loop through the cells and fill the flattened flux array
      call MakeMaskedFlatFluxArr(MaskedObservedDC
     &          ,FlatFlux,nCells)
c       Sort the flattened array
      call indexx(nCells,FlatFlux,FlatIndx)
      MedIndx=nCells/2
c           Now we can get the median, peak, average, and integrated S/N (see Westmeier et al. 2021 and Deg et al. 2022) measures
      ObservedDC%DH%SN_Median=FlatFlux(FlatIndx(MedIndx))
     &          /Noise
      ObservedDC%DH%SN_Peak=SPeak/Noise
      ObservedDC%DH%SN_Avg=SInt/(nCells*Noise)
      ObservedDC%DH%SN_Int=SInt
     &                  /(sqrt(nCells*ObservedBeam%BeamAreaPixels)
     &                  *Noise)

c       Deallocate the flattened arrays
      DEALLOCATE(FlatFlux)
      DEALLOCATE(FlatIndx)

      return
      end subroutine
ccccccc


cccccc
c
      subroutine EstimateFlatVDisp(VDisp)
      implicit none
      real,INTENT(INOUT) :: VDisp

      VDisp=8.

      return
      end subroutine
ccccccc

ccccc
c       This routine takes flat-disk profiles and sets up a
c           tilted ring model
      subroutine ConvertFlatDiskProfilesToTR(nRings,Center,Incl
     &              ,PA,VSys,VDisp,EstimatedRadialProfiles)
      implicit none
      integer, INTENT(IN) :: nRings
      real,INTENT(IN) :: Incl, PA, VSys,Center(0:1)
      real,INTENT(IN) :: EstimatedRadialProfiles(0:2,0:nRings-1)
      real,INTENT(IN) :: VDisp
      integer i

c           Fully initialize the TR model-->src/PreAnalysis/ModellingInitializations.f
      call SetupTiltedRingStructures(nRings)
c           Set all parameters using the estimations done previously
      do i=0, nRings-1
        ModelTiltedRing%R(i)%Rmid=EstimatedRadialProfiles(0,i)
        ModelTiltedRing%R(i)%CentPos(0)=Center(0)
        ModelTiltedRing%R(i)%CentPos(1)=Center(1)
        ModelTiltedRing%R(i)%Inclination=Incl
        ModelTiltedRing%R(i)%PositionAngle=PA
        ModelTiltedRing%R(i)%VSys=VSys
        ModelTiltedRing%R(i)%VRot=EstimatedRadialProfiles(2,i)
        ModelTiltedRing%R(i)%VRad=0.
        ModelTiltedRing%R(i)%VDisp=VDisp
        ModelTiltedRing%R(i)%Vvert=0.
        ModelTiltedRing%R(i)%dvdz=0.
        ModelTiltedRing%R(i)%Sigma=EstimatedRadialProfiles(1,i)
        ModelTiltedRing%R(i)%logSigma
     &                  =log10(EstimatedRadialProfiles(1,i))
        ModelTiltedRing%R(i)%z0=0.
        ModelTiltedRing%R(i)%zGradiantStart=5.*ModelTiltedRing%R(i)%z0

        print*, 'initial SD and VRot',EstimatedRadialProfiles(1,i)
     &              ,EstimatedRadialProfiles(2,i)
c        print*, "hmmm", ModelTiltedRing%R(i)%RMid
cc     &          ,ModelTiltedRing%R(i)%Rwidth
      enddo

      return
      end subroutine
cccccccc


ccccc
      subroutine MakeMaskedFlatFluxArr(MaskCube,FlatFlux,TotCells)
      implicit none
      Type(DataCube),INTENT(IN) :: MaskCube
      integer, INTENT(IN) :: TotCells
      real,INTENT(INOUT) :: FlatFlux(TotCells)
      integer i,j,k,count
      real eps

c       Initialize the counter
      eps=1.e-20
      count=0
      do i=0,MaskCube%DH%nPixels(0)-1
        do j=0,MaskCube%DH%nPixels(1)-1
            do k=0,MaskCube%DH%nChannels-1
                if(abs(MaskCube%Flux(i,j,k)) .gt. eps) then
                    count=count+1
                    FlatFlux(count)=MaskCube%Flux(i,j,k)
                endif
            enddo
        enddo
      enddo
      return
      end subroutine
cccccccccc

      end module
