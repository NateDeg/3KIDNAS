cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the pre-galaxy analysis flow.
c       The goal is to get good initial guesses for the various
c       model parameters.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module PostGalaxyAnalysisMod
      use PipelineGlobals
 
      use MakeMomentMapsMod


      implicit none
      contains

ccccc
      subroutine PostGalaxyAnalysis()
      implicit none

c      integer i, nPartTot
      print*, "Post fitting analysis"


c       After making the best fitting model, generate a set of model maps
      print*, "Making moment maps"
      call RadioMomentMaps(ModelDC,ModelMaps)
c       Also make a set of PV diagrams
      print*, "Making PV diagrams"
      call CalcPVDiagrams(ModelDC,ObservedBeam
     &                  ,ModelTiltedRing%R(0)%PositionAngle
     &                  ,ModelTiltedRing%R(0)%CentPos
     &                  ,ModelPVMaps)


      print*, "Getting distance and mass estimates"
      call DistanceEstimate()
c      call MassFromFlux(MaskedObservedDC)
      call MassFromFlux(ObservedDC)

c       For simplicity and possible future analysis, give the
c           chi^2 statistic to the galaxy dictionary object
      GalaxyDict%GoodnessOfFit=PVModel%BestLike
      GalaxyDict%nCells=ObservedDC%DH%nPixels(0)
     &                      *ObservedDC%DH%nPixels(1)
     &                      *ObservedDC%DH%nChannels
      GalaxyDict%Norm_GoodnessOfFit=GalaxyDict%GoodnessOfFit
     &              /real(GalaxyDict%nCells)

c       Do a quick number of particles test
c      nPartTot=0
c      do i=0,ModelTiltedRing%nRings-1
c        nPartTot=nPartTot+ModelTiltedRing%R(i)%nParticles
c      enddo
c      print*, "Total number of particles in final model", nPartTot

      return
      end subroutine
cccccc


ccccccc
      subroutine DistanceEstimate()
      use CommonConsts
      implicit none

      GalaxyDict%Distance=ModelTiltedRing%R(0)%VSys/H0

      return
      end subroutine
cccccc


cccccc
      subroutine MassFromFlux(DC)
      implicit none
      real FTot, M1
      Type(DataCube), INTENT(IN) :: DC

      FTot=sum(DC%Flux)
      M1=0.236*(GalaxyDict%Distance*1000.)**2.*FTot
      GalaxyDict%logMass=log10(M1
     &              *abs(DC%DH%Channelsize))

      return
      end subroutine
ccccccc


      end module
