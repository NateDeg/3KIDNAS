cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       build a tilted ring model full of particles.
c          The routines assume the tilted ring has been allocated
c           and the individual ring parameters have been set.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingGenerationMod
      use TiltedRingMod
      use SingleRingGenerationMod
      use DataCubeMod
      use BeamMod

      implicit none
      contains

ccccc
      subroutine BuildTiltedRingModel(TR,idum,Noise,DC,BUse)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      integer,INTENT(INOUT) :: idum   !A seed for the random number generation
      Type(DataCube), INTENT(IN) :: DC
      Type(Beam2D), INTENT(IN) :: BUSe
      real, INTENT(IN) ::  Noise
      integer i
      real AvgChanPerPix

c      print*, "Building Tilted Ring"
c      print*, "Surf Dens Sanity", TR%R(1)%Sigma

      do i=0,TR%nRings-1
        call CalcAvgChanPerPix(i,TR,DC,BUSe,AvgChanPerPix)
        call Ring_CalcNumParticles(TR%R(i)
     &         , TR%cmode,TR%CloudBaseSurfDens,Noise
     &         , AvgChanPerPix)
        call Ring_ParticleAllocation(TR%R(i))               !Allocate the particle array
        call Ring_ParticleGeneration(TR%R(i),idum)          !Generate the particles in the ring
      enddo

      return
      end subroutine
cccccc


cccccc
c
      subroutine CalcAvgChanPerPix(ringIndx,TR,DC
     &              ,BUse,AvgChanPerPix)
      implicit none
      integer, INTENT(IN) :: ringIndx
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      Type(DataCube), INTENT(IN) :: DC
      Type(Beam2D), INTENT(IN) :: BUSe
      real, INTENT(INOUT) :: AvgChanPerPix

      real DV, DR,DDisp

    
c
      DV=2.*TR%R(ringIndx)%VRot/abs(DC%DH%ChannelSize)
      DDisp=2.*sqrt(2.)*TR%R(ringIndx)%VDisp/abs(DC%DH%ChannelSize)
      DR=2.*TR%R(ringIndx)%Rmid


c      print*, ringIndx,DR,DV,DDisp
c      print*, "Beam", BUse%BeamMajorAxis

      AvgChanPerPix=(BUse%BeamMajorAxis/DR)*(DV+DDisp)
      if (AvgChanPerPix .lt. DDisp+1) then
        AvgChanPerPix=DDisp+1.
      endif
c      print*, "Avg Chan", AvgChanPerPix
c      AvgChanPerPix=1.0
      return
      end subroutine
c
cccccccc



      end module
