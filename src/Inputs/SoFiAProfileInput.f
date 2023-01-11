cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading a
c         SoFiA catalogue file
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SofiaProfileInputMod

      use PipelineGlobals
      use UnitConvertMod
      use CommonConsts

      implicit none

      contains

ccccccc
c           This routine reads a SoFiA frequency/velocity profile.  It assumes the
c               profile is in frequency space and converts it to velocity space
      subroutine ReadSoFiAVelocityProfile()
      implicit none

      integer i, imax,HeaderLen,BodyLen
      character(1) LeadChar
      integer ChanNum
      real z

c       Figure out the size of the profile
      open(10,file=ProfileName,status='old')
      imax=1e5
      HeaderLen=0
      BodyLen=0
      do i=1,imax
        read(10,*,end=200) LeadChar
        if(LeadChar .eq. '#') then
            HeaderLen=HeaderLen+1
        else
            BodyLen=BodyLen+1
        endif
      enddo
200   continue
c       Now allocate the velocity profile
      ALLOCATE(ObservedVelocityProfile(0:1,0:BodyLen-1))
c       Rewind the file
      rewind(10)
c       Skip through the header
      do i=1, HeaderLen
        read(10,*) LeadChar
      enddo
c       Read in the profile
      do i=1, BodyLen
        read(10,*) ChanNum, ObservedVelocityProfile(0:1,i-1)
c           Since the profile is in frequency, we need to convert to
c           velocities in km/s
        call RedshiftCalc(z,HIRestFreq,ObservedVelocityProfile(0,i-1))
        ObservedVelocityProfile(0,i-1)=z*lightspeed
c        print*, i, HIRestFreq
c        print*,ObservedVelocityProfile(0:1,i-1)
      enddo


      close(10)

      return
      end subroutine
cccccccc


      end module
