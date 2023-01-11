cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for
c       calculating maps from a data cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module MakeMomentMapsMod
      use DataCubeMod
      use CommonConsts


      contains
cccccccc
c       This routine calculates the first three moment maps for a cube
      subroutine MakeMomentMaps(DC,Maps)
      implicit none
      Type(DataCube),INTENT(IN) :: DC
      Type(DataCube),INTENT(INOUT) :: Maps

      integer Moment
      integer i,j
c      print*, "Making moment maps"
c      print*, shape(DC%Flux),idum

      i=50
      j=50
      Moment=0


    
     
      print*, "Make Moment Maps Check", sum(DC%Flux)
      print*, "Cube Channel Width", abs(DC%DH%ChannelSize)
      do i=0,DC%DH%nPixels(0)-1
        do j=0,DC%DH%nPixels(1)-1

            do Moment=0,2
                call CalculateGeneralMoment(DC%DH%nChannels-1
     &          ,DC%Channels(0:DC%DH%nChannels-1)
     &          ,DC%Flux(i,j,0:DC%DH%nChannels-1)
     &          ,Moment,Maps%Flux(i,j,Moment))
                if(Moment .eq. 0) then
                    if(Maps%Flux(i,j,Moment).lt.0.) then
                        Maps%Flux(i,j,Moment)=0.
                    endif
                endif
            enddo
c               Adjust the Moment 0 map flux by the channel width
            Maps%Flux(i,j,0)=Maps%Flux(i,j,0)*abs(DC%DH%ChannelSize)
c           Switch the second moment to a velocity dispersion
            Maps%Flux(i,j,2)=Maps%Flux(i,j,2)-
     &              Maps%Flux(i,j,1)**2.
            Maps%Flux(i,j,2)=sqrt(Maps%Flux(i,j,2))
c            print*, i,j,Maps%Flux(i,j,0:2)
        enddo
      enddo

      return
      end subroutine
cccccc

ccccccc
c
      subroutine CalculateGeneralMoment(n,Arr,WArr,WMom
     &                  ,CalcMoment)
      implicit none
      integer, INTENT(IN) :: n, WMom
      real, INTENT(IN) :: Arr(0:n-1), WArr(0:n-1)
      real, INTENT(OUT) :: CalcMoment
      integer i
      real term,wsum

      CalcMoment=0.
      wsum=0.
c       Loop through all the array elements
      if(WMom .eq. 0) then
        do i=0,n-1
            term=WArr(i)
            CalcMoment=CalcMoment+term
        enddo
      else
        do i=0,n-1
            wsum=wsum+WArr(i)
            term=(Arr(i)**real(WMom))*WArr(i)
            CalcMoment=CalcMoment+term
        enddo
        CalcMoment=CalcMoment/wsum
c        if(wsum .eq. 0) CalcMoment=0.
      endif
c      print*, n,wsum,CalcMoment,WMom



      return
      end subroutine
cccccccc


ccccccc
c       This subroutine makes a 1D velocity profile
      subroutine MakeVelProfile(DC,VelProf)
      implicit none
      Type(DataCube),INTENT(IN) :: DC
      real,ALLOCATABLE, INTENT(INOUT) :: VelProf(:,:)
      integer i, j, k


      Allocate(VelProf(0:1,0:DC%DH%nChannels-1))
      VelProf=0.
      do i=0, DC%DH%nChannels-1
        VelProf(0,i)=DC%Channels(i)
        do j=0,DC%DH%nPixels(0)-1
            do k=0,DC%DH%nPixels(1)-1
            VelProf(1,i)=VelProf(1,i)+DC%Flux(j,k,i)
            enddo
        enddo
c        print*, "VelProfile",i, VelProf(0:1,i)
      enddo

      return
      end subroutine
ccccccc


ccccccc
c
      subroutine RadioMomentMaps(DC,Maps)
      implicit none
      Type(DataCube),INTENT(IN) :: DC
      Type(DataCube),INTENT(INOUT) :: Maps

      integer Moment
      integer i,j,k

      real FluxSum
      real FluxLim

      FluxLim=1.e-8

c           Check if the maps array has already been allocated,
c           allocate the moment maps
      if (ALLOCATED(Maps%Flux) .eqv. .False.) then
        Maps%DH=DC%DH
        Maps%DH%nChannels=3
        call AllocateDataCube(Maps)
      endif

      Maps%Flux=0.          !Set the maps to zero
c       Loope through all pixels
      do i=0,DC%DH%nPixels(0)-1
        do j=0,DC%DH%nPixels(1)-1
c           Do a first loop over all channels
            do k=0,DC%DH%nChannels-1
c               Sum up the fluxes for moment 0
                Maps%Flux(i,j,0)=Maps%Flux(i,j,0)
     &                  +DC%Flux(i,j,k)
c               Sum up the fluxes*velocities for moment 1
                Maps%Flux(i,j,1)=Maps%Flux(i,j,1)
     &                  +DC%Flux(i,j,k)*DC%Channels(k)
            enddo
c               If the total flux is below some limit, set it to zero
            if(abs(Maps%Flux(i,j,0)) .le. FluxLim) then
                Maps%Flux(i,j,0)=0.
            endif
c               Normalize the moment 1 map by the total flux
            Maps%Flux(i,j,1)=Maps%Flux(i,j,1)/Maps%Flux(i,j,0)
c           Now loop through the channels again to get the moment 2 value
            do k=0, DC%DH%nChannels-1
                Maps%Flux(i,j,2)=Maps%Flux(i,j,2)
     &                  +DC%Flux(i,j,k)
     &                  *(DC%Channels(k)-Maps%Flux(i,j,1))**2.
            enddo
c           And normalize by the flux again
c            print*, "Radio check",i,j,Maps%Flux(i,j,0:2)
            Maps%Flux(i,j,2)=sqrt(Maps%Flux(i,j,2)/Maps%Flux(i,j,0))
        enddo
      enddo

      return
      end subroutine
ccccccc



ccccc
c
      subroutine CalcPVDiagrams(DC,Beam,Angle,Center,PVMaps)
      use BeamMod
      implicit none
      Type(DataCube),INTENT(IN) :: DC
      Type(Beam2D),INTENT(IN) :: Beam
      real,INTENT(IN) :: Angle,Center(2)

      Type(DataCube),INTENT(INOUT) :: PVMaps

      integer i,j,k
      real Size
      real AngUse
      real X,Y,XP,YP
      integer ii,jj,kk
      real BeamSize_Pixels
      

      print*, "Calculate a PV diagram"
      print*, "Center point", Center
      print*, "Angle Check", Angle*180./Pi
c       Check whether the PV maps have been allocated
      if (ALLOCATED(PVMaps%Flux) .eqv. .False.) then
        PVMaps%DH=DC%DH
        PVMaps%DH%nChannels=2
c       The Y-axis of the PV diagram should be the total number of channels
        PVMaps%DH%nPixels(1)=DC%DH%nChannels
c       The X-axis fo the PV diagram should be twice the maximum distance of the
c           central pixel to the largest dimension
        Size=0
        do i=0,1
            Size=max(Size,Center(i+1),DC%DH%nPixels(i)-Center(i+1))
        enddo
        PVMaps%DH%nPixels(0)=2*int(Size)+1
c           The reference values of the PV maps also need to be adjusted
        PVMaps%DH%RefLocation(0)=int(Size)
        PVMaps%DH%RefVal(0)=0.
        PVMaps%DH%RefLocation(1)=int(DC%DH%nChannels/2)
        PVMaps%DH%RefVal(1)=DC%Channels(int(PVMaps%DH%RefLocation(1))-1)
        PVMaps%DH%PixelSize(1)=DC%DH%ChannelSize
c        do i=0,1
c            print*,PVMaps%DH%RefLocation(i),PVMaps%DH%RefVal(i)
c     &              ,PVMaps%DH%nPixels(i)
c        enddo

        call AllocateDataCube(PVMaps)
      endif

c           Once the cube is allocated, the 2 PV diagrams can be calculated
c               The first thing is to get the beamsize in pixels
      BeamSize_Pixels=Beam%BeamMajorAxis
      do k=0,1
c           Set the using angle
        if (k .eq. 0) then
            AngUse=Angle
        else
            AngUse=Angle+Pi/2.
        endif
c           Ensure that AngUse is between 0 and 2 pi
100     continue
        if (AngUse .gt. 2.*Pi) then
            AngUse=AngUse-2.*Pi
            goto 100
        elseif(AngUse .lt. 0.) then
            AngUse=AngUse+2.*Pi
            goto 100
        endif
c       Now do a loop through the cube dimensions
        do i=0, DC%DH%nPixels(0)-1
            do j=0,DC%DH%nPixels(1)-1
c               Get the coordinates of the cell relative to the center in pixels
                X=real(i)-Center(1)
                Y=real(j)-Center(2)
c               Rotate this coordinates to the major(k=0) or minor(k=1) axis
c                   The negative angle is becaus this is a rotation back
                XP=X*cos(-AngUse)-Y*sin(-AngUse)
                YP=X*sin(-AngUse)+Y*cos(-AngUse)
c               Adjust these rotated coordinates back to cube indices
                ii=int(XP+Center(1))
                jj=int(YP+Center(2))
c               Now check that YP is within half a beam of the target axis
                if (abs(YP) .lt. BeamSize_Pixels/2.) then
c                   And make sure that ii is inside the PV dimension
                    if(ii .ge. 0 .and.
     &                  ii .le. PVMaps%DH%nPixels(0)-1) then
c                   Loop through the channels and add everything together
                        do kk=0,DC%DH%nChannels-1
c                       But we do need to check if the cube flux is somehow NaN'd
                            if(DC%Flux(i,j,kk)
     &                        .ne. DC%Flux(i,j,kk)+1) then
                                PVMaps%Flux(ii,kk,k)=
     &                                  PVMaps%Flux(ii,kk,k)
     &                                 +DC%Flux(i,j,kk)
                            endif
                        enddo

                    endif
                endif
            enddo
        enddo
      enddo


      return
      end subroutine
ccccccc

      end module
