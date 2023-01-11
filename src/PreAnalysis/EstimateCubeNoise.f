cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       estimate the shape of the object
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module EstimateCubeNoiseMod
      use DataCubeMod
      use CommonConsts

      implicit none
      contains


ccccc
      subroutine EstimateNoise(Cube,RMS)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube
      real,intent(INOUT) :: RMS

      integer count
      integer ii,jj,kk,CountLims(3,2)
      real frac
      integer DimSwitch
      real sum

      print*, "Estimating cube noise"

      frac=0.12

      count=0
      sum=0.
      do ii=1,2
        DimSwitch=1
        call SetLims(Cube,frac,DimSwitch,ii
     &                      ,CountLims(DimSwitch,1:2))
        do jj=1,2
            DimSwitch=2
            call SetLims(Cube,frac,DimSwitch,ii
     &                      ,CountLims(DimSwitch,1:2))
            do kk=1,2
                DimSwitch=3
                call SetLims(Cube,frac,DimSwitch,ii
     &                      ,CountLims(DimSwitch,1:2))
c       Once all the limits are calculated, get square of the flux in each subbox
                call GetSubBox(Cube,CountLims,count,sum)


            enddo
        enddo
      enddo
      print*, "Noise square count", count,sum,sqrt(sum/count)
      RMS=sqrt(sum/count)


      return
      end subroutine
cccccc

ccccc
      subroutine SetLims(Cube,frac,DimSwitch,SideSwitch,Lims)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube
      real,intent(IN) ::frac
      integer,intent(IN) :: DimSwitch,SideSwitch
      integer,INTENT(INOUT):: Lims(2)

      real UpperLim
c       First get the size of the dimension
      if (DimSwitch .eq. 1) then
        UpperLim=Cube%DH%nPixels(0)
      elseif(DimSwitch .eq. 2) then
        UpperLim=Cube%DH%nPixels(1)
      elseif(DimSwitch .eq. 3) then
        UpperLim=Cube%DH%nChannels
      endif

      if(SideSwitch .eq. 1) then
        Lims(1)=0
        Lims(2)=int(frac*UpperLim)
      else
        Lims(1)=(1-frac)*(UpperLim-1)
        Lims(2)=UpperLim-1
      endif
c      print*, Lims,UpperLim

      return
      end subroutine
cccccc

ccccccc
      subroutine GetSubBox(Cube,CountLims,count,sum)
      implicit none
      Type(DataCube), INTENT(IN) :: Cube
      integer, INTENT(IN) :: CountLims(3,2)
      integer, INTENT(INOUT) :: count
      real, INTENT(INOUT) :: sum
      integer i,j,k


      do i=CountLims(1,1),CountLims(1,2)
        do j=CountLims(2,1),CountLims(2,2)
            do k=CountLims(3,1),CountLims(3,2)
                sum=sum+Cube%Flux(i,j,k)**2.
                count=count+1
            enddo
        enddo
      enddo

      return
      end subroutine
cccccc

      end module
