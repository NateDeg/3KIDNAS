cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines to calculate the moments 
c     of some a passed array of n-body particles per 
c     pixel element
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DataCubeMod
      implicit none

      Type DataCubeHeader
        integer nPixels(0:1),nChannels
        real PixelSize(0:1),ChannelSize
        real PixelCent(0:1),ChannelCent
        real Start(0:2)
        character(8) AxisType(0:2),Units(0:2)
        integer DimensionUnitSwitch(0:3)    !0->PixelSize, 1->ChannelSize, 2->RefVal(0:1), 3->RefVal(2)
        integer UncertaintyUnitSwitch
        real Uncertainty
        character(10) FUnit,FType
        real Epoch
        integer PixelCenterIndx(0:1),ChannelCenterIndx
        real RefLocation(0:2), RefVal(0:2)
        real StartFreq,RestFreq,delFreq
        integer nValid  !The number of non-null cells

        logical MaskSwitch      !If the cube is going to be used as a mask, some error checks can be surpressed
        real SN_Peak, SN_Avg, SN_Median, SN_Int  !  A set of S/N measurements
      end Type

      Type DataCube
        Type(DataCubeHeader) DH
        real,dimension(:,:),ALLOCATABLE :: Pixels
        real,dimension(:),ALLOCATABLE :: Channels
        real,dimension(:,:,:),ALLOCATABLE :: Flux
        
        integer,dimension(:),ALLOCATABLE :: FlattendValidIndices

      end Type

      contains
ccccccc
      subroutine AllocateDataCube(DC)
c
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      integer i, j,mPix
      real Delta
      integer nCells

      mPix=maxval(DC%DH%nPixels)
      ALLOCATE(DC%Pixels(0:1,0:mPix-1))
      ALLOCATE(DC%Channels(0:DC%DH%nChannels-1))
      ALLOCATE(DC%Flux(0:DC%DH%nPixels(0)-1,0:DC%DH%nPixels(1)-1
     &          ,0:DC%DH%nChannels-1))

c           Always initialize the number of valid indices to the full cube
c               in case not used elsewhere
      nCells=DC%DH%nPixels(0)*DC%DH%nPixels(1)*DC%DH%nChannels
      DC%DH%nValid=nCells
      ALLOCATE(DC%FlattendValidIndices(0:nCells-1))

      DC%Flux=0.

      do j=0,1
        Delta=(0.0-(DC%DH%RefLocation(j)))
        DC%DH%Start(j)=DC%DH%RefVal(j)+Delta*DC%DH%PixelSize(j)

c        print*, "Allocating DC", j,DC%DH%Start(j),DC%DH%RefLocation(j)
c     &              ,DC%DH%RefVal(j)
c     &              ,DC%DH%nPixels(j),DC%DH%PixelSize(j)
        do i=0,mPix-1
            Delta=(real((i))-(DC%DH%RefLocation(j)))
            DC%Pixels(j,i)=DC%DH%RefVal(j)+Delta*DC%DH%PixelSize(j)
c            print*, "Pix pos", j,i,DC%Pixels(j,i)
        enddo
      enddo

      j=2
      Delta=(0.0-(DC%DH%RefLocation(j)))
      DC%DH%Start(j)=DC%DH%RefVal(j)+(Delta
     &                  *DC%DH%ChannelSize)
      print*, "Vel Start",DC%DH%Start(j),Delta
     &          ,DC%DH%RefLocation(j),DC%DH%ChannelSize

      do i=0,DC%DH%nChannels-1
        Delta=(real((i))-(DC%DH%RefLocation(j)))
        DC%Channels(i)=DC%DH%RefVal(j)+Delta*DC%DH%ChannelSize
c        print*, "Channel Allocation Check", i,DC%Channels(i)
c     &              ,Delta,DC%DH%RefLocation(j)
c     &              ,DC%DH%ChannelSize
      enddo

      DC%DH%PixelCenterIndx=DC%DH%nPixels/2
      DC%DH%ChannelCenterIndx=DC%DH%nChannels/2


      return
      end subroutine
cccccccc

cccccc
      subroutine DeAllocateDataCube(DC)
      implicit none
      Type(DataCube) DC

      DEALLOCATE(DC%Pixels)
      DEALLOCATE(DC%Channels)
      DEALLOCATE(DC%Flux)
      DEALLOCATE(DC%FlattendValidIndices)

      return
      end subroutine
cccccccc


ccccccc
c       This routine calculates the index for a flattened array corresponding
c           to a 3D set of indices
c       It assumes that the loop goes through the x-y-v
      subroutine FlatIndxCalc(i,j,k,DH,l)
      implicit none
      Type(DataCubeHeader), INTENT(IN)::DH
      integer,INTENT(IN):: i,j,k
      integer,INTENT(INOUT) :: l

      l=0
      l=l+k
      l=l+j*DH%nChannels
      l=l+i*DH%nChannels*DH%nPixels(1)

      return
      end subroutine
cccccc

ccccc
c       This routine calculates the 3D indices from a flattened array index
c           It assumes that the loop order when set up is x-y-v
      subroutine ThreeDIndxCalc(l,DH,i,j,k)
      implicit none
      Type(DataCubeHeader),INTENT(IN)::DH
      integer, INTENT(IN) :: l
      integer, INTENT(INOUT) :: i,j,k
      integer ltemp

      ltemp=l
      i=int(ltemp/(DH%nChannels*DH%nPixels(1)))
      ltemp=ltemp-i*DH%nChannels*DH%nPixels(1)
      j=int(ltemp/DH%nChannels)
      ltemp=ltemp-j*DH%nChannels
      k=ltemp

      return
      end subroutine
ccccc

      end module
