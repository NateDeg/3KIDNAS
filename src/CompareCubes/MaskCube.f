cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines for masking
c       some data cube.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module MaskCubeMod
      use DataCubeMod



      contains
cccccccc
c           This routine masks a cube.  It assumes that the mask is set to 1 where
c           there is signal and zero elsewhere
      subroutine MaskCube(Cube,Mask)
      implicit none
      Type(DataCube),INTENT(INOUT) :: Cube
      Type(DataCube),INTENT(IN) :: Mask
      integer i,j,k

c       Loop through all elements
      do i=0, Cube%DH%nPixels(0)-1
        do j=0, Cube%DH%nPixels(1)-1
            do k=0, Cube%DH%nChannels-1
c                print*, i,j,k,Mask%Flux(i,j,k)
c     &                  ,Cube%Flux(i,j,k)
                Cube%Flux(i,j,k)=Cube%Flux(i,j,k)
     &                      *Mask%Flux(i,j,k)
            enddo
        enddo
      enddo

      return
      end subroutine
cccccc


cccccccc
c           This routine makes a basic cube according to some
c               flux tolerance
      subroutine MakeMaskCube(Cube,Mask,tol)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube
      Type(DataCube),INTENT(INOUT) :: Mask
      real,INTENT(IN) :: tol
      integer i,j,k
c       Flatten the flux arrays
      do i=0, Cube%DH%nPixels(0)-1
        do j=0, Cube%DH%nPixels(1)-1
            do k=0, Cube%DH%nChannels-1
                if(Cube%Flux(i,j,k) .lt. tol) then
                    Mask%Flux(i,j,k)=0.
                else
                    Mask%Flux(i,j,k)=1.0
                endif
            enddo
        enddo
      enddo

      return
      end subroutine
cccccc


      end module
