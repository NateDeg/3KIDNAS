c       This file contains routines that are useful for calculating inverse
c           trig functions over a full circle

      module FullCircTrig
      use CommonConsts
      contains
c

ccccccc
c       The full circle arc-tangent function
      subroutine FullCircATan(X,Y,Theta)
      implicit none
      real, INTENT(IN) :: X,Y
      real, INTENT(OUT) :: Theta

      Theta=atan(Y/X)

      if(X .ge. 0. .and. Y .ge. 0.) then
        Theta=Theta
      elseif(X .lt. 0. .and. Y .ge. 0.) then
        Theta=Theta+Pi
      elseif(X .lt. 0. .and. Y .lt. 0.) then
        Theta=Theta+Pi
      elseif(X .ge. 0. .and. Y .lt. 0.) then
        Theta=Theta+2.*Pi
      endif


      return
      end subroutine
cccccccc

      end module

