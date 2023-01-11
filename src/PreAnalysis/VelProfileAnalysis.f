cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       analyze a velocity profile
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module VelProfileAnalysisMod
      use CommonConsts


      implicit none
      contains

cccccc
      subroutine EstimateVSysAndEdges(nChannels,VelProfile
     &              ,FracLimit,VSys,Edges)
      implicit none
      integer, intent(IN) :: nChannels
      real,intent(IN) :: VelProfile(0:1,0:nChannels-1)
      real,intent(IN) :: FracLimit
      real,intent(OUT) :: VSys,Edges(0:1)

      integer nPeaks,PeakIndx(0:1)
      real PeakFlux(0:1),PeakVel(0:1)
      real SmoothFlux(0:nChannels-1),SmoothFluxDeriv(0:nChannels-2)
      integer nBoxCar
      real AvgLim

      integer i,j,k
      real VGuess
      integer EdgeGuessIndx(0:1),VGuessIndx,EdgeIndx(0:1)

c
      nPeaks=0
      SmoothFlux=0.

      nBoxCar=5
      j=nBoxCar/2
      AvgLim=1.
c       Smooth the profile
      call BoxCarAvg(nChannels,VelProfile(1,:),SmoothFlux,nBoxCar)
c           Get the derivative of the smoothed profile
      call FirstDeriv(nChannels,VelProfile(0,:)
     &          ,SmoothFlux,SmoothFluxDeriv)
c       Get the first guess at the Edges using the average flux
      call FluxLims(nChannels,SmoothFlux,AvgLim,EdgeGuessIndx)
      VGuessIndx=int(EdgeGuessIndx(1)+EdgeGuessIndx(0))/2
      VGuess=VelProfile(1,VGuessIndx)


c       Go through the derivs and look for sign changes
      do i=EdgeGuessIndx(0)-1,EdgeGuessIndx(1)
        if(SmoothFluxDeriv(i)*SmoothFluxDeriv(i+1) .lt. 0) then
            nPeaks=nPeaks+1
        endif
      enddo
c           Find the velocity peaks
      call FindPeakFlux(nChannels,nPeaks,VelProfile,EdgeGuessIndx
     &          ,VGuessIndx,PeakIndx,PeakVel,PeakFlux)
c           Use the peaks to get the edges
      call FindEdges(nChannels,nPeaks,VelProfile
     &                  ,PeakIndx,PeakFlux
     &                  ,FracLimit,EdgeIndx,Edges)
c           Use the edges to get VSys
      VSys=(Edges(0)+Edges(1))/2.
c      print*, "Estimated VSys", VSys

      return
      end subroutine
ccccccccc

cccccc
c           This routine implements a basic boxcar average to smooth an array
      subroutine BoxCarAvg(n,Arr,SmoothArr,nBoxCar)
      implicit none
      integer, INTENT(IN) :: n,nBoxCar
      real,intent(IN) :: Arr(0:n-1)
      real,intent(INOUT) :: SmoothArr(0:n-1)

      integer i,j,k,nUse

c       First set the smooth array to zero
      SmoothArr=0.
c       Now get the j index, which is half the number of box car units (should be odd)
      j=nBoxCar/2
c       Loop through all channels
      do i=0,n-1
c           Set nUse=0
        nUse=0
        do k=-j,j
c               Check if we are in the limits
            if((i+k .ge. 0) .and. (i+k) .le. n-1) then
                nUse=nUse+1
                SmoothArr(i)=SmoothArr(i)+Arr(i+k)
            endif
        enddo
        SmoothArr(i)=SmoothArr(i)/real(nUse)
c        print*, i,j,Arr(i),nUse,SmoothArr(i)
      enddo

      return
      end subroutine
ccccccccc

cccccc
c      This routine gets boundaries based on a flux percentage relative to the average flux
c           of the smooth profile
      subroutine FluxLims(n,SmoothFlux,FracLim,LimGuess)
      implicit none
      integer, INTENT(IN) :: n
      real,intent(IN) :: SmoothFlux(0:n-1),FracLim
      integer,intent(OUT):: LimGuess(0:1)

      real AvgFlux,FluxLim
      integer i

      AvgFlux=sum(SmoothFlux)/real(n)
      FluxLim=AvgFlux*FracLim
      LimGuess=-1

      do i=0,n-1
        if(SmoothFlux(i) .gt. FluxLim) then
            LimGuess(0)=i
            goto 100
        endif
      enddo

100   continue
      do i=n-1,0,-1
        if(SmoothFlux(i) .gt. FluxLim) then
            LimGuess(1)=i
            goto 200
        endif
      enddo
200   continue
c      print*, LimGuess

      return
      end subroutine
cccc

cccccc
c       This subroutine calculates the peak flux in some bounds
      subroutine FindPeakFlux(n,nPeaks,VelProfile,EdgeGuessIndx
     &          ,VGuessIndx,PeakIndx,PeakVel,PeakFlux)
      implicit none
      integer, INTENT(IN) :: n,nPeaks,VGuessIndx
      integer, INTENT(IN) :: EdgeGuessIndx(0:1)
      real,intent(IN) :: VelProfile(0:1,0:n-1)
      integer, INTENT(OUT) ::PeakIndx(0:1)
      real,INTENT(OUT) :: PeakFlux(0:1),PeakVel(0:1)
      integer i,j

      PeakFlux=-1.
      do i=VGuessIndx,EdgeGuessIndx(1)+1
        j=VGuessIndx-(i-VGuessIndx)
        if(nPeaks .eq. 1) then
            if(VelProfile(1,i) .gt. PeakFlux(0)) then
                PeakFlux(0)=VelProfile(1,i)
                PeakIndx(0)=i
                PeakVel(0)=VelProfile(0,i)
            endif
            if(VelProfile(1,j) .gt. PeakFlux(0)) then
                PeakFlux(0)=VelProfile(1,j)
                PeakIndx(0)=j
                PeakVel(0)=VelProfile(0,j)
            endif
        else
            if(VelProfile(1,i) .gt. PeakFlux(1)) then
                PeakFlux(1)=VelProfile(1,i)
                PeakIndx(1)=i
                PeakVel(1)=VelProfile(0,i)
            endif
            if(VelProfile(1,j) .gt. PeakFlux(0)) then
                PeakFlux(0)=VelProfile(1,j)
                PeakIndx(0)=j
                PeakVel(0)=VelProfile(0,j)
            endif
        endif
      enddo

      return
      end subroutine
cccccc


ccccc
c       Get the edges from the peak flux
      subroutine FindEdges(n,nPeaks,VelProfile,PeakIndx,PeakFlux
     &                  ,FracLimit,EdgeIndx,Edges)
      implicit none
      integer, INTENT(IN) :: n,nPeaks
      real,intent(IN) :: VelProfile(0:1,0:n-1),FracLimit
      integer, INTENT(IN) ::PeakIndx(0:1)
      real,INTENT(IN) :: PeakFlux(0:1)
      integer, INTENT(OUT) ::EdgeIndx(0:1)
      real, INTENT(OUT) :: Edges(0:1)

      integer i,j,k,l,Dir,Start,Fin
      real FluxTarg

      if(nPeaks .eq. 1) then
        k=0
      else
        k=1
      endif

      do j=0,1
        if(j .eq. 0) then
            Dir=-1
            Start=PeakIndx(0)
            Fin=0
            FluxTarg=PeakFlux(0)*FracLimit
        else
            Dir=1
            Start=PeakIndx(k)
            Fin=n-1
            FluxTarg=PeakFlux(k)*FracLimit
        endif
        do i=Start,Fin,Dir
            if(VelProfile(1,i) .le. FluxTarg) then
                EdgeIndx(j)=i
                if(j .eq. 0) then
                l=i
                else
                    l=i-1
                endif
                call LinearInterp_YTarg(VelProfile(0,l:l+1)
     &                      ,VelProfile(1,l:l+1)
     &                   ,FLuxTarg,Edges(j))
                goto 100
            endif
        enddo
100     continue
      enddo


      return
      end subroutine

cccccc
c       This routine calculates the linear first derivative.
c           It is meant to be used for interpolation only so it just gets
c           the segments between each pair of points
c
      subroutine FirstDeriv(n,X,Y,YP)
      implicit none
      integer, INTENT(IN) :: n
      real,intent(IN) :: X(0:n-1),Y(0:n-1)
      real,intent(INOUT) :: YP(0:n-2)

      integer i
      real dX,dY

      do i=0,n-2
        dX=X(i+1)-X(i)
        dY=Y(i+1)-Y(i)
        YP(i)=dY/dX
      enddo

      return
      end subroutine
cccccc


ccccc
c       This routine does a simple linear interpolation for a pair of
c       points
      subroutine LinearInterp_YTarg(X,Y,YTarg,XFound)
      implicit none
      real, INTENT(IN) :: X(0:1),Y(0:1)
      real, INTENT(IN) :: YTarg
      real, INTENT(OUT) :: XFound

      real dX,dY,delY,delX

      dX=X(1)-X(0)
      dY=Y(1)-Y(0)
      delY=YTarg-Y(0)
      delX=(delY/dY)*dX
      XFound=delX+X(0)


      return
      end subroutine
ccccccc

      end module
