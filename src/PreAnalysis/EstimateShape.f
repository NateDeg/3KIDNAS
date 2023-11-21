cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       estimate the shape of the object
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module EstimateShapeMod
      use DataCubeMod
      use FullCircTrig
      use CommonConsts

      implicit none
      contains


ccccc
      subroutine EstimateShape(Maps,Center
     &              ,incl,Phi)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,intent(OUT) :: Phi,incl
      real,INTENT(IN) :: Center(0:1)

      real Moments(0:2),RMax
      real Q,U,D, ellip
      real A1,Phi1,incTest


      print*, "Estimating shape"

      call CalculateFlatMoments(Maps,Moments,RMax
     &              ,Center)
c       Calculate Q and U from the moments
      Q=Moments(0)-Moments(1)
      U=2*Moments(2)
c      Phi=0.5*atan(U/Q)
      call FullCircATan(Q,U,Phi)
      Phi=0.5*Phi
      D=sqrt(Q*Q+U*U)
      ellip=(1-D)/(1+D)
      incl=acos(ellip)

      if(Phi .le. 0.) Phi=Phi+2.*Pi
      print*, "V1",Q,U,Phi*180./Pi,ellip,incl*180./Pi
      print*, "Astro Phi", Phi*180./Pi-90.


c      call CalculateFourerA2Moments(Maps,A2,Phi2
c     &                  ,RMax,Center)

      call CalculateFourerA1VelMoments(Maps,A1,Phi1
     &                  ,RMax,Center)

      Phi=Phi1
c      incTest=acos(1.-A2)
c      print*, "Inclination Test comp", incTest*180./Pi

c      Q=(Moments(0)-Moments(1))/(Moments(0)+Moments(1))
c      U=2*Moments(2)/(Moments(0)+Moments(1))
c      Phi=0.5*atan(U/Q)
c      ellip=sqrt(Q*U/(cos(2.*Phi)*sin(2.*Phi)))
c      incl=acos(ellip)
c      print*, "V2",Q,U,Phi*180./Pi,ellip,incl*180./Pi


      return
      end subroutine
cccccc


ccccccc
c           This routine estimates the size of the object
      subroutine EstimateSize(Maps,Center,RMax)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(IN) :: Center(0:1)
      real,INTENT(OUT) :: RMax
      integer i,j
      real X,Y,R2
      RMax=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position relative to the center
            X=Maps%Pixels(0,i)-Center(0)
            Y=Maps%Pixels(1,j)-Center(1)
            R2=X*X+Y*Y      !Also get the radius
            if(Maps%Flux(i,j,0) .ne. 0) then
                if(sqrt(R2) .ge. RMax) RMax=sqrt(R2)
            endif
        enddo
      enddo

      return
      end subroutine
cccccccc


cccccc
c           This routine calculates the center using the center of brightness
      subroutine EstimateCenter(Maps,Center)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(OUT) :: Center(0:1)

      integer i, j
      real X,Y,fSum
c
c      print*, "Estimating Center"
      Center=0.
      fSum=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position relative to the center
c            X=Maps%Pixels(0,i)
c            Y=Maps%Pixels(1,j)
            X=i
            Y=j
            Center(0)=Center(0)+X*Maps%Flux(i,j,0)
            Center(1)=Center(1)+Y*Maps%Flux(i,j,0)
            fSum=fSum+Maps%Flux(i,j,0)
c            print*, i,j,X,Y,Maps%Flux(i,j,0),Center,fSum
        enddo
      enddo
      Center=Center/fSum
      print*, "Estimated Center", Center,fsum

      return
      end subroutine
cccccccc


ccccccc
      subroutine Iter_EstimateCenter(Maps,Center,CenterFlag)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(OUT) :: Center(0:1)
      integer,INTENT(INOUT) :: CenterFlag

      real CentTemp(0:1)
      real RLim, RMax,delR

      integer i,imax


      CenterFlag=0
c           First estimate the center using the full map
      call EstimateCenter(Maps,Center)
c       Now set the limiting radius a fraction of the map size
c           First get RMax
      RMax=sqrt(Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.)
      RMax=RMax/2.
c           Now set the radius limit
      RLim=0.5*Rmax
      call RLimited_EstimateCenter(Maps,Center,RLim,CentTemp)
c      print*, "1st RLimited center", CentTemp

c       Now iterate the center calculation until convergence
      imax=10
      do i=1,imax
c           First get the change in the center position
        delR=sqrt((Center(0)-CentTemp(0))**2.
     &          +(Center(1)-CentTemp(1))**2.)
c        print*, "Center CHeck",i, delR, Center, CentTemp,RLim
c       Now set the center to the R-limited value
        Center=CentTemp
c           If delR>0.5 pixels, then calculate the center again
        if(delR .gt. 0.5 .or. i .le.3) then
            call RLimited_EstimateCenter(Maps,Center
     &              ,RLim,CentTemp)
c            print*, "New Center", i,Center, CentTemp
            RLim=Rmax/(i+2)
c           Otherwise, finish with this center
        else
            return
        endif
      enddo
c       If we didn't converge, stick with the last center estimate
      print*, "Warning, intial center estimate did not converge"
      print*, "Using last center estimate",CentTemp
      Center=CentTemp
      CenterFlag=1


      return
      end subroutine
cccccc

cccccc
      subroutine RLimited_EstimateCenter(Maps,Center
     &              ,RLim,CentNew)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(IN) :: Center(0:1)
      real,INTENT(IN) :: RLim
      real,INTENT(OUT) :: CentNew(0:1)

      integer i, j
      real X,Y,R,XX,YY,fSum

      CentNew=0.
      fSum=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position in pixel counts
            X=real(i)
            Y=real(j)
c               Get the position relative to the current center
            XX=X-Center(0)
            YY=Y-Center(1)
c               Get the pixel radius from the current center
            R=sqrt(XX*XX+YY*YY)
c               Check that R is below the limiting radius
c            if (i .eq. 125) then
c                print*, i,j,X,Y,Center,XX,YY,R,RLim
c     &              ,Maps%Flux(i,j,0)
c            endif
            if(R .le. RLim) then
                CentNew(0)=CentNew(0)+X*Maps%Flux(i,j,0)
                CentNew(1)=CentNew(1)+Y*Maps%Flux(i,j,0)
                fSum=fSum+Maps%Flux(i,j,0)
            endif
c            print*, i,j,X,Y,Maps%Flux(i,j,0),Center,fSum
        enddo
      enddo
c      print*, "Rlimited Cent Test", CentNew,fSum,RLim
      CentNew=CentNew/fSum

      return
      end subroutine
cccccccc



cccccc
c
      subroutine CalculateFlatMoments(Maps,Moments,RMax,Center)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(OUT) :: Moments(0:2),RMax
      real,INTENT(IN) :: Center(0:1)
      integer i,j,nNonZero
      real X,Y,R2,RMax2,RLim
      real FTot

      print*, "Calculating moments"
      Moments=0.
      nNonZero=0
      RMax=0.

      RMax2=sqrt(Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.)
      RMax2=RMax2/2.
c           Now set the radius limit
      RLim=0.5*Rmax2

      print*, "Moment shape R max", RMax2,RLim
      FTot=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position relative to the center
c            X=Maps%Pixels(0,i)-Center(0)
c            Y=Maps%Pixels(1,j)-Center(1)

            X=real(i)-Center(0)
            Y=real(j)-Center(1)
            
            R2=X*X+Y*Y      !Also get the radius
            if(Maps%Flux(i,j,0) .ne. 0) then
c               Add a check in case the center is precisely at a pixel
c                   value
                if(R2 .ne. 0.
     &              .and. R2 .le. RLim**2.) then
c               Calculate the xx,yy,and xy 2 moments
                    Moments(0)=Moments(0)+X*X*Maps%Flux(i,j,0)
                    Moments(1)=Moments(1)+Y*Y*Maps%Flux(i,j,0)
                    Moments(2)=Moments(2)+X*Y*Maps%Flux(i,j,0)
                    nNonZero=nNonZero+1
                    FTot=FTot+Maps%Flux(i,j,0)
                    if(sqrt(R2) .ge. RMax) RMax=sqrt(R2)
                endif
c                print*, i,j,X,Y,R2,Moments,Center
            endif
        enddo
      enddo


      print*, "Moment Check", nNonZero,Moments,Center,FTot
      print*, "Multi check",Maps%DH%nPixels(0)*Maps%DH%nPixels(1)
c       Normalize by the number of non zero elements
c      Moments=Moments/real(nNonZero)
      FTot=Moments(0)+Moments(1)
      Moments=Moments/real(FTot)
      print*, "Normalized moments", Moments
      print*, "Moment check", Moments(0)+Moments(1)

      return
      end subroutine
cccccc




cccccc
c
      subroutine CalculateFourerA2Moments(Maps,A2,Phi
     &                  ,RMax,Center)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(OUT) :: A2,Phi,RMax
      real,INTENT(IN) :: Center(0:1)
      integer i,j,nNonZero
      real X,Y,R2,RMax2,RLim,theta
      real aa2,bb2
      real FTot

      print*, "Calculating Fourier moments"

      nNonZero=0
      RMax=0.

      RMax2=sqrt(Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.)
      RMax2=RMax2/2.
c           Now set the radius limit
      RLim=0.5*Rmax2

      FTot=0.
      aa2=0.
      bb2=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position relative to the center
            X=real(i)-Center(0)
            Y=real(j)-Center(1)

            R2=X*X+Y*Y      !Also get the radius
            call FullCircATan(X,Y,theta) !And theta
            if(Maps%Flux(i,j,0) .ne. 0) then
c               Add a check in case the center is precisely at a pixel
c                   value
                if(R2 .ne. 0.
     &              .and. R2 .le. RLim**2.) then
c               Calculate the xx,yy,and xy 2 moments
                    aa2=aa2+Maps%Flux(i,j,0)*cos(2.*theta)
                    bb2=bb2+Maps%Flux(i,j,0)*sin(2.*theta)
                    nNonZero=nNonZero+1
                    FTot=FTot+Maps%Flux(i,j,0)

                endif
c                print*, i,j,X,Y,R2,Moments,Center
            endif
        enddo
      enddo

      A2=sqrt((aa2**2.+bb2**2)/FTot**2.)
      call FullCircATan(aa2,bb2,Phi)
      Phi=0.5*Phi


      print*, "Fourier Mom values", aa2,bb2,FTot
     &          ,A2,Phi*180./Pi

      return
      end subroutine
cccccc


cccccc
c
      subroutine CalculateFourerA1VelMoments(Maps,A1,Phi
     &                  ,RMax,Center)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      real,INTENT(OUT) :: A1,Phi,RMax
      real,INTENT(IN) :: Center(0:1)
      integer i,j,nNonZero
      real X,Y,R2,RMax2,RLim,theta
      real aa1,bb1
      real FTot

      print*, "Calculating Fourier vel moments"

      nNonZero=0
      RMax=0.

      RMax2=sqrt(Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.)
      RMax2=RMax2/2.
c           Now set the radius limit
      RLim=0.5*Rmax2

      FTot=0.
      aa1=0.
      bb1=0.
c       Loop through all pixels
      do i=0,Maps%DH%nPixels(0)-1
        do j=0,Maps%DH%nPixels(1)-1
c               Get the current position relative to the center
            X=real(i)-Center(0)
            Y=real(j)-Center(1)

            R2=X*X+Y*Y      !Also get the radius
            call FullCircATan(X,Y,theta) !And theta
            if(Maps%Flux(i,j,0) .ne. 0) then
c               Add a check in case the center is precisely at a pixel
c                   value
                if(R2 .ne. 0.
     &              .and. R2 .le. RLim**2.) then
c               Calculate the a and b moment 1 components
                    aa1=aa1+Maps%Flux(i,j,1)*cos(1.*theta)
                    bb1=bb1+Maps%Flux(i,j,1)*sin(1.*theta)
                    nNonZero=nNonZero+1
                    FTot=FTot+Maps%Flux(i,j,0)

                    endif
c                print*, i,j,X,Y,R2,Moments,Center
            endif
        enddo
      enddo

      A1=sqrt((aa1**2.+bb1**2)/FTot**2.)
      call FullCircATan(aa1,bb1,Phi)
      Phi=Phi


      print*, "Fourier vel Mom1 values", aa1,bb1,FTot
     &          ,A1,Phi*180./Pi

      return
      end subroutine
cccccc



      end module
