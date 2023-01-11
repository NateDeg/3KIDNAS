cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       generate a bootstrap sample of a cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module PhysCoordMod
      use DataCubeMod
      use BeamMod


      implicit none

      contains


cccc
c       This routine builds an array of coordinates for a cube in
c           R_ellip, Theta, dV about a center and shape
      subroutine BuildPhysCoordsArray(XC,YC,VSys,PA,Inc
     &              ,CubeHeader,CoordArr)
      implicit none
      real,INTENT(IN):: XC,YC,VSys,PA,Inc
      Type(DataCubeHeader), INTENT(IN) :: CubeHeader
      real,ALLOCATABLE, INTENT(INOUT) :: CoordArr(:,:,:,:)

      real PhysCoords(3)
      integer i,j,k, PtIndx(3)

      print*, "Build a full physical coordinates array for a cube"
      ALLOCATE(CoordArr(3,0:CubeHeader%nPixels(0)-1
     &          ,0:CubeHeader%nPixels(1)-1
     &          ,0:CubeHeader%nChannels-1))

      do i=0, CubeHeader%nPixels(0)-1
        PtIndx(1)=i
        do j=0,CubeHeader%nPixels(1)-1
            PtIndx(2)=j
            do k=0,CubeHeader%nChannels-1
                PtIndx(3)=k
                call GetPhysCoords(XC,YC,VSys,PA,Inc
     &                     ,PtIndx,PhysCoords)
                CoordArr(1:3,i,j,k)=PhysCoords(1:3)
c                print*, i,j,k,CoordArr(1:3,i,j,k)
            enddo
        enddo
      enddo

      return
      end subroutine
cccccc


ccccc
c       This routine gets the elliptical coordinates and delta V in pixels,
c           angles, and channels for a pt index about some given center and geometry
      subroutine GetPhysCoords(XC,YC,VSys,PA,Inc
     &              ,PtIndx,PhysCoords)
      use CommonConsts
      implicit none

      real,INTENT(IN):: XC,YC,VSys,PA,Inc
      integer,INTENT(IN) :: PtIndx(3)
      real, INTENT(INOUT) :: PhysCoords(3)

      real X,Y,XRot, YRot

      real Theta, R, REllip, YEllip, Ellip
      
c      print*, "Get physical coords for each cell"

      X=real(PtIndx(1))-XC
      Y=real(PtIndx(2))-YC

c      print*, "X,Y from Cent", X,Y

      XRot=X*cos(-PA)-Y*sin(-PA)
      YRot=X*sin(-PA)+Y*cos(-PA)

c      print*, "Rotated points", XRot,YRot

      Ellip=cos(Inc)
c      print*, "Ellipticity", Ellip
      YEllip=YRot!/Ellip

      REllip=sqrt(XRot**2. + YEllip**2.)
c      print*, "Radius", sqrt(X**2.+Y**2.), REllip

      Theta=atan2(YRot,XRot)
      if(Theta .lt. 0.) Theta=Theta+2.*Pi
      if(Theta .gt. 2.*Pi) Theta=Theta-2.*Pi
c      print*, "Angle", Theta, Theta*180./3.14

      PhysCoords(1)=REllip
      PhysCoords(2)=Theta
c      PhysCoords(3)=(PtIndx(3)-VSys)/cos(Theta)
      PhysCoords(3)=(PtIndx(3)-VSys)

c           THINK ABOUT IF THIS IS GOOD TO USE--MAYBE A PROBLEM
c      if(abs(cos(Theta)) .lt. 0.01) then
c        PhysCoords(3)=(PtIndx(3)-VSys)/cos(Theta+0.01)
c      endif

      return
      end subroutine
cccccc


ccccc
c       This routine gets the elliptical coordinates and delta V in pixels,
c           angles, and channels for a pt index about some given center and geometry
      subroutine GetCubeCoords(XC,YC,VSys,PA,Inc
     &              ,CubePt,PhysCoords)
      implicit none

      real,INTENT(IN):: XC,YC,VSys,PA,Inc
      real, INTENT(IN) :: PhysCoords(3)
      real,INTENT(INOUT) :: CubePt(3)

      real X,Y,XRot, YRot
      real Theta, R, REllip, YEllip, Ellip


      REllip=PhysCoords(1)
      Theta=PhysCoords(2)

c      CubePt(3)=PhysCoords(3)*cos(Theta)+VSys
      CubePt(3)=PhysCoords(3)+VSys


      XRot=REllip*cos(Theta)
      YEllip=REllip*sin(Theta)

      Ellip=cos(Inc)

      YRot=YEllip!*Ellip

      X=XRot*cos(PA)-YRot*sin(PA)
      Y=XRot*sin(PA)+YRot*cos(PA)

      CubePt(1)=X+XC
      CubePt(2)=Y+YC


      return
      end subroutine
cccccc



      end module PhysCoordMod
