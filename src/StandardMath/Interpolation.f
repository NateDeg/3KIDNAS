
      module InterpolateMod
      contains

cccc
      subroutine SimpleInterpolateX(P1,P2,YTarg,XIntVal)
c       This routine does a linear interpolation between 2 points to get the
c       XValue
      implicit none
      real,INTENT(IN) :: P1(2),P2(2), YTarg
      real,INTENT(INOUT) :: XIntVal
      real slope,delY

c      print*, P1, P2,YTarg
      slope=(P2(2)-P1(2))/(P2(1)-P1(1))
      delY=YTarg-P1(2)
      XIntVal=delY/slope+P1(1)
      return
      end subroutine
ccccc

ccc
      subroutine SimpleInterpolateY(P1,P2,XTarg,YIntVal)
c       This routine does a linear interpolation between 2 points to get the
c       YValue
      implicit none
      real,INTENT(IN) :: P1(2),P2(2), XTarg
      real,INTENT(INOUT) :: YIntVal
      real slope,delX

c      print*, "Linear interpolate",P1, P2,XTarg
      slope=(P2(2)-P1(2))/(P2(1)-P1(1))
      delX=XTarg-P1(1)
      YIntVal=delX*slope+P1(2)
      return
      end subroutine
ccccc


ccccc
c
c       Do a bilinear interpolation where PA, PB, PC, and PD are
c           coordinates of a square grid:
c
c               C3           C4
c                   Targ
c               C1           C2
c
cccccc
      subroutine BiLinearInterpolation(PTarg,Corners)
      implicit none
      real,INTENT(INOUT) :: PTarg(3)        !target X,Y,Val
      real,INTENT(IN) :: Corners(4,3) !corner X,Y,Val
    
      real CInterpolate(2,3)     !The interpolated points along the 1-2 and 3-4 lines
      real P1Temp(2),P2Temp(2),CTemp(2)

      integer i, j,k,kk
c       Set the X-target value
      CInterpolate(1:2,1)=PTarg(1)

      CTemp(1)=PTarg(1)
c       Loop across the 1-2 and 3-4 lines
      do i=1,2
c           Keep track of corner indices
        j=(i-1)*2+1
c           Set the two corner points to (X,Val)
        do k=1,2
            kk=(k-1)*2+1
            P1Temp(k)=Corners(j,kk)
            P2Temp(k)=Corners(j+1,kk)
        enddo
c        Interpolate along the line to the target x value
        call SimpleInterpolateY(P1Temp,P2Temp,CTemp(1),CTemp(2))
c        Store the 2D point cooridnates to the CInterploate array
        CInterpolate(i,2)=Corners(j,2)
        CInterpolate(i,3)=CTemp(2)
c        print*, "Interpolated points", CInterpolate(i,1:3)
      enddo
c       Now intepolate across the 2 intepolated points
      P1Temp(1:2)=CInterpolate(1,2:3)
      P2Temp(1:2)=CInterpolate(2,2:3)
      CTemp(1)=PTarg(2)
      call SimpleInterpolateY(P1Temp,P2Temp,CTemp(1),CTemp(2))
      PTarg(3)=CTemp(2)


      return
      end subroutine
cccccc


ccccc
c
c       Do a trilinear interpolation where corners are
c           coordinates of a  ube grid:
c
c                  C7                C8
c               C5              C6
c                       Targ
c                  C3               C4
c               C1               C2
c
cccccc

      subroutine TriLinearInterpolation(PTarg,Corners)
      implicit none
      real,INTENT(INOUT) :: PTarg(4)        !target X,Y,Z,Val
      real,INTENT(IN) :: Corners(8,4) !corner X,Y,Val

      real CornersTemp(4,3)
      real SurfacePoints(2,3)
      real PTemp(2,2),CTemp(2)

      integer i,j,k, ii,jj
c       Set the Surface points X-Y coordinates to that of the target
c           X-Y values
      SurfacePoints(1:2,1)=PTarg(1)
      SurfacePoints(1:2,2)=PTarg(2)

c       Loop through the two surfaces and do the bilinear interpolation
c           to get the center points of the surfaces
c           Set the corners to the lower Z and then the upper Z pairs
      jj=0
      do i=1,2
c        print*, "Doing surface", i
        ii=0    !   Set the surface corner counter to zero
        do j=1,2
            do k=1,2
                ii=ii+1 !   Increase the surface corner value
                jj=jj+1 !   Increase the cubic corner value
                CornersTemp(ii,1:2)=Corners(jj,1:2) !   Set the corner X-Y values
                CornersTemp(ii,3)=Corners(jj,4) !Set the values of the temp corners
c                print*, CornersTemp(ii,1:3)
            enddo
        enddo
c           Once the corners are set, run the bilinear interpolation on the surface
        call BiLinearInterpolation(SurfacePoints(i,1:3), CornersTemp)
c       Check the interpolated value
c        print*, "Surface points",SurfacePoints(i,1:3)
      enddo
c       Now do a linear interpolation between the 2 surface points
      CTemp(1)=PTarg(3)
   
      do i=1,2
        ii=(i-1)*4+1
        PTemp(i,1)=Corners(ii,3)
        PTemp(i,2)=SurfacePoints(i,3)
c        print*, "Temp points", PTemp(i,1:2)
      enddo
      call SimpleInterpolateY(PTemp(1,1:2),PTemp(2,1:2)
     &          ,CTemp(1),CTemp(2))
      PTarg(4)=CTemp(2)

    

      return
      end subroutine
cccccccc

      end module

