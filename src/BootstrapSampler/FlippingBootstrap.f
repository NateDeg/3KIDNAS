cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       generate a bootstrap sample of a cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FlippingBootstrapMod
      use DataCubeMod
      use BeamMod
      use CubeDiffMod

      use DataCubeOutputsMod

      use PhysCoordMod
      use GenBootstrapMod

      implicit none

      contains

ccccccc
c       This is the main routine for generating a bootstrap sample.
c           It assumess that the NewCube hasn't been allocated.
      subroutine GenFlipBootstrapSample()
      use BootstrapGlobals
c      character(100) SampleFile
      implicit none

      real, ALLOCATABLE :: CoordArr(:,:,:,:)

      integer FlipType
 

c      print*, "Chan Block", VelBlockSize
c      print*, "Make a flipping bootstrap sample"
c       First make a difference cube using the observed cube and model cube
      call ConstructDiffCube(ObservedCube,ModelCube
     &          ,DifferenceCube)

c      SampleFile="Flipped_DifferenceCube.fits"
c      call WriteDataCubeToFITS(DifferenceCube
c     &          ,ObservedBeam,trim(SampleFile),"Test")



      call BuildPhysCoordsArray(BS_Cent%CentX, BS_Cent%CentY
     &                  ,BS_Cent%CentV, BS_Cent%PA, BS_Cent%Inc
     &                  ,ObservedCube%DH,CoordArr)
c       Copy the difference cube flux into the bootstrap cube
      BootstrapCube=DifferenceCube
c           Go through all channels and randomly flip across the major axis
      FlipType=1
      call AxisFlip(BootstrapCube,CoordArr
     &                  ,BS_Cent%CentX, BS_Cent%CentY
     &                  ,BS_Cent%CentV, BS_Cent%PA, BS_Cent%Inc
     &                  ,FlipType
     &                  ,VelBlockSize)

c      SampleFile="FlippedMajorAxis_DifferenceCube.fits"
c      call WriteDataCubeToFITS(BootstrapCube
c     &          ,ObservedBeam,trim(SampleFile),"Flipped")


      FlipType=2
      call AxisFlip(BootstrapCube,CoordArr
     &                  ,BS_Cent%CentX, BS_Cent%CentY
     &                  ,BS_Cent%CentV, BS_Cent%PA, BS_Cent%Inc
     &                  ,FlipType
     &                  ,VelBlockSize)

c      SampleFile="DoubleFlipped_DifferenceCube.fits"
c      call WriteDataCubeToFITS(BootstrapCube
c     &          ,ObservedBeam,trim(SampleFile))


      BootstrapCube%Flux=BootstrapCube%Flux+ModelCube%Flux

      DEALLOCATE(CoordArr)



      return
      end subroutine
cccccccc



ccccc
c       This routine goes through all the channels
c           and randomly chooses to flip across the major axis
      subroutine AxisFlip(Cube
     &              ,CoordArr
     &              ,XC,YC,VSys,PA,Inc
     &              ,FlipType
     &              ,VBlockSize)
      implicit none
      Type(DataCube), INTENT(INOUT) :: Cube
      real,INTENT(IN):: XC,YC,VSys,PA,Inc
      real,INTENT(IN) :: CoordArr(3,0:Cube%DH%nPixels(0)-1
     &                      ,0:Cube%DH%nPixels(1)-1
     &                      ,0:Cube%DH%nChannels-1)
      integer, INTENT(IN) :: FlipType
      real,INTENT(IN) :: VBlockSize
      integer VBlockInt, count

      integer i,j
      real RandVal
      Type(DataCube) TempCube

      TempCube=Cube
      VBlockInt=int(VBlockSize)
      if(VBlockInt .lt. 1) VBlockInt=1

      print*, "Flipping axis",FlipType

      count=0
      do i=0,Cube%DH%nChannels-1
c        print*, "hmm", count
c        do i=30,30
        call RANDOM_NUMBER(RandVal)
c        print*, i, RandVal
        if(RandVal .ge. 0.5) then
c            print*, "Flip Channel spatially"
            do j=1,VBlockInt
            
            call FlipChannelSpatial_Spectral_Dir(
     &                  Cube,TempCube
     &                  ,CoordArr,count
     &                  ,XC,YC,VSys,PA,Inc
     &                  ,FlipType)
            count=count+1
            if(count .ge. Cube%DH%nChannels) goto 100
            enddo
        endif
        count=count+1
        if(count .ge. Cube%DH%nChannels) goto 100
      enddo

100   continue
      Cube%Flux=TempCube%Flux
      call DeAllocateDataCube(TempCube)

      return
      end subroutine
ccccc


cccc
c
c       This routine flips a slice spatially
      subroutine FlipChannelSpatial_Spectral_Dir(
     &              Cube,TempCube
     &              ,CoordArr,ChanID
     &              ,XC,YC,VSys,PA,Inc,FlipType)
      use CommonConsts
      implicit none
      Type(DataCube), INTENT(IN) :: Cube
      Type(DataCube), INTENT(INOUT):: TempCube
      real,INTENT(IN):: XC,YC,VSys,PA,Inc
      real,INTENT(IN) :: CoordArr(3,0:Cube%DH%nPixels(0)-1
     &                      ,0:Cube%DH%nPixels(1)-1
     &                      ,0:Cube%DH%nChannels-1)
      integer,INTENT(IN):: ChanID
      integer,INTENT(IN):: FlipType

      real,ALLOCATABLE :: FlippedMap(:,:)
      integer i,j
      real ThetaNew,PhysCoordFlip(3), CubeCoordFlip(3)

      logical BoundCheck
      real Flux

c
c      print*, "Flip TYpe", FlipType
      ALLOCATE(FlippedMap(0:Cube%DH%nPixels(0)-1
     &              ,0:Cube%DH%nPixels(1)-1))

      do i=0,Cube%DH%nPixels(0)-1
        do j=0, Cube%DH%nPixels(1)-1
            if(FlipType .eq. 1) then
c            print*, i,j,CoordArr(1:3,i,j,ChanID)
                ThetaNew=2.*Pi-CoordArr(2,i,j,ChanID)
                PhysCoordFlip(1)=CoordArr(1,i,j,ChanID)
                PhysCoordFlip(3)=CoordArr(3,i,j,ChanID)
                PhysCoordFlip(2)=ThetaNew
            elseif(FlipType .eq. 2) then
                ThetaNew=CoordArr(2,i,j,ChanID)+Pi
                PhysCoordFlip(3)=-CoordArr(3,i,j,ChanID)
                if(ThetaNew .gt. 2.*Pi) ThetaNew=ThetaNew-2.*Pi
                PhysCoordFlip(1)=CoordArr(1,i,j,ChanID)
                PhysCoordFlip(2)=ThetaNew
            endif

            call GetCubeCoords(XC,YC,VSys,PA,Inc
     &              ,CubeCoordFlip,PhysCoordFlip)
c            print*, "Coord Comp",i,j,ChanID
c     &              ,CubeCoordFlip
c
            call SimpleBoundCheck_Real(Cube
     &                  ,CubeCoordFlip(1),CubeCoordFlip(2)
     &                  ,CubeCoordFlip(3),BoundCheck)
            if(BoundCheck) then
                call GetFluxAtPoint(Cube,CubeCoordFlip,Flux)
                FlippedMap(i,j)=Flux
            else
                FlippedMap(i,j)=Cube%Flux(i,j,ChanID)
            endif
        enddo
      enddo

      do i=0,Cube%DH%nPixels(0)-1
        do j=0, Cube%DH%nPixels(1)-1
            TempCube%Flux(i,j,ChanID)=FlippedMap(i,j)
        enddo
      enddo

      DEALLOCATE(FlippedMap)

      return
      end subroutine
ccccc




      end module FlippingBootstrapMod
