cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       generate a bootstrap sample of a cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GenBootstrapMod
      use DataCubeMod
      use BeamMod
      use CubeDiffMod

      use DataCubeOutputsMod

      use PhysCoordMod

      implicit none

      contains

ccccccc
c       This is the main routine for generating a bootstrap sample.
c           It assumess that the NewCube hasn't been allocated.
      subroutine GenBootstrapSample()
      use BootstrapGlobals
c      character(100) SampleFile
      implicit none
 
      integer TestPt(3)
      real CentX, CentY, CentV, PA, Inc
      real PhysCoords(3)
      real CubePt(3)


      real,ALLOCATABLE :: CoordArr(:,:,:,:)


      print*, "Full Bootstrap resample"
c       First make a difference cube using the observed cube and model cube
      call ConstructDiffCube(ObservedCube,ModelCube
     &          ,DifferenceCube)

      SampleFile="Symmetric_Asymmetric_DifferenceCube.fits"
      call WriteDataCubeToFITS(DifferenceCube
     &          ,ObservedBeam,trim(SampleFile)
     &          ,"TestDifference")





      call BuildPhysCoordsArray(BS_Cent%CentX, BS_Cent%CentY
     &                  ,BS_Cent%CentV, BS_Cent%PA, BS_Cent%Inc
     &                  ,ObservedCube%DH,CoordArr)
      

c       Now resample the difference cube and place it into the bootstrap cube
c      call BlockResampleCube(BootstrapCube,DifferenceCube
c     &                      ,MaskCube,ObservedBeam
c     &                      ,SpatialBlackSize,VelBlockSize)

      call BlockResampleCube_Phys2(BootstrapCube,DifferenceCube
     &                      ,CoordArr,ObservedBeam
     &                      ,SpatialBlackSize,VelBlockSize
     &                      ,BS_Cent%CentX, BS_Cent%CentY
     &                      ,BS_Cent%CentV, BS_Cent%PA, BS_Cent%Inc)

      SampleFile="Resampled_DifferenceCube_1x1x2.fits"
      call WriteDataCubeToFITS(BootstrapCube
     &          ,ObservedBeam,trim(SampleFile)
     &              ,"TestResample")

c       Now add the model back to the resampled cube
      BootstrapCube%Flux=BootstrapCube%Flux+ModelCube%Flux

c       When finished, deallocate the coordinate array
c      print*, "Coord array shape", shape(CoordArr)
      DEALLOCATE(CoordArr)

      return
      end subroutine
cccccccc


cccc
c       This routine resamples some cube
c
      subroutine BlockResampleCube(NewCube,BaseCube
     &              ,CubeMask,Beam
     &              ,SblockSize,VblockSize)
      implicit none
      Type(DataCube),INTENT(IN) :: BaseCube,CubeMask
      Type(Beam2D), INTENT(IN) :: Beam
      real,INTENT(IN) ::SblockSize,VblockSize
      Type(DataCube), INTENT(INOUT) :: NewCube

      integer bSizePix, bSizeChan

      integer nblocks_Dim(0:2)
      integer i,j,k, blockID(3)

      real,ALLOCATABLE :: DataBlock(:,:,:)

      print*, "Resampling cube"

c           First initialize the new cube that will be filled with block samples
c       Set the header to match the
      NewCube%DH=BaseCube%DH
      call AllocateDataCube(NewCube)        !/src/ObjectDefinitions/DataCube.f
c       Next set up the block size in pixels/channels
      bSizeChan=NINT(VblockSize)
      bSizePix=NINT(SblockSize*Beam%BeamMajorAxis)
      print*, "Block size", bSizePix,bSizeChan
c           Now we can allocate the data block array that will be filled with the base cube's data
      ALLOCATE(DataBlock(bSizePix,bSizePix,bSizeChan))

      do i=0,1
        nblocks_Dim(i)=BaseCube%DH%nPixels(i)/bSizePix
      enddo
      nblocks_Dim(2)=BaseCube%DH%nChannels/bSizeChan
      print*, "Number of blocks per dim", nblocks_Dim


c       Loop through all the blocks
      do i=0,nblocks_Dim(0)
c      do i=0,0
        blockID(1)=i
        do j=0,nblocks_Dim(1)
c        do j=0,0
            blockID(2)=j
c            do k=0,0
            do k=0,nblocks_Dim(2)
                blockID(3)=k
c                   For each block, select a subsection of the base cube
                call SelectDataBlock(BaseCube,DataBlock
     &              ,bSizePix,bSizeChan,CubeMask)
c                   Now that we have a data block, we can fill in the appropriate section of the data cube
                call FillInCubeByBlock(NewCube,DataBlock
     &              ,bSizePix,bSizeChan,blockID)

            enddo
        enddo
      enddo

c       Get rid of the data block memory
      DEALLOCATE(DataBlock)

      return
      end subroutine
ccccccc


ccccc
c           This routine selects a subsection of data from the
c               main cube that may be correlated
      subroutine SelectDataBlock(BaseCube,DataBlock
     &              ,bSizePix,bSizeChan,CubeMask)
      implicit none

      integer, INTENT(IN) :: bSizePix,bSizeChan
      Type(DataCube),INTENT(IN) :: BaseCube,CubeMask
      real, INTENT(INOUT) :: DataBlock(bSizePix,bSizePix,bSizeChan)

      integer i,j,k,ii,jj,kk
      integer LowCut(3),HighCut(3)
      integer RanSize(3)
      integer Cent(3)
      real RandVal
c
c           Set up the limits for the random selections to make sure everything stays inside the cube
      do i=1,2
        RanSize(i)=BaseCube%DH%nPixels(i-1)-bSizePix
      enddo
      RanSize(3)=BaseCube%DH%nChannels-bSizeChan
c       Nose select the low limits randomly and set the upper limits accordingly
100   continue
      do i=1,3
        call RANDOM_NUMBER(RandVal)
        LowCut(i)=int(RanSize(i)*RandVal)
        if(i .le. 2) then
            HighCut(i)=LowCut(i)+bSizePix-1
            Cent(i)=LowCut(i)+((bSizePix-1)/2)
        else
            HighCut(i)=LowCut(i)+bSizeChan-1
            Cent(i)=LowCut(i)+((bSizeChan-1)/2)
        endif
      enddo
c      print*, LowCut, Cent
c           Check that the midpoint of the block is inside the mask
      if(CubeMask%Flux(Cent(1),Cent(2),Cent(3)) .lt. 1) goto 100

      DataBlock=0.
c       Now fill in the data block
      ii=0
      do i=LowCut(1),HighCut(1)
        ii=ii+1
        jj=0
        do j=LowCut(2),HighCut(2)
            jj=jj+1
            kk=0
            do k=LowCut(3),HighCut(3)
                kk=kk+1
                DataBlock(ii,jj,kk)=BaseCube%Flux(i,j,k)
            enddo
        enddo
      enddo
                

      return
      end subroutine
ccccc


ccccc
c       This subroutine fills in a section of the new cube
c           using the filled data block
      subroutine FillInCubeByBlock(NewCube,DataBlock
     &              ,bSizePix,bSizeChan
     &              ,BlockID)
      implicit none

      integer, INTENT(IN) :: bSizePix,bSizeChan
      real, INTENT(IN) :: DataBlock(bSizePix,bSizePix,bSizeChan)
      integer,INTENT(IN):: BlockID(3)
      Type(DataCube),INTENT(INOUT) :: NewCube
c
      integer CubeStart(3)
      integer i,j,k
      integer ii,jj,kk
      logical BoundCheck

c      print*, "Filling in block section", BlockID
      do i=1,2
        CubeStart(i)=BlockID(i)*bSizePix
      enddo
      CubeStart(3)=BlockID(3)*bSizeChan

c      print*, "Cube Start", CubeStart
      do i=1,bSizePix
        ii=CubeStart(1)+i-1
        do j=1,bSizePix
            jj=CubeStart(2)+j-1
            do k=1,bSizeChan
                kk=CubeStart(3)+k-1
c                print*, i,j,k,ii,jj,kk
                call SimpleBoundCheck(NewCube,ii
     &                  ,jj,kk,BoundCheck)
                if(BoundCheck) then
                    NewCube%Flux(ii,jj,kk)=DataBlock(i,j,k)
                endif
            enddo
        enddo
      enddo


      return
      end subroutine
ccccccc


ccccc
c       This subroutine does a very simple check to make sure all indices are in bounds
c
      subroutine SimpleBoundCheck(Cube,i,j,k,BoundCheck)
      implicit none
      Type(DataCube), INTENT(IN) :: Cube
      integer,INTENT(IN) :: i,j,k
      logical, INTENT(INOUT) :: BoundCheck

      BoundCheck=.True.
      if( i .ge. Cube%DH%nPixels(0).or. i .lt. 0) then
        BoundCheck=.False.
      elseif( j .ge. Cube%DH%nPixels(1).or. i .lt. 0) then
        BoundCheck=.False.
      elseif( k .ge. Cube%DH%nChannels.or. i .lt. 0) then
        BoundCheck=.False.
      endif

      return
      end subroutine
cccccc

ccccc
c       This subroutine does a very simple check to make sure all indices are in bounds
c
      subroutine SimpleBoundCheck_Real(Cube,i,j,k,BoundCheck)
      implicit none
      Type(DataCube), INTENT(IN) :: Cube
      real,INTENT(IN) :: i,j,k
      logical, INTENT(INOUT) :: BoundCheck

      BoundCheck=.True.
      if( i .ge. real(Cube%DH%nPixels(0)).or. i .lt. 0) then
        BoundCheck=.False.
      elseif( j .ge. real(Cube%DH%nPixels(1)).or. j .lt. 0) then
        BoundCheck=.False.
      elseif( k .ge. real(Cube%DH%nChannels).or. k .lt. 0) then
        BoundCheck=.False.
      endif

      return
      end subroutine
cccccc




cccc
c       This routine resamples some cube
c
      subroutine BlockResampleCube_Phys(NewCube,BaseCube
     &              ,CoordArr,Beam
     &              ,SblockSize,VblockSize)
      implicit none
      Type(DataCube),INTENT(IN) :: BaseCube
      real,INTENT(IN) :: CoordArr(3,0:BaseCube%DH%nPixels(0)-1
     &                      ,0:BaseCube%DH%nPixels(1)-1
     &                      ,0:BaseCube%DH%nChannels-1)
      Type(Beam2D), INTENT(IN) :: Beam
      real,INTENT(IN) ::SblockSize,VblockSize
      Type(DataCube), INTENT(INOUT) :: NewCube

      real DeltaRange(3)

      integer bSizePix, bSizeChan

      integer nblocks_Dim(0:2)
      integer i,j,k, blockID(3),CentIndx(3)

      real,ALLOCATABLE :: DataBlock(:,:,:)
      

      print*, "Resampling cube"

c           First initialize the new cube that will be filled with block samples
c       Set the header to match the
      NewCube%DH=BaseCube%DH
      call AllocateDataCube(NewCube)        !/src/ObjectDefinitions/DataCube.f
c       Next set up the block size in pixels/channels
      bSizeChan=NINT(VblockSize)
      bSizePix=NINT(SblockSize*Beam%BeamMajorAxis)
      print*, "Block size", bSizePix,bSizeChan
c           Now we can allocate the data block array that will be filled with the base cube's data
      ALLOCATE(DataBlock(bSizePix,bSizePix,bSizeChan))

      do i=0,1
        nblocks_Dim(i)=BaseCube%DH%nPixels(i)/bSizePix
      enddo
      nblocks_Dim(2)=BaseCube%DH%nChannels/bSizeChan
      print*, "Number of blocks per dim", nblocks_Dim


      DeltaRange(1)=5.
      DeltaRange(2)=2.*3.14
      DeltaRange(3)=3.

c       Loop through all the blocks
      do i=0,nblocks_Dim(0)-1
c      do i=0,0
        blockID(1)=i
        CentIndx(1)=(i+0.5)*(bSizePix)
        do j=0,nblocks_Dim(1)-1
c        do j=0,0
            blockID(2)=j
            CentIndx(2)=(j+0.5)*(bSizePix)
c            do k=0,0
            do k=0,nblocks_Dim(2)-1
                blockID(3)=k
                CentIndx(3)=(k+0.5)*(bSizeChan)
c                   For each block, select a subsection of the base cube
                call SelectDataBlock_Phys(BaseCube
     &              ,DataBlock
     &              ,bSizePix,bSizeChan, CentIndx
     &              ,CoordArr
     &              ,DeltaRange)
c                   Now that we have a data block, we can fill in the appropriate section of the data cube
                call FillInCubeByBlock(NewCube,DataBlock
     &              ,bSizePix,bSizeChan,blockID)

            enddo
        enddo
      enddo

c       Get rid of the data block memory
      DEALLOCATE(DataBlock)

      return
      end subroutine
ccccccc



ccccc
c           This routine selects a subsection of data from the
c               main cube that may be correlated
      subroutine SelectDataBlock_Phys(BaseCube,DataBlock
     &              ,bSizePix,bSizeChan, PtIndx
     &              ,CoordArr
     &              ,DeltaRange)
      implicit none

      integer, INTENT(IN) :: bSizePix,bSizeChan
      Type(DataCube),INTENT(IN) :: BaseCube
      real,INTENT(IN) :: CoordArr(3,0:BaseCube%DH%nPixels(0)-1
     &                      ,0:BaseCube%DH%nPixels(1)-1
     &                      ,0:BaseCube%DH%nChannels-1)
      real, INTENT(INOUT) :: DataBlock(bSizePix,bSizePix,bSizeChan)

      integer,INTENT(IN) :: PtIndx(3)
      real,INTENT(IN):: DeltaRange(3)

      integer i,j,k,ii,jj,kk
      integer LowCut(3),HighCut(3)
      integer RanSize(3)
      integer Cent(3)
      real RandVal
      real Block_RandDiffs(3)
      integer Dir, Temp
c
c           Set up the limits for the random selections to make sure everything stays inside the cube
      do i=1,2
        RanSize(i)=BaseCube%DH%nPixels(i-1)-bSizePix
      enddo
      RanSize(3)=BaseCube%DH%nChannels-bSizeChan
c       Nose select the low limits randomly and set the upper limits accordingly
100   continue
      do i=1,3
        call RANDOM_NUMBER(RandVal)
        LowCut(i)=int(RanSize(i)*RandVal)
        if(i .le. 2) then
            HighCut(i)=LowCut(i)+bSizePix-1
            Cent(i)=LowCut(i)+((bSizePix-1)/2)
        else
            HighCut(i)=LowCut(i)+bSizeChan-1
            Cent(i)=LowCut(i)+((bSizeChan-1)/2)
        endif
      enddo
c      print*, "Blck Center Coords", CoordArr(1:3
c     &          ,PtIndx(1),PtIndx(2),PtIndx(3))
c      print*, "Rand Center Coords",CoordArr(1:3
c     &          ,Cent(1),Cent(2),Cent(3))
      do i=1,3
        if(i .le. 2) then
            if(Cent(i) .ge. BaseCube%DH%nPixels(i-1)) goto 100
        else
            if(Cent(i) .ge. BaseCube%DH%nChannels) goto 100
        endif
      enddo

      do i=1,3
        Block_RandDiffs(i)=CoordArr(i
     &          ,PtIndx(1),PtIndx(2),PtIndx(3))
     &          -CoordArr(i
     &          ,Cent(1),Cent(2),Cent(3))
c        print*, i, Block_RandDiffs(i)
        if(abs(Block_RandDiffs(i)) .gt. DeltaRange(i)) goto 100
      enddo
      print*, "Diff", Block_RandDiffs
c      print*, LowCut, Cent
c           Check that the midpoint of the block is inside the mask
c      if(CubeMask%Flux(Cent(1),Cent(2),Cent(3)) .lt. 1) goto 100

      Dir=1
c      do i=1,3
c        call RANDOM_NUMBER(RandVal)
c        if(RandVal .lt. 0.5) then
c            Dir=1
c        else
c            Dir=-1
c            Temp=LowCut(i)
c            LowCut(i)=HighCut(i)
c            HighCut(i)=Temp
c        endif
c      enddo

      DataBlock=0.
c       Now fill in the data block
      ii=0
      do i=LowCut(1),HighCut(1),Dir
        ii=ii+1
        jj=0
        do j=LowCut(2),HighCut(2),Dir
            jj=jj+1
            kk=0
            do k=LowCut(3),HighCut(3),Dir
                kk=kk+1
                DataBlock(ii,jj,kk)=BaseCube%Flux(i,j,k)
            enddo
        enddo
      enddo
                

      return
      end subroutine
ccccc



cccc
c       This subroutine initializes the quantities needed for block resampling
      subroutine BlockResample_Ini(NewCube,BaseCube
     &              ,Beam
     &              ,SblockSize,VblockSize
     &              ,bSizeChan, bSizePix, nblocks_Dim
     &              ,DataBlock)
      implicit none
      Type(DataCube),INTENT(IN) :: BaseCube
      Type(Beam2D), INTENT(IN) :: Beam
      real,INTENT(IN) ::SblockSize,VblockSize
      Type(DataCube), INTENT(INOUT) :: NewCube

      real,ALLOCATABLE, INTENT(INOUT) :: DataBlock(:,:,:)

      integer,INTENT(INOUT) :: bSizeChan, bSizePix, nblocks_Dim(0:2)
      integer i



c           First initialize the new cube that will be filled with block samples
c       Set the header to match the
      NewCube%DH=BaseCube%DH
      call AllocateDataCube(NewCube)        !/src/ObjectDefinitions/DataCube.f
c       Next set up the block size in pixels/channels
      bSizeChan=NINT(VblockSize)
      bSizePix=NINT(SblockSize*Beam%BeamMajorAxis)
      print*, "Block size", bSizePix,bSizeChan
c           Now we can allocate the data block array that will be filled with the base cube's data
      ALLOCATE(DataBlock(bSizePix,bSizePix,bSizeChan))

      do i=0,1
        nblocks_Dim(i)=BaseCube%DH%nPixels(i)/bSizePix
      enddo
      nblocks_Dim(2)=BaseCube%DH%nChannels/bSizeChan
      print*, "Number of blocks per dim", nblocks_Dim


      return
      end subroutine
ccccccccc



cccc
c       This routine resamples some cube
c
      subroutine BlockResampleCube_Phys2(NewCube,BaseCube
     &              ,CoordArr,Beam
     &              ,SblockSize,VblockSize
     &              ,XC,YC,VSys,PA,Inc)
      implicit none
      Type(DataCube),INTENT(IN) :: BaseCube
      real,INTENT(IN) :: CoordArr(3,0:BaseCube%DH%nPixels(0)-1
     &                      ,0:BaseCube%DH%nPixels(1)-1
     &                      ,0:BaseCube%DH%nChannels-1)
      Type(Beam2D), INTENT(IN) :: Beam
      real,INTENT(IN) ::SblockSize,VblockSize
      Type(DataCube), INTENT(INOUT) :: NewCube

      real,INTENT(IN):: XC,YC,VSys,PA,Inc

      real DeltaRange(3)

      integer bSizePix, bSizeChan

      integer nblocks_Dim(0:2)
      integer i,j,k, blockID(3),CentIndx(3)

      real,ALLOCATABLE :: DataBlock(:,:,:)
      

      print*, "Resampling cube"

c           First initialize the new cube that will be filled with block samples
c       Set the header to match the
      call BlockResample_Ini(NewCube,BaseCube
     &              ,Beam
     &              ,SblockSize,VblockSize
     &              ,bSizeChan, bSizePix, nblocks_Dim
     &              ,DataBlock)


      DeltaRange(1)=0.
      DeltaRange(2)=2.*3.14
      DeltaRange(3)=0.0

c       Loop through all the blocks
      do i=0,nblocks_Dim(0)-1
c      do i=0,0
        blockID(1)=i
        CentIndx(1)=(i+0.5)*(bSizePix)
        do j=0,nblocks_Dim(1)-1
c        do j=8,8
            blockID(2)=j
            CentIndx(2)=(j+0.5)*(bSizePix)
c            do k=0,0
            do k=0,nblocks_Dim(2)-1
c                print*, i,j,k
                blockID(3)=k
                CentIndx(3)=(k+0.5)*(bSizeChan)

                call Build_DataBlock_PhysSelect(BaseCube
     &              ,DataBlock
     &              ,bSizePix,bSizeChan, CentIndx
     &              ,CoordArr
     &              ,DeltaRange
     &              ,XC,YC,VSys,PA,Inc)

c                   For each block, select a subsection of the base cube
c                call SelectDataBlock_Phys(BaseCube
c     &              ,DataBlock
c     &              ,bSizePix,bSizeChan, CentIndx
c     &              ,CoordArr
c     &              ,DeltaRange)
c                   Now that we have a data block, we can fill in the appropriate section of the data cube
                call FillInCubeByBlock(NewCube,DataBlock
     &              ,bSizePix,bSizeChan,blockID)

            enddo
        enddo
      enddo

c       Get rid of the data block memory
      DEALLOCATE(DataBlock)

      return
      end subroutine
ccccccc




ccccc
c           This routine selects a subsection of data from the
c               main cube that may be correlated
      subroutine Build_DataBlock_PhysSelect(BaseCube,DataBlock
     &              ,bSizePix,bSizeChan, PtIndx
     &              ,CoordArr
     &              ,DeltaRange
     &              ,XC,YC,VSys,PA,Inc)
      implicit none

      integer, INTENT(IN) :: bSizePix,bSizeChan
      Type(DataCube),INTENT(IN) :: BaseCube
      real,INTENT(IN) :: CoordArr(3,0:BaseCube%DH%nPixels(0)-1
     &                      ,0:BaseCube%DH%nPixels(1)-1
     &                      ,0:BaseCube%DH%nChannels-1)
      real, INTENT(INOUT) :: DataBlock(bSizePix,bSizePix,bSizeChan)

      integer,INTENT(IN) :: PtIndx(3)
      real,INTENT(IN):: DeltaRange(3)

      real,INTENT(IN):: XC,YC,VSys,PA,Inc

      integer i,j,k,ii,jj,kk

      real Delta(3), RandVal
      real AdjustedPt(3), AdjustedCubePt(3), Flux
      logical BoundCheck
c
c           First randomly select the delta values
100   continue
      do i=1,3
        call RANDOM_NUMBER(RandVal)
        Delta(i)=(2.*RandVal-1) *DeltaRange(i)
      enddo
c      print*, "Shift in physical coordinates", Delta,PtIndx
c      Delta(2)=1.0*2.*3.14

c           Loop through the datablock
      do i=1,bSizePix
c       Get the index for the coordinate array
        ii=PtIndx(1)-bSizePix/2+i-1
        do j=1, bSizePix
            jj=PtIndx(2)-bSizePix/2+j-1
            do k=1, bSizeChan
                kk=PtIndx(3)-bSizeChan/2+k-1
            
c                print*, i,j,k,ii,jj,kk,PtIndx
c                print*, "Phys Pt Coords", CoordArr(1:3,ii,jj,kk)
                AdjustedPt(1:3)=CoordArr(1:3,ii,jj,kk)+Delta(1:3)
c                   Make sure the new value of R is positive
                if(AdjustedPt(1) .lt. 0.) goto 100

c                print*, "Adjusted Pt Coords", AdjustedPt
                call GetCubeCoords(XC,YC,VSys,PA,Inc
     &              ,AdjustedCubePt,AdjustedPt)
c                print*, "Adjusted Cube Coords", AdjustedCubePt
c                   Check if the points are inbounds
                call SimpleBoundCheck_Real(BaseCube
     &                  ,AdjustedCubePt(1),AdjustedCubePt(2)
     &                  ,AdjustedCubePt(3),BoundCheck)
                if(BoundCheck .eqv. .False.) goto 100

                call GetFluxAtPoint(BaseCube,AdjustedCubePt,Flux)
                DataBlock(i,j,k)=Flux
c                print*, i,j,k,Flux

            enddo
        enddo
      enddo
c      print*, "Final Delta Selection", Delta


      return
      end subroutine
ccccc



ccccc
c           This routine gets the flux at a specific point via interpolation
      subroutine GetFluxAtPoint(Cube,Pt,Flux)
      use InterpolateMod
      implicit none

      Type(DataCube),INTENT(IN) :: Cube
      real, INTENT(IN) :: Pt(3)
      real, INTENT(INOUT) :: Flux

      integer i,j,k,ii,jj,kk,ll,CurrIndx(3)
      real CornerPts(8,4),InterpolatePt(4)
      logical BoundCheck


      do i=1,3
        CurrIndx(i)=int(Pt(i))
      enddo

c       Set the interpolated point position
      InterpolatePt(1:3)=Pt(1:3)

      do i=1,2
        ii=CurrIndx(3)+(i-1)
        do j=1,2
            jj=CurrIndx(2)+(j-1)
            do k=1,2
                kk=CurrIndx(1)+(k-1)
c                   Get the count for the corners
                ll=k+(j-1)*2+(i-1)*2*2
c               Set the corner point values
                CornerPts(ll,1)=real(kk)
                CornerPts(ll,2)=real(jj)
                CornerPts(ll,3)=real(ii)
c                print*, kk,jj,ii,ll
                call SimpleBoundCheck(Cube,kk,jj,ii,BoundCheck)
c                print*, "Corner Pts", ll,CornerPts(ll,1:3)
c     &                  ,BoundCheck,Size
                if(BoundCheck) then
                    CornerPts(ll,4)=Cube%Flux(kk,jj,ii)
c                    print*,"Corner checck", DC%Flux(kk,jj,ii)
                else
                    CornerPts(ll,4)=0.
                endif
c                print*, "Corner Pts", ll,CornerPts(ll,1:4)
            enddo
        enddo
      enddo
c           Use trilinear interpolation to get the flux at the specified point
      call TriLinearInterpolation(InterpolatePt,CornerPts) !/src/StandardMath/Interpolation.f
c      print*, "Interpolated flux", InterpolatePt
      Flux=InterpolatePt(4)

      return
      end subroutine
ccccc
      end module GenBootstrapMod
