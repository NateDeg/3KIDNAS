cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for
c       comparing the flux contained in 2 different cubes.
c        The cubes must have the same dimensions.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CubeCompareMod
      use LikelihoodMod
      use DataCubeMod



      contains
cccccccc
c
      subroutine CubeCompare(Cube1,Cube2,chi2,sigma)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube1,Cube2
      real,INTENT(OUT) :: chi2
      real,INTENT(IN) :: sigma
      integer nSpaxils              !The number of spaxils in the cubes
      real,dimension(:),ALLOCATABLE :: FlattenedF1,FlattenedF2
      real,dimension(:),ALLOCATABLE :: Uncertainties
      integer i,j,k,l,lUse

c      print*, "Comparing Cubes",sigma
      nSpaxils=Cube1%DH%nPixels(0)*Cube1%DH%nPixels(1)
     &              *Cube1%DH%nChannels
c      print*, "The number of elements is", nSpaxils
c     &          ,Cube2%DH%nPixels(0)*Cube2%DH%nPixels(1)
c     &              *Cube2%DH%nChannels
c      print*, "in Cube Compare", sum(Cube1%Flux), sum(Cube2%Flux)

c      ALLOCATE(FlattenedF1(0:nSpaxils-1))
c      ALLOCATE(FlattenedF2(0:nSpaxils-1))
c      ALLOCATE(Uncertainties(0:nSpaxils-1))


c       Flatten the flux arrays
c      do i=0, Cube1%DH%nPixels(0)-1
c        do j=0, Cube1%DH%nPixels(1)-1
c            do k=0, Cube1%DH%nChannels-1
c                l=i*Cube1%DH%nChannels*Cube1%DH%nPixels(1)
c     &              +j*Cube1%DH%nChannels+k
c                FlattenedF1(l)=Cube1%Flux(i,j,k)
c                FlattenedF2(l)=Cube2%Flux(i,j,k)
c                Uncertainties(l)=sigma
c            enddo
c        enddo
c      enddo

c      call LikePoint(chi2,nSpaxils,FlattenedF1
c     &                  ,FlattenedF2,Uncertainties)

c       Only use the non-nan cells for comparisons
      ALLOCATE(FlattenedF1(0:Cube1%DH%nValid-1))
      ALLOCATE(FlattenedF2(0:Cube1%DH%nValid-1))
      ALLOCATE(Uncertainties(0:Cube1%DH%nValid-1))


c       Loop through the valid cells
      do l=0, Cube1%DH%nValid-1
c           Get the specific flattened index
        lUse=Cube1%FlattendValidIndices(l)
c           Get the 3D indices
        call ThreeDIndxCalc(lUse,Cube1%DH,i,j,k) !/src/ObjectDefinitions/DataCube.f
c           Assign them to the flattened arrays
            FlattenedF1(l)=Cube1%Flux(i,j,k)
            FlattenedF2(l)=Cube2%Flux(i,j,k)
            Uncertainties(l)=sigma
c            print*, "Flux/index check", l,i,j,k
c     &          ,FlattenedF1(l),FlattenedF2(l)
c     &          ,Uncertainties(l)
      enddo

c      print*, "Compare Sanity Check"
c     &              , sum(FlattenedF1),sum(FlattenedF2)

c           Call the generic likelihood function
      call LikePoint(chi2,Cube1%DH%nValid,FlattenedF1
     &                  ,FlattenedF2,Uncertainties)
c       When finished, deallocate the flattendarrays
      DEALLOCATE(FlattenedF1)
      DEALLOCATE(FlattenedF2)
      DEALLOCATE(Uncertainties)

      return
      end subroutine

cccccc

      end module
