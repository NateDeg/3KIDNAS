cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines MPI routines
c           that string names of the catalogue
c
c       This file is necessary due to some compilers not
c          recognizing that Scatterv is a generic function that can
c           take both reals and characters and integers
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module MPICharSplitMod
      use SoFiACatalogueMod
      use mpi

      contains

      subroutine ScatterStr(CatFull,CatLocal,Primary
     &              ,PID,nProcessors,ierr,nSplit,nPerProc
     &              ,SendA,TempStrPrimary,TempStr)
      implicit none

      Type(SoFiACat),INTENT(IN) :: CatFull
      Type(SoFiACat),INTENT(INOUT) :: CatLocal
      integer ,INTENT(IN) :: Primary,PID,nProcessors
      integer ierr,nSplit

      integer nPerProc(nProcessors), SendA(nProcessors)


      character TempStrPrimary(0:CatLocal%nObjects-1)
      character TempStr(0:CatLocal%nObjects-1)


      call MPI_Scatterv(TempStrPrimary(0:CatFull%nObjects-1)
     &          ,nPerProc,SendA,MPI_CHAR
     &          ,TempStr(0:CatLocal%nObjects-1)
     &          ,nPerProc(PID+1)
     &          ,MPI_CHAR
     &          ,Primary,MPI_COMM_WORLD,ierr)


      end subroutine
cccccc

      end module
