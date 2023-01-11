cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines MPI routines
c           that integers of the catalogue
c
c       This file is necessary due to some compilers not
c          recognizing that Scatterv is a generic function that can
c           take both reals and characters and integers
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module MPIIntSplitMod
      use SoFiACatalogueMod
      use mpi

      contains

      subroutine ScatterInt(CatFull,CatLocal,Primary
     &              ,PID,nProcessors,ierr,nSplit,nPerProc
     &              ,SendA,TempStrPrimary,TempStr)
      implicit none

      Type(SoFiACat),INTENT(IN) :: CatFull
      Type(SoFiACat),INTENT(INOUT) :: CatLocal
      integer ,INTENT(IN) :: Primary,PID,nProcessors
      integer ierr,nSplit

      integer nPerProc(nProcessors), SendA(nProcessors)


      character TempStrPrimary(*),TempStr(*)



c           Scatter the object ids
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%ObID
     &          ,nPerProc,SendA,MPI_INTEGER
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%ObID
     &          ,nPerProc(PID+1)
     &          ,MPI_INTEGER
     &          ,Primary,MPI_COMM_WORLD,ierr)


c       Scatter the number of channels
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%nChannels
     &          ,nPerProc,SendA,MPI_INTEGER
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%nChannels
     &          ,nPerProc(PID+1)
     &          ,MPI_INTEGER
     &          ,Primary,MPI_COMM_WORLD,ierr)
c       Scatter the number of LOS pixels
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%n_LOS
     &          ,nPerProc,SendA,MPI_INTEGER
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%n_LOS
     &          ,nPerProc(PID+1)
     &          ,MPI_INTEGER
     &          ,Primary,MPI_COMM_WORLD,ierr)


      end subroutine
cccccc

      end module
