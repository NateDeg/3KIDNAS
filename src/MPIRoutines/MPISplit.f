cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines MPI routines
c           that split the catalogue
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module MPISplitMod
      use SoFiACatalogueMod
      use MPICharSplitMod
      use MPIIntSplitMod
      use mpi

      contains
cccccccc
c       This routine splits the catalogue into sub-catalogues in each processor
      subroutine SplitCatalogue(CatFull,CatLocal,Primary
     &                  ,PID,nProcessors)
      implicit none
      Type(SoFiACat),INTENT(IN) :: CatFull
      Type(SoFiACat),INTENT(INOUT) :: CatLocal
      integer ,INTENT(IN) :: Primary,PID,nProcessors
      integer ierr,nSplit

      integer nPerProc(nProcessors), SendA(nProcessors)
      integer i

      character,ALLOCATABLE :: TempStrPrimary(:),TempStr(:)

      print*, "Splitting catalogue across processors",PID
c           Get the number of items per processor.  This is calculated
c           in the primary processor as it is the only one with the full catalogue
      if(PID .eq. Primary) then
        nSplit=CatFull%nObjects/nProcessors
        print*, "Number to split", nSplit
        do i=1,nProcessors-1
            nPerProc(i)=nSplit
        enddo
        i=nProcessors
        nPerProc(i)=CatFull%nObjects-sum(nPerProc(1:nProcessors-1))
        print*, "Number per processor", nPerProc
      endif

c       Now broadcast that number across all processors
      call MPI_BCast(nPerProc,nProcessors,MPI_INTEGER,Primary
     &          ,MPI_COMM_WORLD,ierr)
      print*, "NUmber of objects in this processor", PID,nPerProc(PID+1)

c       Allocate size of the local catalogues
      CatLocal%nObjects=nPerProc(PID+1)
      call SoFiACatalogueAllocation(CatLocal)

c       Set up the sending arrays needed
      do i=1, nProcessors
        if(i .eq. 1) then
            SendA(i)=0
        else
            SendA(i)=sum(nPerProc(1:i-1))
        endif
      enddo


c       Scatter the IDs
      call ScatterInt(CatFull,CatLocal,Primary
     &              ,PID,nProcessors,ierr,nSplit,nPerProc
     &              ,SendA,TempStrPrimary,TempStr)

c       Scatter the Centers
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%RA
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%RA
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)


      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%DEC
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%DEC
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)

      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%ZCent
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%ZCent
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)

c       Scatter the Mins and Maxs
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%XMin
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%XMin
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%XMax
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%XMax
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%YMin
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%YMin
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%YMax
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%YMax
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%ZMin
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%ZMin
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%ZMax
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%ZMax
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
c       Scatter the ellipse parameters
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%EllipseMaj
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%EllipseMaj
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%EllipseMin
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%EllipseMin
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%EllipsePA
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%EllipsePA
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%kinPA
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%kinPA
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
c           Scatter the RMS
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%RMS
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%RMS
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
c           Scatter the w20, w50, and central frequency arrays
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%w20
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%w20
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(CatFull%Objects(0:CatFull%nObjects-1)%w50
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%w50
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(
     &          CatFull%Objects(0:CatFull%nObjects-1)%centralFreq
     &          ,nPerProc,SendA,MPI_REAL
     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%centralFreq
     &          ,nPerProc(PID+1)
     &          ,MPI_REAL
     &          ,Primary,MPI_COMM_WORLD,ierr)



c               Splitting the object names is tricky
c           The general idea is to make a temp array of characters in each position and
c               pass that array.  Then reassamble it to the full name
c       Allocate temporary strings for the names
      ALLOCATE(TempStrPrimary(0:CatFull%nObjects-1))
      ALLOCATE(TempStr(0:CatLocal%nObjects-1))

c       Loop through all the characters in a name
      do i=1,len(CatFull%Objects(0)%ObjName)
c           Have the primary array set up the string to be passed
        if(PID .eq. Primary) then
c            Put the i'th character from each name into an array
            TempStrPrimary(0:CatFull%nObjects-1)=
     &          CatFull%Objects(0:CatFull%nObjects-1)%ObjName(i:i)
        endif
c           Scatter that array
      call ScatterStr(CatFull,CatLocal,Primary
     &              ,PID,nProcessors,ierr,nSplit,nPerProc
     &              ,SendA,TempStrPrimary,TempStr)

c        call MPI_Scatterv(TempStrPrimary(0:CatFull%nObjects-1)
c     &          ,nPerProc,SendA,MPI_CHAR
c     &          ,TempStr(0:CatLocal%nObjects-1)
c     &          ,nPerProc(PID+1)
c     &          ,MPI_CHAR
c     &          ,Primary,MPI_COMM_WORLD,ierr)
c           Reassemble the character name
        CatLocal%Objects(0:CatLocal%nObjects-1)%ObjName(i:i)=
     &          TempStr(0:CatLocal%nObjects-1)
      enddo
c           Deallocate the temporary strings
      DEALLOCATE(TempStrPrimary)
      DEALLOCATE(TempStr)


c      print*, "Local IDs", PID
c     &          ,CatLocal%Objects(0:CatLocal%nObjects-1)%ObID
c     &          ,trim(CatLocal%Objects(0)%ObjName)



      return
      end subroutine
cccccc


      end module
