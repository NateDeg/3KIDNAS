cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading in a
c       datacube.  It also contains routines for getting
c       a text version of a data cube object header.
c       Since the beam parameters are usually contained in the
c       datacube headers, this file also gets those in combined
c       input files
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DataCubeInputMod

      use DataCubeMod
      use BeamMod
      use CommonConsts
      use CalcBeamKernelMod
      use InputUnitConversionsMod

      implicit none

      contains

ccccccc
c           This routine reads in a text version of a data cube object header
      subroutine ReadTextDataCubeHeader(DC,Beam,HeaderFileName)
      implicit none
      Type(DataCube) DC
      Type(Beam2D) Beam
      character(*) HeaderFileName
      logical FileCheck

      print*, "Getting the data cube input file"
      INQUIRE(file=HeaderFileName,EXIST=FileCheck)
      if(FileCheck .eqv. .False.) then
        print*, "Data cube header file does not exist "
     &          , trim(HeaderFileName)
        stop
      endif

      open(10,file=trim(HeaderFileName),status='old')
c       Get the number of pixels and channels
      read(10,*)
      read(10,*) DC%DH%nPixels(0:1), DC%DH%nChannels
c       Get the units for the datacube sizes, reference positions, beam axis, and beam angle
      read(10,*)
      read(10,*) DC%DH%DimensionUnitSwitch(0:3)
     &                  ,Beam%BeamUnitsSwitch(0:1)
c       Get the dimensions of the pixels and channels
      read(10,*)
      read(10,*) DC%DH%PixelSize(0:1), DC%DH%ChannelSize
c       Get the reference locations of each dimension
      read(10,*)
      read(10,*) DC%DH%RefLocation(0:2)
c       Get the reference values for each dimension
      read(10,*)
      read(10,*) DC%DH%RefVal(0:2)
c      DC%DH%RefVal(0:1)=DC%DH%RefVal(0:1)*3600
c       Get the beam size/shape
      read(10,*)
      read(10,*)  Beam%BeamMajorAxis, Beam%BeamMinorAxis
     &              ,Beam%BeamPositionAngle

c      Beam%BeamPositionAngle=Beam%BeamPositionAngle*Pi/180.


c       Get the number of sigma lengths to calculate the beam to
      read(10,*)
      read(10,*) Beam%SigmaLengths
c       Get the type of velocity smoothing being used
      read(10,*)
      read(10,*) Beam%VelocitySmoothSwitch
c       If using Gaussian smoothing read in the PSF
      if(Beam%VelocitySmoothSwitch .eq. 1) then
        read(10,*)
        read(10,*) Beam%VelocitySmoothSigma
      endif
      close(10)


c       Do the datacube and beam unit conversions
      call TextDCHeader_AngularConversions(DC,Beam)


c       Finally allocate the data cube and the beam
      call AllocateDataCube(DC)     !/src/ObjectDefinitions/DataCube.f
c               It is necessary to set the beam pixel size to datacube pixel size
      Beam%PixelSize=DC%DH%PixelSize
      call Allocate_Beam2D(Beam,DC%DH%nPixels)   !/src/ObjectDefinitions/Beam.f
c       Also calculate the real kernel for the beam
c      call Calculate2DBeamKernel(Beam,Beam%PixelSize)

      return
      end subroutine
cccccccc


cccccc
c
      subroutine ReadFullDataCube(DC,Beam,FileName,MaskUseSwitch)
      implicit none
      integer,INTENT(IN) :: MaskUseSwitch
      Type(DataCube), INTENT(INOUT) :: DC
      character(*),INTENT(IN) :: FileName
      Type(Beam2D), INTENT(INOUT) :: Beam

      logical FileCheck

c      print*, "Reading full data cube ",trim(FileName)
c      print*, trim(FileName)

      INQUIRE(file=FileName,EXIST=FileCheck)
      if(FileCheck .eqv. .False.) then
        print*, "Data cube file does not exist "
     &          , trim(FileName)
        stop
      endif

c           Set the datacube mask switch first
      if(MaskUseSwitch .eq. 0) then
        DC%DH%MaskSwitch = .False.
      else
        DC%DH%MaskSwitch =.True.
      endif
c       Read in the header
      call ReadDCHeader(FileName,DC,Beam)
c       Allocate the data cube
      call AllocateDataCube(DC)

      call readDCFile(FileName, DC)

      print*, "input DC flux", sum(DC%flux)


      return
      end subroutine

cccccccc



ccccccc
c       This routine is meant to read the header for a datacube object
c           It is based on the fitsio routines.  It also includes the beam
c           object header input
      subroutine ReadDCHeader(FileName,DC,Beam)
      implicit none
      character(*),INTENT(IN) :: FileName
      Type(DataCube), INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(INOUT) :: Beam


      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character record*80

      integer naxes
      integer nfound
      integer k
      integer StartVal

      integer*8 nInts(3)
      real vReal(3)


c      print*, "Reading DC Header", trim(FileName)
C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)


C     open the FITS file, with read-only access.  The returned BLOCKSIZE
C     parameter is obsolete and should be ignored.
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

      j = 0
100   continue
      j = j + 1


C  Get the data cube key words
      call ftghsp(unit,nkeys,nspace,status)
c      print*, "Get key words", unit, nkeys,nspace,status


      naxes=3
      StartVal=1
      call ftgknk(unit,'NAXIS',StartVal,naxes,nInts,nfound,status)
      DC%DH%nPixels=nInts(1:2)
      DC%DH%nChannels=nInts(3)
      print*, "Axies",nInts!,status,nfound

      call ftgkne(unit,'CDELT',StartVal,naxes,vReal,nfound,status)
      DC%DH%PixelSize=vReal(1:2)
      DC%DH%ChannelSize=vReal(3)
c      print*, "huh CDelt", vReal,nfound,status
      call ftgkne(unit,'CRPIX',1,naxes,vReal,nfound,status)
      DC%DH%RefLocation(0:2)=vReal(1:3)
      DC%DH%RefLocation(0:2)=DC%DH%RefLocation(0:2)-1.
c      print*, "huh CRPIX", vReal,nfound,status
c      DC%DH%PixelCenterIndx(0:1)=int(vReal(1:2))
c      DC%DH%ChannelCenterIndx=int(vReal(3))

      call ftgkne(unit,'CRVAL',1,naxes,vReal,nfound,status)
      DC%DH%RefVal(0:2)=vReal(1:3)
c      DC%DH%PixelCent=vReal(1:2)*3600.
c      DC%DH%ChannelCent=vReal(3)/1000.
c      print*, "huh CRVAL", vReal

      call ftgkns(unit,'CTYPE',1,naxes,DC%DH%AxisType,nfound,status)
c      print*, "Type?", nfound

      call ftgkns(unit,'CUNIT',1,naxes,DC%DH%Units(0:2),nfound,status)
c      print*, "units?",DC%DH%Units(0:2)
c           If CUNIT isn't found, set to defaults
      if(nfound .eq. 0) then
        print*, "Unit definitions not found.  Setting to "
     &          //"degrees and m/s (or Hz)"
        DC%DH%Units(0)="DEGREE"
        DC%DH%Units(1)="DEGREE"
        if(DC%DH%AxisType(2) .eq. "VELO-LSR"
     &          .or. DC%DH%AxisType(2) .eq. "VELO") then
            DC%DH%Units(2)="m/s"
        elseif(DC%DH%AxisType(2) .eq. "FREQ") then
            DC%DH%Units(2)="HZ"
        endif
      endif

c           If using frequency, load in the rest frequency
      if(DC%DH%AxisType(2) .eq. "FREQ") then
        call ftgkye(unit,'RESTFRQ',DC%DH%RestFreq,record,status)
      endif

      call FTGKYE(unit,'BMAJ',Beam%BeamMajorAxis,record,status)
      Beam%BeamMajorAxis=Beam%BeamMajorAxis
c      print*, "Beam Major Axis in '' ", Beam%BeamMajorAxis,status

      call FTGKYE(unit,'BMIN',Beam%BeamMinorAxis,record,status)
      Beam%BeamMinorAxis=Beam%BeamMinorAxis
c      print*, "Beam Minor Axis in '' ",Beam%BeamMinorAxis

      call FTGKYE(unit,'BPA',Beam%BeamPositionAngle,record,status)
c      print*, "Beam Position Angle Axis is '' "
c     &                  , Beam%BeamPositionAngle

      call FTGKYS(unit,'BUNIT',DC%DH%FUnit,record,status)
c      print*, "Data cube flux units",status,DC%DH%FUnit

      call FTGKYE(unit,'EPOCH',DC%DH%Epoch,record,status)
c      print*, "epoch"

c      call ftgkns(unit,'BTYPE',DC%DH%FType,record,status)
c      print*, "bType"

c      print*, "Datacube read Check", DC%DH%PixelSize,Beam%BeamMajorAxis
c     &              , Beam%BeamMajorAxis/DC%DH%PixelSize(0)

c      print*, "Done reading header"

      call DataCubeUnitConversions(DC,Beam)



C  The FITS file must always be closed before exiting the program.
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)


      end subroutine
C *************************************************************************








*************************************************************************
      subroutine readDCFile(filename,DC)
      implicit none

C  Read a FITS image and determine the minimum and maximum pixel value.
C  Rather than reading the entire image in
C  at once (which could require a very large array), the image is read
C  in pieces, 100 pixels at a time.

      Type(DataCube), INTENT(INOUT) :: DC
      integer status,unit,readwrite,blocksize,nfound
      integer i,j,k,l
      integer*8 group,firstpix,nbuffer,nCells,nPixels,naxes(3)
      real datamin,datamax,nullval,buffer(100)
      logical anynull
      character(*) filename
      real,ALLOCATABLE :: FlatArr(:)
      integer nNan,ltest,itest,jtest,ktest

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file previously created by WRITEIMAGE
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)


C  Initialize variables
      nCells=DC%DH%nPixels(0)*DC%DH%nPixels(1)
     &              *DC%DH%nChannels
      group=1
      firstpix=1
      nullval=-1010
c      nullval=0.
      datamin=1.0E30
      datamax=-1.0E30

      ALLOCATE(FlatArr(nCells))
      j=0
      nPixels=nCells
      do while (npixels .gt. 0)
C         read up to 100 pixels at a time
        nbuffer=min(100,npixels)

        call ftgpve(unit,group,firstpix,nbuffer,nullval,
     &            buffer,anynull,status)

C         find the min and max values
        do i=1,nbuffer
            j=j+1
c            datamin=min(datamin,buffer(i))
c            datamax=max(datamax,buffer(i))
            FlatArr(j)=buffer(i)
        end do

C         increment pointers and loop back to read the next group of pixels
        npixels=npixels-nbuffer
        firstpix=firstpix+nbuffer
      end do

c      print *
c      print *,'Min and max image pixels = ',datamin,datamax

      l=0
      nNan=0
      do k=0, DC%DH%nChannels-1
        do j=0, DC%DH%nPixels(1)-1
            do i=0, DC%DH%nPixels(0)-1
                l=l+1
c                   Do a nan check and sum up the total Nan cells
                if(FlatArr(l) .eq. nullval) then
                    nNan=nNan+1
                endif
                DC%Flux(i,j,k)=FlatArr(l)
            enddo
        enddo
      enddo
      
c       Allocate the number of valid cells
      DC%DH%nValid=nCells-nNan
c           This first requires deallocating the default allocation
      DEALLOCATE(DC%FlattendValidIndices)
      ALLOCATE(DC%FlattendValidIndices(0:DC%DH%nValid-1))
      

      l=0
      do i=0, DC%DH%nPixels(0)-1
        do j=0, DC%DH%nPixels(1)-1
            do k=0, DC%DH%nChannels-1
                call FlatIndxCalc(i,j,k,DC%DH,ltest)    !src/ObjectDefinitions/DataCube.f
    
                if(DC%Flux(i,j,k) .eq. nullval) then
                    DC%Flux(i,j,k)=0.
                else
                    DC%FlattendValidIndices(l)=ltest
                    l=l+1
c                    print*, "hmm",i,j,k,ltest,DC%Flux(i,j,k)
                endif
            enddo
        enddo
      enddo
      print*, "Data cube input test", sum(DC%Flux)



C  The FITS file must always be closed before exiting the program.
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
c       Also, deallocate the flattened array
      DEALLOCATE(FlatArr)

      end subroutine
ccccccccccccccccccccccccccccccc



ccccccc
c       This routine is meant to read the header for a map header
c           It is based on the fitsio routines.  It also includes the beam
c           object header input
      subroutine ReadMapHeader(FileName,DC)
      implicit none
      character(*),INTENT(IN) :: FileName
      Type(DataCube), INTENT(INOUT) :: DC

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character record*80

      integer naxes,nfound
      integer k

      integer nInts(2)
      real vReal(2)

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)


C     open the FITS file, with read-only access.  The returned BLOCKSIZE
C     parameter is obsolete and should be ignored.
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
      j = 0
100   continue
      j = j + 1


C  Get the data cube key words
      call ftghsp(unit,nkeys,nspace,status)


      naxes=2
      call ftgknj(unit,'NAXIS',1,naxes,nInts,nfound,status)
      DC%DH%nPixels=nInts(1:2)

      call ftgkne(unit,'CDELT',1,naxes,vReal,nfound,status)
      DC%DH%PixelSize=vReal(1:2)
c      print*, "huh CDelt", vReal,nfound,status
      call ftgkne(unit,'CRPIX',1,naxes,vReal,nfound,status)
      DC%DH%RefLocation(0:1)=vReal(1:2)
      DC%DH%RefLocation(0:1)=DC%DH%RefLocation(0:1)

      call ftgkne(unit,'CRVAL',1,naxes,vReal,nfound,status)
      DC%DH%RefVal(0:1)=vReal(1:2)

      call ftgkns(unit,'CTYPE',1,naxes,DC%DH%AxisType,nfound,status)

      call ftgkns(unit,'CUNIT',1,naxes,DC%DH%Units(0:1),nfound,status)
c      print*, "units?",DC%DH%Units(0:2)
c           If CUNIT isn't found, set to defaults
      if(nfound .eq. 0) then
        print*, "Unit definitions not found.  Setting to "
     &          //"degrees and m/s (or Hz)"
        DC%DH%Units(0)="DEGREE"
        DC%DH%Units(1)="DEGREE"
      endif


C  The FITS file must always be closed before exiting the program.
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)


      end subroutine
C *************************************************************************



*************************************************************************
      subroutine readImgAsDCSlice(filename,DCChan,DC)
      implicit none

C  Read a FITS image and determine the minimum and maximum pixel value.
C  Rather than reading the entire image in
C  at once (which could require a very large array), the image is read
C  in pieces, 100 pixels at a time.

      Type(DataCube), INTENT(INOUT) :: DC
      integer, INTENT(IN) :: DCChan
      integer status,unit,readwrite,blocksize,nfound
      integer i,j,k,l
      integer*8 group,firstpix,nbuffer,nCells,nPixels,naxes(3)
      real datamin,datamax,nullval,buffer(100)
      logical anynull
      character(*), INTENT(IN) :: filename
      real,ALLOCATABLE :: FlatArr(:)

      real sumTest

      print*, "Loading moment", trim(filename)

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file previously created by WRITEIMAGE
      readwrite=0
      call ftopen(unit,trim(filename),readwrite,blocksize,status)
      print*, "file opened", status


C  Initialize variables
      nCells=DC%DH%nPixels(0)*DC%DH%nPixels(1)
      group=1
      firstpix=1
      nullval=-999
      datamin=1.0E30
      datamax=-1.0E30
      sumTest=0.

      ALLOCATE(FlatArr(nCells))
      j=0
      nPixels=nCells
      do while (npixels .gt. 0)
C         read up to 100 pixels at a time
        nbuffer=min(100,npixels)
        call ftgpve(unit,group,firstpix,nbuffer,nullval,
     &            buffer,anynull,status)

C         find the min and max values
        do i=1,nbuffer
            j=j+1
            FlatArr(j)=buffer(i)
            sumTest=sumTest+FlatArr(j)
        end do
C         increment pointers and loop back to read the next group of pixels
        npixels=npixels-nbuffer
        firstpix=firstpix+nbuffer
      end do
      l=0
      k=DCChan  !Load the flat array to the target datacube channel
      do j=0, DC%DH%nPixels(1)-1
        do i=0, DC%DH%nPixels(0)-1
            l=l+1
            DC%Flux(i,j,k)=FlatArr(l)
        enddo
      enddo
      print*, "Sum Test", sumTest, sum(DC%Flux), nCells


C  The FITS file must always be closed before exiting the program.
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

      return
      end subroutine
ccccccccccccccccccccccccccccccc



      end module
