ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This is the main routine for the the WALLABY Resolved Kinematic Pipeline.
c   It is open source, but requires the openmpi and fftw3 libraries.
c
c   The code is meant to take in a catalogue of resolved detections and
c   generate a rotating disk model for each object.
c
c   The code is focused on low-resolution images, but it is meant to be
c   flexible and easy for the user to modify.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program WallabyKinematicPipeline

      use TiltedRingMod
      use SingleRingGenerationMod
      use TiltedRingGenerationMod
      use TiltedRingOutputsMod
      use DataCubeMod
      use FillDataCubeWithTiltedRingMod
      use DataCubeOutputsMod
      use BeamMod
      use CalcBeamKernelMod
      use TwoDConvolutionMod
      use CubeKernelConvolutionMod
      use CubeCompareMod
      use LikelihoodMod
      use ParameterVectorMod
      use ParameterVectorToTiltedRingMod
      use MaskCubeMod

      use DataCubeInputMod
      use FullGalaxyInputMod

      use PipelineGlobals
      use FullModelComparisonMod

      use DownhillSimplexMod

      use MakeMomentMapsMod

      use PreGalaxyAnalysisMod

      use SofiaInputMod

      use PipelineInputMod

      use GalaxyFitMod
      use BBaroloFitMod

      use FitOutputMod

      use MPISplitMod
      use mpi


      implicit none

      Type(TiltedRingModel) TR1,TR2
      Type(TiltedRingFittingOptions) TRFO1

      Type(DataCube) TestDC,TestDC2,TestMask

      Type(Beam2D) Test2DBeam

      Type(ParameterVector) PV1, PV2

      real Inclination, PositionAngle

      real RotationMatrix(0:2,0:2)

      real pixelarea

      integer i,j, k

      character(100) TestFile

      integer SizeDC(2),SizeK(2),SizeP(2)


      real chi2

      integer nRingFree,nRingFixed

      integer TestArrSize(2)
      real,ALLOCATABLE :: TestArr(:,:)


      real,ALLOCATABLE :: paramGuesses(:,:),chiArray(:)
      integer iter
      real lambda


      integer ierr,rc
c      integer nProcessors,PID
c      integer Primary
c      parameter(Primary = 0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "WallabyKinematicPipeline"

c       Start by setting up the processors
      call MPI_INIT(ierr)
      if (ierr .ne. MPI_SUCCESS) then
        print *,'Error starting MPI program. Terminating.'
        call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
      end if
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcessors, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, PID, ierr)

c       Read in the main pipeline inputs and set the points
      call PipelineIn(PID,Primary)             !Inputs/PipelineRuntimeInputs.f
c       Have the primary processor read in the full catalogue
      if(PID .eq. Primary) then
        call GeneralSoFiAIn()       !Inputs/SofiaInputs.f
      endif

c       Split the catalogue across processors
      call SplitCatalogue(SCatalogue,SCatLocal,Primary,PID,nProcessors) !MPIRoutines/MPISplit.f

      GalaxyFit=>GalaxyFit_Simple
      OutputFit=>OutputBestFit_Simple
c      GalaxyFit=>BBaroloFit


c           Loop through all local catalogue items
c      print*, "number of objects", PID,SCatLocal%nObjects
c      do i=1,1
c      do i=0,SCatLocal%nObjects-1
      do i=0,0
        print*, "Object in Loop",PID, i, SCatLocal%Objects(i)%ObID
     &          ,SCatLocal%nObjects
cc           Read in the specific galaxy
        call FullGalaxyIn(i)        !Inputs/FullGalaxyInputs.f
cc        print*, ObservedBeam%BeamMajorAxis
c        print*, "After inputs test", PID,i,sum(ObservedDC%Flux)
        call PreGalaxyAnalysis(SCatLocal%Objects(i))
cc        print*, ObservedBeam%BeamMajorAxis
ccc        print*, "After pre-analysis test", PID,i,sum(ObservedDC%Flux)
        call GalaxyFit(SCatLocal%Objects(i))
cc        print*, ObservedBeam%BeamMajorAxis
cc        print*, "After fit test", PID,i,sum(ObservedDC%Flux)
        call OutputFit(SCatLocal%Objects(i))

c        print*, "Done Main and about to deallocate"
c        call DeAllocateDataCube(ObservedDC)
cc        call DeAllocateDataCube(ModelDC)
c        call DeAllocateDataCube(DataCubeMask)
c        call DeAllocate_Beam2D(ObservedBeam)
c        call DeAllocateDataCube(ObservedMaps)
c        call TiltRing_DeAllocateStruct(ModelTiltedRing)
c        call TiltRingFittingOptions_DeAllocate(TR_FittingOptions)
c        call DeAllocateParamVector(PVIni)
ccc        call DeAllocateParamVector(PVModel)
c        print*, " finished step", i, " in PID", PID
c        print*, " "

      enddo
      print*, "DOne Loop", PID,i,SCatLocal%nObjects-1
      print*, "MPI status", PID,ierr


      call MPI_FINALIZE(ierr)
      print*, "MPI status", PID,ierr


      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




