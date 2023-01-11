ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This is the main routine for a code to generate multiple bootstrap
c       samples from some input cube.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program CubeBootStrapGenerator
      use BootstrapInputMod
      use BootstrapGlobals
      use DataCubeInputMod
      use GenBootstrapMod
      use DataCubeMod
      use DataCubeOutputsMod
      use CalcBeamKernelMod
      use BeamMod

      use CubeDiffMod
      use DataCubeOutputsMod

      use GenBootstrapMod

      use FlippingBootstrapMod


      implicit none

      character(500) OutName
      integer i

      real IncTemp,PATemp,Center(2),VSysTemp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "Boot Strap Sampler"
c           Get the runtime inputs
      call BootstrapIn()
c       Load in the relevant cubes
      call LoadCubesForBootstrap()

c       Make a bootstrap sample
c      call GenBootstrapSample()
      call GenFlipBootstrapSample()


c       Output the resampled cube
      BootstrapCube%DH%PixelSize=
     &              BootstrapCube%DH%PixelSize/3600.
c      BootstrapCube%DH%RefLocation(0:2)=
c     &                  BootstrapCube%DH%RefLocation(0:2)+1

      SampleFile=trim(BaseOutName)//".fits"
      call WriteDataCubeToFITS(BootstrapCube
     &          ,ObservedBeam,SampleFile
     &          ,trim(BaseOutName))


      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



