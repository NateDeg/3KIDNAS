cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading a
c         SoFiA catalogue file
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FittingOptionsMod

      use PipelineGlobals
      use LikelihoodMod
      use ParameterVectorToTiltedRingMod
      use GalaxyFitMod
      use BBaroloFitMod
      use FitOutputMod

      implicit none

      contains

ccccccc
c           This routine selects the specific format of the SoFiA catalogue file
      subroutine FittingOptionsIn()
      implicit none
      integer ParamInterfaceSwitch,i
      character(20),ALLOCATABLE:: RadStr(:)

      print*, "Reading in the various fitting options"

      open(10,file=FittingOptionsFile,status='old')

c       Get the switch for the algorithm engine
      read(10,*)
      read(10,*) PFlags%CoreCodeSwitch
c           Set the fitting algorithm
      if(PFlags%CoreCodeSwitch .eq. 1) then
        GalaxyFit=>GalaxyFit_Simple
        OutputFit=>OutputBestFit_Simple
      elseif(PFlags%CoreCodeSwitch .eq. 2) then
        GalaxyFit=>BBaroloFit
        OutputFit=>OutputBestFit_Simple
      else
        print*, "No valid code driver selected"
        stop
      endif

c       Get the switch for the likelihood function
      read(10,*)
      read(10,*) PFlags%LikelihoodSwitch
c           Set the likelihood function pointer
      if(PFlags%LikelihoodSwitch .eq. 1) then
        LikePoint => Chi2Calc
      elseif(PFlags%LikelihoodSwitch .eq. 2) then
        LikePoint => Chi2Calc_logElements
      elseif(PFlags%LikelihoodSwitch .eq. 3) then
        LikePoint => LogChi2Calc
      else
        print*, "No valid likelihood function selected"
        stop
      endif

c       Get the switch for the parameter interface switch
      read(10,*)
      read(10,*) ParamInterfaceSwitch
c           Set the interface pointer switch
      if(ParamInterfaceSwitch .eq. 1) then
        ParamToTiltedRing=> GeneralizedParamVectorToTiltedRing
c      elseif(ParamInterfaceSwitch .eq. 2) then
c        ParamToTiltedRing=> SimpleFlatDiskParamToTiltedRing
      else
        print*, "No valid parameter-tilted ring interface selected"
        stop
      endif


c       Get the switch for how the moment maps are obtained
      read(10,*)
      read(10,*) MomentMapSwitch
c           Check if the center switch is valid
      if(MomentMapSwitch .eq. 3) then
        read(10,*)
        read(10,'(a)') SecondaryCubeFile
        read(10,'(a)') SecondaryMaskFile
      elseif(MomentMapSwitch .eq. 2) then
        read(10,*)
        read(10,'(a)') MapNames(0)
        read(10,'(a)') MapNames(1)
        read(10,'(a)') MapNames(2)
        read(10,'(a)') ProfileName
      endif
      if(MomentMapSwitch .ge. 4
     &          .or. MomentMapSwitch .le. -1) then
        print*, "No valid method of obtaining "
     &          //"moment maps selected"
        stop
      endif

c       Get the switch for how the center of the image estimate is to be found
      read(10,*)
      read(10,*) PFlags%CenterSource
c           Check if the center switch is valid
      if(PFlags%CenterSource .ge. 3
     &          .or. PFlags%CenterSource .le. -1) then
        print*, "No valid method of estimating the "
     &          //"image centers selected"
        stop
      endif

c       Get the switch for how the image shape estimate is found
      read(10,*)
      read(10,*) PFlags%ShapeSource
c           Check if the center switch is valid
      if(PFlags%ShapeSource .ge. 2
     &          .or. PFlags%ShapeSource .le. -1) then
        print*, "No valid method of estimating the "
     &          //"image shape selected"
        stop
      endif

c       Get the switch for how the estimated image size is found
      read(10,*)
      read(10,*) PFlags%SizeSource
c           Check if the center switch is valid
      if(PFlags%SizeSource .ge. 2
     &          .or. PFlags%SizeSource .le. -1) then
            print*, "No valid method of estimating the "
     &          //"image size selected"
        stop
      endif


c       Get the switch for how the systemic velocity is obtained
      read(10,*)
      read(10,*) PFlags%VSysSource
c           Check if the center switch is valid
      if( PFlags%VSysSource.ge. 2
     &          .or. PFlags%VSysSource .le. -1) then
        print*, "No valid method of obtaining "
     &          //"the systemic velocity selected"
            stop
      endif

c       Get the switch whether using linear or logarithmic surface density
      read(10,*)
      read(10,*) PFlags%Linear_Log_SDSwitch
c           Check if the switch is valid
      if( PFlags%Linear_Log_SDSwitch.ge. 2
     &          .or. PFlags%Linear_Log_SDSwitch .le. -1) then
        print*, "No valid linear/logarithmic "
     &          //"fitting option selected."
        stop
      endif

c       Get the cloud mode being used
      read(10,*)
      read(10,*) ModelTiltedRing%cmode
      if(ModelTiltedRing%cmode .gt. 2
     &          .or. ModelTiltedRing%cmode .le. -1) then
        print*, "No valid tilted ring cloud "
     &          //"mode selected"
        stop
      endif



c       Get the base cloud surface density
      read(10,*)
      read(10,*) ModelTiltedRing%CloudBaseSurfDens

c       Get the number of sigma's to use for the beam kernal calculation/convolution
      read(10,*)
      read(10,*) ObservedBeam%SigmaLengths

c       Get the noise limit to reach for figuring out the number of rings
      read(10,*)
      read(10,*) NoiseSigmaLim

c       Get the total number of rings to use.  If the total number of rings is -1, then calculate the number automatically
      read(10,*)
      read(10,*) TR_FittingOptions%nTargRings
      if (TR_FittingOptions%nTargRings .lt. -1
     &              .or. TR_FittingOptions%nTargRings .eq. 0) then
        print*, "The number of rings must be either -1 for automated"
     &      //" calculation or a positive number"
        print*, "Current value",TR_FittingOptions%nTargRings
        stop
      endif
c       If the number of rings is specified, read in the specific radii to be used
      if(TR_FittingOptions%nTargRings .gt. 0) then
c           Start by allocating the ring and string arrays
        ALLOCATE(RadStr(TR_FittingOptions%nTargRings))
        ALLOCATE(RadGrid(TR_FittingOptions%nTargRings))
        read(10,*)
c           Read all the radii as strings
        read(10,*) RadStr(1:TR_FittingOptions%nTargRings)
c           Convert the strings to reals
        do i=1,TR_FittingOptions%nTargRings
            read(RadStr(i) , '(f10.3)') RadGrid(i)
        enddo
        DEALLOCATE(RadStr)
      endif

c       Get the number of rings/beam
      read(10,*)
      read(10,*) TR_FittingOptions%nRingsPerBeam
      if(TR_FittingOptions%nRingsPerBeam .lt. 1
     &          .or. TR_FittingOptions%nRingsPerBeam .ge. 5) then
        print*, "No valid number of rings/beam "
     &          //"selected"
        stop
      endif


c       Read in all the tilted ring fitting options
      do i=0, 12
        read(10,*)
        read(10,*) TR_FittingOptions%ConstParams(i)
     &                  ,TR_FittingOptions%FixedParams(i)

      enddo


      close(10)

      return
      end subroutine
cccccccc






      end module
