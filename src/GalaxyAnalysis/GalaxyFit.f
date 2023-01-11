cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the general procedures for modeling
c       a galaxy once the pre-analysis is completed.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GalaxyFitMod
      use ParameterVectorMod
      use DataCubeMod
      use SoFiACatalogueMod
      use BeamMod
      use CalcBeamKernelMod
      use DownhillSimplexMod
      use FullModelComparisonMod
      implicit none

      PROCEDURE(GeneralFitInterface),POINTER :: GalaxyFit =>null()

      ABSTRACT INTERFACE
        subroutine GeneralFitInterface(CatItem)
            import :: CatalogueItem
            Type(CatalogueItem),INTENT(IN) :: CatItem
        END subroutine GeneralFitInterface
      END INTERFACE

      contains
cccccccc
c
      subroutine GalaxyFit_Simple(CatItem)
      use PipelineGlobals
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      real,ALLOCATABLE :: paramGuesses(:,:),chiArray(:)
      real chi2
      integer i,iter,StrictEstimate

      logical FitFlag

      character(100) FullEnsambleOutput

      print*, "Fitting Galaxy",PID

      print*, "Initial PV"
      do i=0, PVIni%nParams-1
        print*, i,PVini%Param(i)
      enddo
c
c       Set up the beam
c
c      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
      call Calculate2DBeamKernel(ObservedBeam,ObservedDC%DH%PixelSize)

c       Allocate a model datacube with the same dimensions as the observed cube
      ModelDC%DH=ObservedDC%DH

      call AllocateDataCube(ModelDC)
c       Copy the initial parameter guess into the model parameter vector
      PVModel%nParams=PVIni%nParams
      call AllocateParamVector(PVModel)
      PVModel%Param=PVIni%Param
      PVModel%ParamLowerLims=PVIni%ParamLowerLims
      PVModel%ParamUpperLims=PVIni%ParamUpperLims
      PVModel%CyclicSwitch=PVIni%CyclicSwitch
      PVModel%ParamRange=PVIni%ParamRange

c           Set up the array of parameter guesses and chi^2 values needed
      ALLOCATE(paramGuesses(PVModel%nParams+1,
     &          PVModel%nParams))
      ALLOCATE(chiArray(PVModel%nParams+1))
c           Create an array of parameter guesses


      call TiltedRingModelComparison(PVModel%Param,chi2)
      print*, "Initial model fit", chi2
      IniGuessWidth=1.
      ftol=0.005

c           Make the set of initial guesses
      StrictEstimate=1
      call MakeParamGuessArray(PVModel,ParamGuesses
     &          ,PVIni%nParams,idum,IniGuessWidth,StrictEstimate)

c      FullEnsambleOutput="AllInitialParamVectors.txt"
c      open(10,file=FullEnsambleOutput,status='replace')
c      do i=1, PVModel%nParams+1
c        write(10,*) ParamGuesses(i,1:PVModel%nParams)
c      enddo
c      close(10)

      call DownhillSimplexRun(PVModel%nParams
     &                  ,paramGuesses,chiArray)


      PV_FirstFit%nParams=PVIni%nParams
      call AllocateParamVector(PV_FirstFit)
      PV_FirstFit%Param(0:PVModel%nParams-1)=
     &              ParamGuesses(1,1:PVModel%nParams)
      PV_FirstFit%BestLike=chiArray(1)

      IniGuessWidth=0.5
      ftol=ftol/5.

c           Make the set of initial guesses
      StrictEstimate=0
      call MakeParamGuessArray(PVModel,ParamGuesses
     &          ,PVIni%nParams,idum,IniGuessWidth,StrictEstimate)

c      print*, "PV Model",PVModel%Param

      call DownhillSimplexRun(PVModel%nParams
     &                  ,paramGuesses,chiArray)



c      print*, PVModel%Param
c           TEMPORARY WORK
c           convert the best model to a TR model
      call ParamToTiltedRing(PVModel,ModelTiltedRing
     &          ,TR_FittingOptions)
c      do i=0, ModelTiltedRing%nRings-1
c        print*, i, ModelTiltedRing%R(i)%VRot
c     &              , ModelTiltedRing%R(i)%VDisp
c     &              , ModelTiltedRing%R(i)%Inclination*180./Pi
c     &              , ModelTiltedRing%R(i)%PositionAngle*180./Pi
c      enddo

c      FullEnsambleOutput="AllFinalParamVectors.txt"
c      open(10,file=FullEnsambleOutput,status='replace')
c      do i=1, PVModel%nParams+1
c        write(10,*) ParamGuesses(i,1:PVModel%nParams)
c      enddo
c      close(10)


c       Deallocate the model vector at the end
c      call DeAllocateParamVector(PVModel)
      DEALLOCATE(chiArray)
      DEALLOCATE(paramGuesses)
      return
      end subroutine
ccccccccc



ccccccc
c
      subroutine DownhillSimplexRun(nParams,paramGuesses,chiArray)
      use PipelineGlobals
      implicit none
      integer,INTENT(IN) :: nParams
      real,INTENT(INOUT):: paramGuesses(nParams+1,nParams)
      real,INTENT(INOUT) :: chiArray(nParams+1)
      integer i,iter
      real chi2
      logical FitFlag


      print*, "ftol ChecK", ftol, ObservedDC%DH%Uncertainty

c           Make the set of initial guesses
c      call MakeParamGuessArray(PVModel,ParamGuesses
c     &          ,PVIni%nParams,idum,IniGuessWidth)

c       Get the goodness of fit for each guess
      chiArray=0.
      do i=1,PVModel%nParams+1
        PVModel%Param(0:PVModel%nParams-1)=
     &                  ParamGuesses(i,1:PVModel%nParams)
        call TiltedRingModelComparison(PVModel%Param,chi2)
        chiArray(i)=chi2
        print*, PID,i,paramGuesses(i,:),chi2
      enddo

c       Run the downhill simplex
      call amoeba(paramGuesses,chiArray
     &                  ,PVModel%nParams+1
     &                  ,PVModel%nParams
     &                  ,PVModel%nParams,ftol
     &                  ,TiltedRingModelComparison,iter,PID
     &                  ,FitFlag)

c       Store the best model in the PVModel object
      PVModel%BestLike=chiArray(1)
      PVModel%Param(0:PVModel%nParams-1)=
     &              ParamGuesses(1,1:PVModel%nParams)


      return
      end subroutine
cccccccc

ccccccc
      subroutine MakeParamGuessArray(PredictedPV,ParamGuesses
     &                      ,ndim,idum,lambda,StrictEstimate)

      use BasicRanNumGen
      implicit none

      integer idum
      integer ndim
      Type(ParameterVector) PredictedPV
      real ParamGuesses(ndim+1,ndim)
      real lambda,lambdaPar

      integer i,j,k
      integer counter
      integer StrictEstimate,ParCounter,AcceptBeyondLimits

      print*, "Param Limits Check",idum
      do i=0, ndim-1
        print*, i, PredictedPV%ParamLowerLims(i)
     &          ,PredictedPV%ParamUpperLims(i)
     &          ,PredictedPV%Param(i)
     *          ,PredictedPV%ParamRange(i)
      enddo

      counter=0
      print*, "Param Guess Array Creation"
      do k=1,ndim+1
        if(k .eq. 1) then
            ParamGuesses(k,1:ndim)=PredictedPV%Param(0:ndim-1)
        else
            do i=1,ndim
                ParCounter=0
                AcceptBeyondLimits=0
                j=i-1
c100             lambdaPar=lambda*(PredictedPV%ParamUpperLims(j)
c     &                      -PredictedPV%ParamLowerLims(j))
 100            lambdaPar=lambda*PredictedPV%ParamRange(j)
                lambdaPar=(2*ran2(idum)-1.)*lambdaPar
                ParamGuesses(k,i)=PredictedPV%Param(j)+lambdaPar
                counter=counter+1
                ParCounter=ParCounter+1
c                   If the strict estimate is relaxed, accept the parameter regardless of the limits after 200 tries
                if(StrictEstimate .eq. 0) then
                    if(ParCounter .ge. 200) AcceptBeyondLimits=1
                endif
c       print*, i,ParamGuesses(k,i),PredictedPV%ParamLowerLims(j)
c     &                  ,PredictedPV%ParamUpperLims(j),counter
c     &              ,j,PredictedPV%CyclicSwitch(j)
c               For cyclic parameters, adjust them to be inside the range
200             continue
                if(PredictedPV%CyclicSwitch(j) .eq. 1) then
                    if (ParamGuesses(k,i) .lt.
     &                   PredictedPV%ParamLowerLims(j)) then
                            ParamGuesses(k,i)=ParamGuesses(k,i)
     &                          +PredictedPV%ParamUpperLims(j)
                        goto 200
                    elseif(ParamGuesses(k,i) .gt.
     &                   PredictedPV%ParamUpperLims(j)) then
                            ParamGuesses(k,i)=ParamGuesses(k,i)
     &                          -PredictedPV%ParamUpperLims(j)
                        goto 200
                    endif
                endif

                if(AcceptBeyondLimits .eq. 0) then
                    if(ParamGuesses(k,i) .lt.
     &                   PredictedPV%ParamLowerLims(j)) goto 100
                    if(ParamGuesses(k,i) .gt.
     &                   PredictedPV%ParamUpperLims(j)) goto 100
                endif
            enddo
        endif
c        print*, k,ParamGuesses(k,:)
      enddo
      print*, "Done param guess array creation"

      return
      end subroutine
cccccc

      end module
