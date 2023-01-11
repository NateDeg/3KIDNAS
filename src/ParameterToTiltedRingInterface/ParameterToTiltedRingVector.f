cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
c       All the routines assume that some basic parameters
c       of the tilted ring model have been set:
c           nRings, cmode, CloudBaseSurfDens, Rmid, Rwidth
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ParameterVectorToTiltedRingMod
      use TiltedRingMod
      use ParameterVectorMod

      implicit none

      Procedure(P_To_TR_Interface),POINTER:: ParamToTiltedRing=>null()

      ABSTRACT INTERFACE
        subroutine P_To_TR_Interface(PV,TR,TRFO)
            import :: ParameterVector
            import :: TiltedRingModel
            import :: TiltedRingFittingOptions
            implicit none
            Type(ParameterVector), INTENT(IN) :: PV
            Type(TiltedRingModel), target, INTENT(INOUT) :: TR
            Type(TiltedRingFittingOptions), target,INTENT(IN) :: TRFO
        end subroutine P_To_TR_Interface
      END INTERFACE


      contains



cccccccc
c       This routine does a generalized conversion of a parameter vector to a tilted ring model
      subroutine GeneralizedParamVectorToTiltedRing(PV,TR,TRFO)
      implicit none
      Type(ParameterVector), INTENT(IN) :: PV
      Type(TiltedRingModel), target,INTENT(INOUT) :: TR
      Type(TiltedRingFittingOptions), target, INTENT(IN) :: TRFO
      integer i, CurrParam
      real,pointer,dimension(:) :: TempVector,VariableVector

c      print*, "in generalize param vector"

      CurrParam=0
      ALLOCATE(TempVector(0:TR%nRings-1))
      ALLOCATE(VariableVector(0:TR%nRings-1))

c      print*, "Param Vec", PV%Param(0:PV%nParams-1), PV%nParams

c           Loop through all the specifc parameters
      do i=0,12
c           Associate the temporary pointers with the correct array
        if(i .eq.0) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(0)
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%CentPos(0)
        elseif(i .eq.1) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(1)
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%CentPos(1)
        elseif(i .eq.2) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Inclination
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%Inclination
        elseif(i .eq.3) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%PositionAngle
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%PositionAngle
        elseif(i .eq.4) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VSys
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VSys
        elseif(i .eq.5) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRot
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VRot
        elseif(i .eq.6) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRad
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VRad
        elseif(i .eq.7) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VDisp
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VDisp
        elseif(i .eq.8) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Vvert
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%Vvert
        elseif(i .eq.9) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%dvdz
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%dvdz
        elseif(i .eq.10) then
            VariableVector=>
     &                TRFO%RadialProfiles(0:TRFO%nRings-1)%SigUse
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%SigUse

        elseif(i .eq.11) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%z0
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%z0
        elseif(i .eq.12) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%zGradiantStart
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%zGradiantStart
c
        endif
c               Pass all the vectors to the vector assignment
        call SetSpecificVector(PV%nParams,TR%nRings,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%Param,VariableVector,TempVector)

c        print*, "Fit Check",TRFO%ConstParams(i),TRFO%FixedParams(i)
c        print*, "vari", VariableVector
c        print*, i, TempVector

c        DEALLOCATE(VariableVector)
c        DEALLOCATE(TempVector)
      enddo
c       Calculate the various important derivatives on both sides of each point
c      call TiltedRingDerivativeCalc(TR)



      return
      end subroutine
ccccccc

cccccc
c           This routine sets a specifc vector based on the fitting options
      subroutine SetSpecificVector(nParam,nRings,Const,Fixed,CurrParam
     &                  ,Param,AlternateSourceVector
     &                  ,TargVec)
      implicit none
      integer,INTENT(IN) :: nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(IN) :: Param(0:nParam-1)
      real, INTENT(IN) :: AlternateSourceVector(0:nRings-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(INOUT) :: TargVec(0:nRings-1)

      integer i

c      print*, "Setting Vec", CurrParam, Param(Currparam),Const,Fixed
      if(Const) then
        if(Fixed) then
            TargVec(0:nRings-1)=AlternateSourceVector(0)
        else
            TargVec(0:nRings-1)=Param(CurrParam)
            CurrParam=CurrParam+1
        endif
      else
        if(Fixed) then
            TargVec(0:nRings-1)=AlternateSourceVector(0:nRings-1)
        else
            TargVec(0:nRings-1)=Param(CurrParam:CurrParam+nRings-1)
            CurrParam=CurrParam+nRings
        endif
      endif
c      print*,TargVec(0:nRings-1)


      return
      end subroutine
cccccccc




ccccccc
c       This routine does a generalized conversion of a parameter vector to a tilted ring model
      subroutine TiltedRingOptionsToPV(PV,TRFO)
      implicit none
      Type(ParameterVector), INTENT(INOUT) :: PV
      Type(TiltedRingFittingOptions), target, INTENT(IN) :: TRFO
      integer i, CurrParam
      real,pointer,dimension(:) :: VariableVector

c      print*, "in generalize param vector"

      CurrParam=0
      ALLOCATE(VariableVector(0:TRFO%nRings-1))
c           Loop through all the specifc parameters
      do i=0,12
c           Associate the temporary pointers with the correct array
        if(i .eq.0) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(0)
        elseif(i .eq.1) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(1)
        elseif(i .eq.2) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Inclination
        elseif(i .eq.3) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%PositionAngle
        elseif(i .eq.4) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VSys
        elseif(i .eq.5) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRot
        elseif(i .eq.6) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRad
        elseif(i .eq.7) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VDisp
        elseif(i .eq.8) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Vvert
        elseif(i .eq.9) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%dvdz
        elseif(i .eq.10) then
            VariableVector=>
     &                 TRFO%RadialProfiles(0:TRFO%nRings-1)%SigUse
        elseif(i .eq.11) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%z0
        elseif(i .eq.12) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%zGradiantStart
        endif
c               Pass all the vectors to the vector assignment

c           Pass the lower limits to the parameter vector object
        call SetParamLimsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%ParamLowerLims
     &          ,TRFO%ParamLowerLims)

c           Pass the upper limits to the parameter vector object
        call SetParamLimsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%ParamUpperLims
     &          ,TRFO%ParamUpperLims)

c           Pass the parameter ranges to the parameter vector object
        call SetParamLimsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%ParamRange
     &          ,TRFO%ParamRange)



      call SetParamCyclicsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%CyclicSwitch
     &          ,TRFO%CyclicSwitch)

        call SetParamFromFittingOptions(PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%Param,VariableVector)

c        print*, "Fit Check Ini",i,TRFO%ConstParams(i)
c     &                  ,TRFO%FixedParams(i)
c        print*, VariableVector

      enddo
c      print*, "Final Parameter vector", PV%Param(0:PV%nParams-1)
c     &              ,PV%nParams


      return
      end subroutine
ccccccc


cccccc
c           This routine steps through a parameter vector and sets the
c           values from a target vector based on the logical switchs
      subroutine SetParamLimsFromFittingOptions(ParamNum,nParam,nRings
     &                  ,Const,Fixed,CurrParam
     &                  ,ParamLims
     &                  ,TargLims)
      implicit none
      integer,INTENT(IN) :: ParamNum,nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(INOUT) :: ParamLims(0:nParam-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(IN) :: TargLims(0:12)

      integer i

      if(Fixed .eqv. .False.) then
        if(Const) then
            ParamLims(CurrParam)=TargLims(ParamNum)
        else
            ParamLims(CurrParam:CurrParam+nRings-1)=TargLims(ParamNum)
        endif
      endif
      

      return
      end subroutine
cccccccc




cccccc
c           This routine steps through a parameter vector and sets the
c           values from a target vector based on the logical switchs
      subroutine SetParamFromFittingOptions(nParam,nRings
     &                  ,Const,Fixed,CurrParam
     &                  ,Param
     &                  ,TargVec)
      implicit none
      integer,INTENT(IN) :: nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(INOUT) :: Param(0:nParam-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(INOUT) :: TargVec(0:nRings-1)

      integer i

      if(Fixed .eqv. .False.) then
        if(Const) then
            Param(CurrParam)=TargVec(0)
            CurrParam=CurrParam+1
        else
            Param(CurrParam:CurrParam+nRings-1)=TargVec(0:nRings-1)
            CurrParam=CurrParam+nRings
        endif
      endif

      return
      end subroutine
cccccccc





cccccc
c           This routine steps through a parameter vector and sets the
c           cyclic switches
      subroutine SetParamCyclicsFromFittingOptions(ParamNum,nParam
     &              ,nRings
     &                  ,Const,Fixed,CurrParam
     &                  ,ParamCyclics
     &                  ,TargCyclics)
      implicit none
      integer,INTENT(IN) :: ParamNum,nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      integer,INTENT(INOUT) :: ParamCyclics(0:nParam-1)
      integer,INTENT(INOUT) :: CurrParam
      integer,INTENT(IN) :: TargCyclics(0:12)

      integer i


      if(Fixed .eqv. .False.) then
          if(Const) then
            ParamCyclics(CurrParam)=TargCyclics(ParamNum)
          else
            ParamCyclics(CurrParam:CurrParam+nRings-1)
     &                  =TargCyclics(ParamNum)
          endif
      endif

      return
      end subroutine
cccccccc


      end module
