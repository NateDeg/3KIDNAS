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

      module InputUnitConversionsMod

      use DataCubeMod
      use BeamMod
      use TiltedRingMod
      use CalcBeamKernelMod

      use CommonConsts
      use UnitConvertMod

      implicit none

      contains



ccccc
c     This routine converts angles to radians
      subroutine GeneralAngularConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to radians
        call DegreesToRadians(A,A)
      endif

      return
      end subroutine
cccccccc

cccccccc
c       This routine converts distances (angular) to arcseconds
      subroutine GeneralDistanceConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to arcseconds
        call DegreesToArcSec(A,A)
      endif
      return
      end subroutine
cccccccc

cccccc
c       This routine does velocity conversions
      subroutine GeneralVelocityConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

c      print*, "Vel Conversion", A, Switch
      if(Switch .eq. 0) then        !0=m/s to km/s
        call MToKM(A,A)
      endif
      return
      end subroutine
ccccccc
ccccccc
c     This routine does brightness conversion
      subroutine GeneralBrightnessConversion(Switch,A,ChannelSize
     &              ,BeamArea,PixelSize)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A
      real, INTENT(IN) :: ChannelSize,BeamArea,PixelSize

c       The surface brightness needs to be converted to Jy pixel^2
      if(Switch .eq. 0) then        !0=Jy km/s arcsec^-2 to Jy arcsec^-2
        A=A/abs(ChannelSize)
        A=A*abs(PixelSize)**2.           !Then Jy arcsec^2 to Jy pixel
      elseif(Switch .eq. 1) then    !1=mJy/beam to Jy/Beam
        A=A/1000.
      elseif(Switch .eq. 2) then    !2=Jy arcsec^-2 to Jy/Beam
        A=A*BeamArea
      endif
      return
      end subroutine
cccccc








ccccc
c
c       This routine converts the tilted ring
c           center coordinates to a pixel value
c
c       For this routine to work, StartVal must be the
c       pixel location (0) and have units of arcseconds
      subroutine PixelPositionConversion(Switch,A
     &                  ,PixelSize,StartVal)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(IN) :: StartVal,PixelSize
      real, INTENT(INOUT) :: A
c
      if(Switch .eq. 2) then        !2==pixels already
        A=A
      elseif(Switch .eq. 1) then    !1==arcseconds to pixels
        A=(A-StartVal)/PixelSize
      elseif(Switch .eq. 0) then    !0==degrees to pixels
        A=A*3600.
        A=(A-StartVal)/PixelSize
      endif

      return
      end subroutine
ccccccc


cccccc
c       This routine does the angular conversions for tilted ring parameters
      subroutine TR_UnitConversions(TR,ChannelSize,BeamArea
     &                  ,DC)
      implicit none
      Type(TiltedRingModel),INTENT(INOUT) :: TR
      Type(DataCube),INTENT(IN) :: DC
      real, INTENT(IN) :: ChannelSize,BeamArea
      integer i
c       Loop through all rings
      do i=0, TR%nRings-1
c       Do the central position conversions first
        call PixelPositionConversion(TR%UnitSwitchs(0)
     &          ,TR%R(i)%CentPos(0)
     &          ,DC%DH%PixelSize(0),DC%DH%Start(0))
        call PixelPositionConversion(TR%UnitSwitchs(1)
     &          ,TR%R(i)%CentPos(1)
     &          ,DC%DH%PixelSize(1),DC%DH%Start(1))
c       Now do the angular conversions
c       Do the geometric conversions first
        call GeneralAngularConversion(TR%UnitSwitchs(1)
     &          ,TR%R(i)%Inclination)
        call GeneralAngularConversion(TR%UnitSwitchs(1)
     &          ,TR%R(i)%PositionAngle)
c           Once the position angle is in radians, add a 90^degree rotation to get the standard astronomy
c               definition for position angle (North-South with approaching velocity in the North)
c                   Note that in outputs we'll need to undo this rotation to record the position angle
        TR%R(i)%PositionAngle=TR%R(i)%PositionAngle-Pi/2.
c
c       Next do all the velocities
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VSys)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VRot)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VRad)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%Vvert)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VDisp)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%dvdz)
c       Do the surface brightness conversion
        call GeneralBrightnessConversion(TR%UnitSwitchs(3)
     &          ,TR%R(i)%Sigma,ChannelSize,BeamArea
     &          ,DC%DH%PixelSize(0))
c       Finally convert the radii to pixels
        call Length_PixelConversion(TR%R(i)%Rmid
     &              ,TR%R(i)%Rmid,DC%DH%PixelSize(0))
        call Length_PixelConversion(TR%R(i)%Rwidth
     &              ,TR%R(i)%Rwidth,DC%DH%PixelSize(0))
c       Also make sure to convert z0 and z0Grad to pixels
        call Length_PixelConversion(TR%R(i)%z0
     &              ,TR%R(i)%z0,DC%DH%PixelSize(0))
        call Length_PixelConversion(TR%R(i)%zGradiantStart
     &              ,TR%R(i)%zGradiantStart,DC%DH%PixelSize(0))
      enddo


      return
      end subroutine
ccccccccc



ccccc
c       This routine does the angular conversions for the data cube parameters
      subroutine TextDCHeader_AngularConversions(DC,Beam)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(INOUT) :: Beam

c       Do the data cube unit conversions
c           Start with the pixel sizes
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(0)
     &          ,DC%DH%PixelSize(0))    !This File
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(0)
     &          ,DC%DH%PixelSize(1))    !This File
c           Next the channel size
      call GeneralVelocityConversion(DC%DH%DimensionUnitSwitch(2)
     &          ,DC%DH%ChannelSize)     !This File
c       Now do up the reference values
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(1)
     &          ,DC%DH%RefVal(0))       !This File
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(1)
     &          ,DC%DH%RefVal(1))       !This File
      call GeneralVelocityConversion(DC%DH%DimensionUnitSwitch(3)
     &          ,DC%DH%RefVal(2))       !This File

c      Next do the Beam unit conversions
c           Want to have the beam in pixel dimensions
      call Length_PixelConversion(Beam%BeamMajorAxis
     &              ,Beam%BeamMajorAxis,DC%DH%PixelSize(0))
      call Length_PixelConversion(Beam%BeamMinorAxis
     &              ,Beam%BeamMinorAxis,DC%DH%PixelSize(0))
      call GeneralAngularConversion(Beam%BeamUnitsSwitch(1)
     &          ,Beam%BeamPositionAngle)    !This File

      return
      end subroutine
ccccccc



cccccccc
c
      subroutine DataCubeUnitConversions(DC,Beam)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(INOUT) :: Beam
      real z
      real TempFreq,V1,V2
      integer i

c       Do the velocity conversion first --- goal is to have the units in km/s
      call ChannelVelocityConversion(DC)

c       Do the angule conversion next -- initially want units of arcseconds
      if(trim(DC%DH%Units(0)) .eq. 'DEGREE'
     &              .or. trim(DC%DH%Units(0)) .eq. 'DEGREES'
     &              .or. trim(DC%DH%Units(0)) .eq. 'deg') then
        call DegreesToArcSec(Beam%BeamMajorAxis,Beam%BeamMajorAxis)
        call DegreesToArcSec(Beam%BeamMinorAxis,Beam%BeamMinorAxis)
        do i=0,1
            call DegreesToArcSec(DC%DH%PixelSize(i),DC%DH%PixelSize(i))
            call DegreesToArcSec(DC%DH%RefVal(i),DC%DH%RefVal(i))
        enddo
        Beam%PixelSize=DC%DH%PixelSize
      endif
c       Finally, put the beam units into pixels as the TR portion of the
c           code will work in pixel space
      call Length_PixelConversion(Beam%BeamMajorAxis
     &                  ,Beam%BeamMajorAxis,DC%DH%PixelSize(0))
      call Length_PixelConversion(Beam%BeamMinorAxis
     &                  ,Beam%BeamMinorAxis,DC%DH%PixelSize(0))

      return
      end subroutine
cccccccc



ccccccc
c     This routine does brightness conversion for a data cube to
c       get final units of Jy/pixel
      subroutine DCBrightnessConversion(DC,Beam)
      implicit none
      Type(Beam2D),INTENT(IN) :: Beam
      Type(DataCube),INTENT(INOUT) :: DC

      real BeamArea

c       The cell brighness needs to be in Jy/pixel
c      print*, "DataCube Units ", DC%DH%FUnit
      if(DC%DH%FUnit .eq. 'Jy/beam'
     &       .or. DC%DH%FUnit .eq. 'Jy/Beam') then
c              The beam area in pixels should already be
c                   calculated in the beam allocation.
        BeamArea=Beam%BeamAreaPixels
        print*, "Input beam area", BeamArea,sum(DC%Flux)
     &              ,sum(DC%Flux/BeamArea)
        DC%Flux=DC%Flux/BeamArea

      endif

      return
      end subroutine
ccccccc


cccccc
      subroutine ChannelVelocityConversion(DC)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC

      logical AcceptableType, AcceptableUnits
      character(8) PossibleChannelKeys(8)
      character(8) PossibleVelKeys(2)
      character(8) PossibleFreqKeys(2)
      integer i,j,k

      real FreqTemp(2),VelTemp(2),Redshifts(2)

      AcceptableType=.False.
      AcceptableUnits=.False.
      PossibleChannelKeys(1)='VELO-LSR'
      PossibleChannelKeys(2)='VELOHEL '
      PossibleChannelKeys(3)='VOPT'
      PossibleChannelKeys(4)='VELO'
      PossibleChannelKeys(5)='FREQ'

      PossibleVelKeys(1)="m/s"
      PossibleVelKeys(2)="km/s"

      PossibleFreqKeys(1)="Hz"


c       First loop through all possible keys to see if the axis type is allowed
      do i=1,4
        if(trim(DC%DH%AxisType(2)) .eq.
     &          trim(PossibleChannelKeys(i))) then
            AcceptableType=.True.
c           If it is acceptable, check on the units
c               1 or 2 is velocity
            if( i.le. 4) then
                do j=1,2
                    if(trim(DC%DH%Units(2))
     &                  .eq. trim(PossibleVelKeys(j))) then
                        AcceptableUnits=.True.
c                           For m/s, convert all units to km/s
                        if(j .eq. 1) then
                            call MToKM(DC%DH%ChannelSize
     &                              ,DC%DH%ChannelSize)
                            call MToKM(DC%DH%ChannelCent
     &                              ,DC%DH%ChannelCent)
                            call MToKM(DC%DH%RefVal(2)
     &                              ,DC%DH%RefVal(2))
                        endif
                    endif
                enddo
            elseif(i .eq. 5) then       ! 4 is frequency
c               Currently frequency is not supported in the code
                AcceptableType=.False.
                do j=1,1
                    if(trim(DC%DH%Units(2))
     &                  .eq. trim(PossibleFreqKeys(j))) then
                        AcceptableUnits=.True.
                    endif
                enddo
            endif
        endif
      enddo


      if(AcceptableType .eqv. .False.) then
c           The code only
        if(DC%DH%MaskSwitch .eqv. .False.) then
            print*, "Channel type is not supported for "
     &              , "kinematic models "
     &              , trim(DC%DH%AxisType(2))
            print*, "Acceptable types are: "
            do i=1,4
                print*,trim(PossibleChannelKeys(i))
            enddo
            stop
        endif
      endif

      if(AcceptableUnits .eqv. .False.) then
        if(DC%DH%MaskSwitch .eqv. .False.) then
            print*, "Channel units are not supported for "
     &              , "kinematic models "
     &              , trim(DC%DH%Units(2))
            print*, "Acceptable units are: "
            do i=1,2
                print*, trim(PossibleVelKeys(i))
            enddo
            stop
        endif
      endif
 


      return
      end subroutine
cccccc





      end module
