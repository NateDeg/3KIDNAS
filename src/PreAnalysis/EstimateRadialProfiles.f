cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       estimate the profiles across the major axis of
c       the galaxy
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module EstimateProfileMod
      use CommonConsts
      use DataCubeMod
      use BeamMod
      use TiltedRingMod
      use InterpolateMod

      implicit none
      contains


ccccc
      subroutine EstimateProfiles(Maps,Beam,Center
     &              ,incl,PA,Size
     &              ,EstimatedProfile,RadialProfile
     &              ,VSys,FittingOptions,nRings
     &              ,SDLims,VLims,NoiseFrac)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      Type(Beam2D),INTENT(IN) :: BEAM
      real,intent(IN) :: incl,Size,VSys
      real,INTENT(IN) :: Center(0:1)
      real,intent(INOUT) :: PA
      Type(TiltedRingFittingOptions) FittingOptions
      integer,intent(INOUT) :: nRings
      real,INTENT(IN) :: SDLims(2),VLims(2)
      real,INTENT(IN) :: NoiseFrac



      real XTest,PosRot(0:1),Area,XTest2
      integer i,j,k
      integer nRingsMax

      real,ALLOCATABLE,INTENT(INOUT) :: EstimatedProfile(:,:)
      real,ALLOCATABLE,INTENT(INOUT) :: RadialProfile(:,:)
      real,ALLOCATABLE :: TempRadialProfile(:,:)
      real,ALLOCATABLE :: ProjectedProfile(:,:)
      real Avg(2),t1,t2

      integer,ALLOCATABLE :: ModelerableRingSwitch(:)
      integer,ALLOCATABLE :: RadialProfileModelSwitch(:)

      real RTest
      real Noise_SDLim

      real LeadingVelProfile


      real AvgSD, AvgRC, AvgR
      real sumX,sumX2,sumY,sumXY,a,b,nUse
      integer LastPoint

      real RCorr,dR
      integer RIndx
      real SDL(2),SDH(2),RCL(2),RCH(2)
      real,ALLOCATABLE :: SDCorr(:),RCCorr(:)

cccc

c
c       First set the maximum number of rings possible based on the map shape
c

      RTest=Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.
      RTest=sqrt(RTest)/2.
      RTest=RTest/Beam%BeamMajorAxis
      nRings=int(RTest*FittingOptions%nRingsPerBeam)+1
c
      nRingsMax=nRings
      print*, "Map shape", Maps%DH%nPixels(0),Maps%DH%nPixels(1),RTest

      print*, "Mom0 flux", sum(Maps%Flux(:,:,0))

      print*, "Estimating shape check", PA*180./Pi, incl*180./Pi
     &                  ,VSys, Center

c       Allocate the full profile --- due to the inclusion of zero, there is a plus one in the allocation
      ALLOCATE(EstimatedProfile(0:2,-nRings+1:nRings))
c       Allocate the radial profile arrays
      ALLOCATE(TempRadialProfile(0:2,0:nRings-1))
      ALLOCATE(ProjectedProfile(0:2,0:nRings-1))
c       Allocate the number ring switch array
      ALLOCATE(ModelerableRingSwitch(0:nRings-1))


c       Estimate the profiles by slicing across the moment maps
      call SliceMap(Maps,Beam,Center
     &              ,incl,PA,Size, nRings
     &              ,EstimatedProfile
     &              ,FittingOptions%nRingsPerBeam)
c       Now build the initial profiles by averaging the sliced profile
        call BuildIniProfile(nRings, VSys,Incl
     &              ,ModelerableRingSwitch
     &              ,EstimatedProfile
     &              ,TempRadialProfile
     &              ,ProjectedProfile)
c       Next go through the initial profiles and decide whether each ring can be modelled
c           Start by setting the noise limit
      Noise_SDLim=NoiseFrac*Maps%DH%Uncertainty*abs(Maps%DH%ChannelSize)
c       Now go through everything on a point by point basis
      do i=1, nRings
        call CheckProfilePointModelability(SDLims
     &              ,VLims,Noise_SDLim
     &              ,TempRadialProfile(0:2,i-1)
     &              ,ProjectedProfile(0:2,i-1)
     &              ,ModelerableRingSwitch(i-1))
      enddo
c       Now it's necessary to figure out the number of rings to model

      print*, "Total potentially modelable rings", nRings
c       When not using an input profile, use the modelable switch
      if(FittingOptions%nTargRings .eq. -1) then
        nRings=sum(ModelerableRingSwitch)
c       Otherwise use the number of input rings
      else
        nRings=FittingOptions%nTargRings
      endif
      print*, "Total modelable number of rings", nRings,nRingsMax
      if(nRings .eq. 0) then
        print*, "No modelable rings have been found"
        stop
      endif
c       Allocate the radial profiles
      ALLOCATE(RadialProfile(0:2,0:nRings-1))
      ALLOCATE(RadialProfileModelSwitch(0:nRings-1))
c       Assign the 'successful' radial profile estimates to the profile array
c           If the code is automatically assigning the rings, just go through the main loop and use the modelerable ring switch
      call AssignProfile(nRings,nRingsMax
     &              ,FittingOptions
     &              ,RadialProfile,TempRadialProfile
     &              ,ModelerableRingSwitch
     &              ,RadialProfileModelSwitch,Maps)
c       When using a pre-defined grid, it's possible that a ring might be empty that needs an initial estimate as it may have been selected as not modelable

c       Now fill in all the ring values
      call FillInMisingRings(nRings
     &          ,RadialProfile
     &          ,RadialProfileModelSwitch)

c       The next step is to do beam corrections
      call BeamCorrectProfile(nRings
     &      ,RadialProfile
     &      ,Beam
     &      ,SDLims,VLims)


c       Finally do a check on the profile to make sure that we have the PA in the right direction
      LeadingVelProfile=sum(EstimatedProfile(2,1:nRings)-VSys)

      print*, "PA Ini",PA
      if(LeadingVelProfile .le. 0.) then
        PA=PA+Pi
        print*, "PA Adjust 1",PA
        if(PA .ge. 2.*pi) PA=PA-2.*pi    !Keep the PA between 0<PA<Phi
      endif
      print*, "Final PA Check", PA


c       And deallocate the temporary profiles
      DEALLOCATE(TempRadialProfile)
      DEALLOCATE(ModelerableRingSwitch)
      DEALLOCATE(RadialProfileModelSwitch)

      return
      end subroutine
cccccc

cccccc
      subroutine SliceMap(Maps,Beam,Center
     &              ,incl,Phi,Size, nRings
     &              ,AvgQuantities,nRingsPerBeam)

      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      Type(Beam2D),INTENT(IN) :: BEAM
      real,intent(IN) :: incl,Size
      real,INTENT(IN) :: Center(0:1)
      real,intent(INOUT) :: Phi
      integer,intent(INOUT) :: nRings,nRingsPerBeam
      real,intent(INOUT) :: AvgQuantities(0:2,-nRings+1:nRings)

      real BeamPix,CentPix(0:1)
      real SliceLimit
      integer i,j,k,l
      real x,y,xPrime,yPrime

      real AvgX(-nRings:nRings)
      integer nBin(-nRings:nRings)

c      print*, "Slicing the map to get profile estimates"
c       Set the limit on the slice thickness in beams (right now testing with 0.5 beams)
      SliceLimit=0.5
c       Get the thickness of the beam in pixels
      BeamPix=abs(Beam%BeamMajorAxis)!/(Maps%DH%PixelSize(0)))
c       Get the center position in pixel space
      CentPix=(Center(0:1))!-Maps%DH%Start(0:1))
c     &          /(Maps%DH%PixelSize(0:1))
c      print*, "Center in pixels", CentPix,Phi

c       Initialize all the arrays to zero
      do i=-nRings+1,nRings
        AvgQuantities(0:2,i)=0.
        nBin(i)=0
        AvgX(i)=0.
      enddo



c       Now loop through all pixels in the map!
      do  i=0, Maps%DH%nPixels(0)-1
        do j=0, Maps%DH%nPixels(1)-1
c               Get the pixel positions relative to the center pixel
            x=real(i-CentPix(0))
            y=real(j-CentPix(1))
c               Rotate the positions by the position angle
            xPrime=x*cos(-Phi)-y*sin(-Phi)
            yPrime=x*sin(-Phi)+y*cos(-Phi)
c               Convert the x value to units of rings
            xPrime=xPrime/(BeamPix/nRingsPerBeam)
c               Convert the y value to units of beams
            yPrime=yPrime/BeamPix
c           Round these points to the nearest ring
            k=nint(xPrime+0.5)
            l=nint(yPrime)
c           Only use those points within the target thickness in terms of beams
            if(abs(yPrime) .le. SliceLimit) then

c           Check if within nRings of the center
c                if(abs(k) .le. nRings) then
                if(k .ge. -nRings+1 .and. k .le. nRings) then
c               Check that the pixel has flux
                    if(Maps%Flux(i,j,0) .ge. 0.
     &                  .and. isnan(Maps%Flux(i,j,1))
     &                  .eqv. .False.) then
c               Add in the moment map values for this 'beam' averaged quantity
                        AvgQuantities(1,k)=AvgQuantities(1,k)
     &                      +Maps%Flux(i,j,0)
                        AvgQuantities(2,k)=AvgQuantities(2,k)
     &                      +Maps%Flux(i,j,1)
                        nBin(k)=nBin(k)+1
                        AvgX(k)=AvgX(k)+xPrime

                    endif

                endif
            endif
        enddo
      enddo

c       Now loop through all the beam bins and get the average values
      do i=-nRings+1,nRings
c        print*, "pre average", i,AvgQuantities(0:2,i),nBin(i)
        AvgQuantities(0,i)=Beam%BeamMajorAxis*(real(i)-0.5)
     &                      /nRingsPerBeam
        AvgQuantities(1,i)=AvgQuantities(1,i)/real(nBin(i))
        AvgQuantities(1,i)=AvgQuantities(1,i)
     &          /abs(Maps%DH%ChannelSize)
        AvgQuantities(2,i)=AvgQuantities(2,i)/real(nBin(i))
        AvgX(i)=AvgX(i)/real(nBin(i))
        print*, "Averaged rings", i, AvgQuantities(0:2,i),nBin(i)
     &          ,AvgX(i)
     &          ,AvgQuantities(0,i)*abs(Maps%DH%PixelSize(0))
      enddo


      return
      end subroutine
ccccccc


cccc
      subroutine AssignProfile(nRings,nRingsMax
     &              ,FittingOptions
     &              ,RadialProfile,TempRadialProfile
     &              ,ModelerableRingSwitch
     &              ,RadialProfileModelSwitch
     &              ,Maps)
      use PipelineGlobals
      implicit none
      integer,INTENT(IN) :: nRings, nRingsMax
      Type(TiltedRingFittingOptions),INTENT(IN):: FittingOptions
      real,INTENT(IN) :: TempRadialProfile(0:2,0:nRingsMax-1)
      real,INTENT(INOUT) :: RadialProfile(0:2,0:nRings-1)
      integer,INTENT(IN) :: ModelerableRingSwitch(0:nRingsMax-1)
      integer, INTENT(INOUT) :: RadialProfileModelSwitch(0:nRings-1)
      Type(DataCube),INTENT(IN) :: Maps

      integer i,j
      real PixSize


c       Assign the 'successful' radial profile estimates to the profile array
c           If using a supplied grid, adjust for arcsec's to pixels
      if(FittingOptions%nTargRings .ge. 1) then
        PixSize=abs(Maps%DH%PixelSize(0))
        RadGrid=RadGrid/PixSize
      endif
c           If the code is automatically assigning the rings, just go through the main loop and use the modelerable ring switch
      j=0
      do i=0,nRingsMax-1
        if(FittingOptions%nTargRings .eq. -1) then
c           Without a grid, check on the whether the ring is modelable
            if(ModelerableRingSwitch(i) .eq. 1) then
                RadialProfile(0:2,j)=TempRadialProfile(0:2,i)
                RadialProfileModelSwitch(j)=ModelerableRingSwitch(i)
                j=j+1
            endif
        else
c           With a grid, check that the radius matches
            if( abs(TempRadialProfile(0,i)
     &             - RadGrid(j+1)) .lt. 0.01) then
                RadialProfile(0:2,j)=TempRadialProfile(0:2,i)
                RadialProfileModelSwitch(j)=ModelerableRingSwitch(i)
                j=j+1
            endif
        endif

        if(j .gt. 0) then
            print*, "Profile assigned", i,j-1
     &                  ,RadialProfile(0:2,j-1)
     &                  ,TempRadialProfile(0:2,i)
     &                  ,ModelerableRingSwitch(i)
     &                  ,RadialProfileModelSwitch(j-1)
        endif
c        When at the end of assigning the profile, end the loop
        if(j .eq. nRings) return
      enddo



      return
      end subroutine
ccccc


ccccc
c
      subroutine BuildIniProfile(nRings, VSys,Incl
     &              ,ModelerableRingSwitch
     &              ,EstimatedProfile
     &              ,RadialProfile
     &              ,ProjectedProfile)
      implicit none
      integer, INTENT(IN) :: nRings
      integer, INTENT(INOUT) :: ModelerableRingSwitch(0:nRings-1)
      real, INTENT(IN) :: VSys,Incl
      real, INTENT(IN) :: EstimatedProfile(0:2,-nRings+1:nRings)
      real, INTENT(INOUT) :: RadialProfile(0:2,0:nRings-1)
      real, INTENT(INOUT) :: ProjectedProfile(0:2,0:nRings-1)

      real InclU

      integer i, j
      real t1,t2

      if(Incl .gt. 70./180.*Pi) then
        InclU=70./180.*Pi
      else
        InclU=Incl
      endif

      do i=1,nRings
c           Assume that the ring is fine
        ModelerableRingSwitch(i-1)=1

        j=-i+1
c           Get the absolute difference from Vsys
          t1=abs(EstimatedProfile(2,i)-VSys)
          t2=abs(EstimatedProfile(2,j)-VSys)
c           It's possible that one of the bins might end up empty due to holes and other issues
c               In that case set the velocity difference to the paired value
        if(EstimatedProfile(2,i) .le. 0) t1=t2
        if(EstimatedProfile(2,j) .le. 0) t2=t1

c       Now average the velocity terms together
        ProjectedProfile(2,i-1)=(t1+t2)/2.
c           Next average the surface density terms together
        ProjectedProfile(1,i-1)=(EstimatedProfile(1,i)
     &              +EstimatedProfile(1,j))/2.
c           Place the averages into the radial profile structure
        RadialProfile(0,i-1)=abs(EstimatedProfile(0,j))
c       For the surface density, normalize by the cosine
        RadialProfile(1,i-1)=ProjectedProfile(1,i-1)*cos(InclU)
c           For the velocity, normalize by the inclination to get a better guess
        RadialProfile(2,i-1)=ProjectedProfile(2,i-1)/sin(Incl)
      enddo

      return
      end subroutine
ccccccc


cccccc
      subroutine CheckProfilePointModelability(SDLims
     &              ,VLims,Noise_SDLim
     &              ,ProfilePt,ProjectedProfilePt
     &              ,ModelSwitch)

      implicit none
      real,INTENT(IN):: SDLims(2), VLims(2),Noise_SDLim
      real,INTENT(INOUT):: ProfilePt(0:2),ProjectedProfilePt(0:2)
      integer,INTENT(INOUT):: ModelSwitch


c               A SET OF PRINT OUTS FOR CHECKING THINGS
      print*,"Estimated Radial Profile", ProfilePt(0:2)
     *          ,ProjectedProfilePt(1),Noise_SDLim

c       Do a check to see if the are problems with the ring
c           First check on negative SDs
      if (ProfilePt(1).le. 0.) then
c            print*, "negative SD value", i, TempRadialProfile(0:2,i-1)
            ModelSwitch = 0
      endif
c           Do a second check for NaN's
      if (ProfilePt(1).ne. ProfilePt(1)) then
c            print*, "NaN SD value", i, TempRadialProfile(0:2,i-1)
            ModelSwitch = 0
      endif
c           Next check that the SD value isn't below the limit implied by the noise level -- note that here we use the non-inclination corrected profile
      if (ProjectedProfilePt(1) .le. Noise_SDLim) then
c            print*, "SD initial estimate below noise"
c     &          ,i, TempRadialProfile(1,i-1),Noise_SDLim,Avg(1)
            ModelSwitch = 0
      endif

c           Do a third check for low and high SD values (and replace them with the limits)
c      if (ProfilePt(1) .le. SDLims(1)) then
cc      if (ProjectedProfilePt(1) .le. SDLims(1)) then
c        print*, "Low SD initial estimate"
c     &                  ,ProfilePt(0:2),SDLims(1)
c        ProfilePt(1)=SDLims(1)
c        ModelSwitch = 0      !TEMPORARY FOR TESTING PURPOSES
c      endif
      if (ProfilePt(1) .ge. SDLims(2)) then
        print*, "High SD initial estimate"
        ProfilePt(1)=SDLims(2)
      endif
c           Do the same check for the rotation velocity
      if (ProfilePt(2) .le. VLims(1)) then
        print*, "Low rotation velocity initial estimate"
        ProfilePt(2)=VLims(1)
      endif
      if (ProfilePt(2) .ge. VLims(2)) then
        print*, "High VRot initial estimate"
        ProfilePt(2)=VLims(2)
      endif
      return
      end subroutine
ccccc

cccc
c
      subroutine FillInMisingRings(nRings
     &          , RadialProfile
     &          , ModelerableRingSwitch)
      implicit none
      integer, INTENT(IN) :: nRings
      real,INTENT(INOUT) :: RadialProfile(0:2,0:nRings-1)
      integer,INTENT(IN) :: ModelerableRingSwitch(0:nRings-1)

      integer i,j,k

      do i=0, nRings-1
        print*, "Radial profiles before",i,RadialProfile(0:2,i)
     &              ,ModelerableRingSwitch(i)
c           Check if we have a value for the ring -- if not, then we'll fill it in
        if(ModelerableRingSwitch(i) .eq. 0) then
c           Now check if we are the first ring
            if(i .eq. 0) then
c           For this ring loop to the right
                do k=i, nRings-1
                    if(ModelerableRingSwitch(k) .eq. 1) then
                        print*, "Fill in ring by", i,k
     &                      ,RadialProfile(0,i)
     &                      ,RadialProfile(1:2,k)
                        RadialProfile(1:2,i)=RadialProfile(1:2,k)
                        exit
                    endif
                enddo
c           If there is no ring with a viable value is found, then there is no point in running the fit, so we'll stop here
                print*, "No ring with an acceptable"
     &                  //" surface density found"
                stop
c           After the first ring is done, there must be at least one value to the left of the empty ring.  We'll use that value for the current ring
        else
                do k=i,0,-1
                    if(ModelerableRingSwitch(k) .eq. 1) then
                        RadialProfile(1:2,i)=RadialProfile(1:2,k)
                        exit
                    endif
                enddo
            endif
        endif
c       Add a print out here to check if the profiles make sense
      print*, "Radial profiles after",i,RadialProfile(0:2,i)
      enddo


      return
      end subroutine
cccccc

ccccccc
      subroutine BeamCorrectProfile(nRings
     &      ,RadialProfile
     &      ,Beam
     &      ,SDLims,VLims)
      implicit none
      integer, INTENT(IN) :: nRings
      real,INTENT(INOUT) :: RadialProfile(0:2,0:nRings-1)
      Type(Beam2D),INTENT(IN) :: BEAM
      real,INTENT(IN):: SDLims(2), VLims(2)

      real,ALLOCATABLE :: SDCorr(:),RCCorr(:)
      real dR,RCorr
      integer i,k,RIndx
      real SDL(2),SDH(2),RCL(2),RCH(2)

      print*, "Beam Corr", shape(RadialProfile)
c
c       Only do corrections if there are more than one ring
      if (nRings .lt. 2) return
c
      ALLOCATE(SDCorr(0:nRings-1))
      ALLOCATE(RCCorr(0:nRings-1))
c       Now try to correct the radial profile for the beam
      dR=RadialProfile(0,1)-RadialProfile(0,0)
c       It is necessary to start RIndx at 0 and then add from there
      do i=0,nRings-1
        print*, i, RadialProfile(0,i), RadialProfile(1,i)
     &              , RadialProfile(2,i)
        k=2
c           Get the corrected radius....note that this may need some adjustment in the innermost region
100     RCorr=RadialProfile(0,i)**2.-(Beam%BeamMajorAxis/k)**2.
        if (RCorr .lt. 0.) then
            k=k+1
            goto 100
        endif
c        print*, RCorr,sqrt(RCorr),Beam%BeamMajorAxis,dR
        RCorr=sqrt(RCorr)
c           Get the index for the corrected radius
c               Do this by slowly increasing RIndx
 200    continue
        if (RCorr .gt. RadialProfile(0,RIndx)) then
            RIndx=RIndx+1
            goto 200
        endif

c        RIndx=(RCorr-RadialProfile(0,0))/dR
        if(RIndx .eq. nRings-1) then
            RIndx=RIndx-1
        endif
        SDL(1)=RadialProfile(0,RIndx)
        SDL(2)=RadialProfile(1,RIndx)
        SDH(1)=RadialProfile(0,RIndx+1)
        SDH(2)=RadialProfile(1,RIndx+1)
        call SimpleInterpolateY(SDL,SDH,RCorr,SDCorr(i))

        RCL(1)=RadialProfile(0,RIndx)
        RCL(2)=RadialProfile(2,RIndx)
        RCH(1)=RadialProfile(0,RIndx+1)
        RCH(2)=RadialProfile(2,RIndx+1)
        call SimpleInterpolateY(RCL,RCH,RCorr,RCCorr(i))

        if(RCCorr(i) .gt. VLims(2)) then
            RCCorr(i)=VLims(2)
        elseif(RCCorr(i) .lt. VLims(1)) then
            RCCorr(i)=VLims(1)
        endif

        if(SDCorr(i) .gt. SDLims(2)) then
            SDCorr(i)=SDLims(2)
        elseif(SDCorr(i) .lt. SDLims(1)) then
            SDCorr(i)=SDLims(1)
        endif

c        print*, i, RadialProfile(0,i),RCorr,RadialProfile(1,i),SDCorr(i)
c     &          ,RadialProfile(2,i),RCCorr(i)

      enddo
      RadialProfile(1,0:nRings-1)=SDCorr(0:nRings-1)
      RadialProfile(2,0:nRings-1)=RCCorr(0:nRings-1)
      DEALLOCATE(SDCorr)
      DEALLOCATE(RCCorr)

      return
      end subroutine
ccccc


      end module
