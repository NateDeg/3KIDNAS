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
     &              ,SDLims,VLims)
      implicit none
      Type(DataCube),INTENT(IN) :: Maps
      Type(Beam2D),INTENT(IN) :: BEAM
      real,intent(IN) :: incl,Size,VSys
      real,INTENT(IN) :: Center(0:1)
      real,intent(INOUT) :: PA
      Type(TiltedRingFittingOptions) FittingOptions
      integer,intent(INOUT) :: nRings
      real,INTENT(IN) :: SDLims(2),VLims(2)



      real XTest,PosRot(0:1),Area,XTest2
      integer i,j,k
      integer nRingsMax

      real,ALLOCATABLE,INTENT(INOUT) :: EstimatedProfile(:,:)
      real,ALLOCATABLE,INTENT(INOUT) :: RadialProfile(:,:)
      real,ALLOCATABLE :: TempRadialProfile(:,:)
      real Avg(2),t1,t2

      integer,ALLOCATABLE :: ModelerableRingSwitch(:)

      real RTest
      real Noise_SDLim,NoiseFrac

      real LeadingVelProfile


      real AvgSD, AvgRC, AvgR
      real sumX,sumX2,sumY,sumXY,a,b,nUse
      integer LastPoint

      real RCorr,dR
      integer RIndx
      real SDL(2),SDH(2)
      real,ALLOCATABLE :: SDCorr(:)

cccc

c
c       First set the maximum number of rings possible based on the map shape
c

      RTest=Maps%DH%nPixels(0)**2.+Maps%DH%nPixels(1)**2.
      RTest=sqrt(RTest)/2.
      RTest=RTest/Beam%BeamMajorAxis
      nRings=int(RTest*FittingOptions%nRingsPerBeam)+1

      if(FittingOptions%nTargRings .eq. -1) then
        nRings=int(RTest*FittingOptions%nRingsPerBeam)+1
      else
        nRings=FittingOptions%nTargRings
      endif

c      nRings=20
      nRingsMax=nRings
      print*, "Map shape", Maps%DH%nPixels(0),Maps%DH%nPixels(1),RTest

      print*, "Mom0 flux", sum(Maps%Flux(:,:,0))

      print*, "Estimating shape check", PA*180./Pi, incl*180./Pi
     &                  ,VSys, Center

c       Set the local nRings to that given in the fitting options
c      nRings=FittingOptions%nRings
      print*, 'nBeam', nRings
c       Allocate the full profile --- due to the inclusion of zero, there is a plus one in the allocation
      ALLOCATE(EstimatedProfile(0:2,-nRings+1:nRings))
c       Allocate the radial profile arrays
      ALLOCATE(TempRadialProfile(0:2,0:nRings-1))
c       Allocate the number ring switch array
      ALLOCATE(ModelerableRingSwitch(0:nRings-1))


c       Estimate the profiles by slicing across the moment maps
      call SliceMap(Maps,Beam,Center
     &              ,incl,PA,Size, nRings
     &              ,EstimatedProfile
     &              ,FittingOptions%nRingsPerBeam)

c       By definition the profile should go from low to high velocity.  If
c           that's not true, adjust the position angle



c       Now average everything the profile together for all rings
c
      NoiseFrac=1.
      Noise_SDLim=NoiseFrac*Maps%DH%Uncertainty*abs(Maps%DH%ChannelSize)
      


      print*, "Channelsize",abs(Maps%DH%ChannelSize)
     &          ,NoiseFrac,Maps%DH%Uncertainty
     &          ,Noise_SDLim, SDLims
      do i=1,nRings
c           Assume that the ring is fine
        ModelerableRingSwitch(i-1)=1

        j=-i+1
c           Get the absolute difference from Vsys
        t1=abs(EstimatedProfile(2,i)-VSys)
        t2=abs(EstimatedProfile(2,j)-VSys)
c        print*, "Rot check", t1, t2, VSys
c           It's possible that one of the bins might end up empty due to holes and other issues
c               In that case set the velocity difference to the paired value
        if(EstimatedProfile(2,i) .le. 0) t1=t2
        if(EstimatedProfile(2,j) .le. 0) t2=t1

c       Now average the velocity terms together
        Avg(2)=(t1+t2)/2.
c           Next average the surface density terms together
        Avg(1)=(EstimatedProfile(1,i)+EstimatedProfile(1,j))/2.
c           Place the averages into the radial profile structure
        TempRadialProfile(0,i-1)=abs(EstimatedProfile(0,j))
c       For the surface density, normalize by the cosine
        TempRadialProfile(1,i-1)=Avg(1)*cos(Incl)
c           For the velocity, normalize by the inclination to get a better guess
        TempRadialProfile(2,i-1)=Avg(2)/sin(Incl)


        print*, "SD consistency checks", i, Avg(1)
     &              ,EstimatedProfile(1,i),EstimatedProfile(1,j)
        print*,"Estimated Radial Profile", i,TempRadialProfile(0:2,i-1)
     &          ,VSys
     *          ,Noise_SDLim
c       Do a check to see if the are problems with the ring
c           First check on negative SDs
        if (TempRadialProfile(1,i-1).le. 0.) then
c            print*, "negative SD value", i, TempRadialProfile(0:2,i-1)
            ModelerableRingSwitch(i-1) = 0
        endif
c           Do a second check for NaN's
        if (TempRadialProfile(1,i-1).ne. TempRadialProfile(1,i-1)) then
c            print*, "NaN SD value", i, TempRadialProfile(0:2,i-1)
            ModelerableRingSwitch(i-1) = 0
        endif
c           Next check that the SD value isn't below the limit implied by the noise level
        if (TempRadialProfile(1,i-1) .le. Noise_SDLim) then
            print*, "SD initial estimate below noise"
     &          ,i, TempRadialProfile(1,i-1),Noise_SDLim
            ModelerableRingSwitch(i-1) = 0
        endif

c           Do a third check for low and high SD values (and replace them with the limits)
        if (TempRadialProfile(1,i-1) .le. SDLims(1)) then
            print*, "Low SD initial estimate"
     &                  ,TempRadialProfile(1,i-1),SDLims(1)
            TempRadialProfile(1,i-1)=SDLims(1)
            ModelerableRingSwitch(i-1) = 0      !TEMPORARY FOR TESTING PURPOSES
        endif
        if (TempRadialProfile(1,i-1) .ge. SDLims(2)) then
            print*, "High SD initial estimate"
            TempRadialProfile(1,i-1)=SDLims(2)
        endif
c           Do the same check for the rotation velocity
        if (TempRadialProfile(2,i-1) .le. VLims(1)) then
            print*, "Low rotation velocity initial estimate"
            TempRadialProfile(2,i-1)=VLims(1)
        endif
        if (TempRadialProfile(2,i-1) .ge. VLims(2)) then
            print*, "High VRot initial estimate"
            TempRadialProfile(2,i-1)=VLims(2)
        endif

      enddo


      print*, "Total modelerable rings", nRings
      if(FittingOptions%nTargRings .eq. -1) then
        nRings=sum(ModelerableRingSwitch)
      else
        nRings=FittingOptions%nTargRings
      endif
      print*, "Recalculated number of rings", nRings,nRingsMax
      ALLOCATE(RadialProfile(0:2,0:nRings-1))
c       Assign the 'successful' radial profile estimates to the profile array
      j=0
      do i=0,nRingsMax-1
        if(ModelerableRingSwitch(i) .eq. 1) then
            RadialProfile(0:2,j)=TempRadialProfile(0:2,i)
            j=j+1
        endif
      enddo

c       Now fill in all the ring values

      do i=0, nRings-1
        print*, "Radial profiles before",i,RadialProfile(0:2,i)
c           Check if we have a value for the ring -- if not, then we'll fill it in
        if(ModelerableRingSwitch(i) .eq. 0) then
c           Now check if we are the first ring
            if(i .eq. 0) then
c           For this ring loop to the right
                do k=i, nRings-1
                    if(ModelerableRingSwitch(k) .eq. 1) then
                        RadialProfile(1:2,i)=TempRadialProfile(1:2,k)
                        RadialProfile(0,i)=TempRadialProfile(0,i)
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
                        RadialProfile(1:2,i)=TempRadialProfile(1:2,k)
                        RadialProfile(0,i)=TempRadialProfile(0,i)
                        exit
                    endif
                enddo
            endif
        endif
c       Add a print out here to check if the profiles make sense
        print*, "Radial profiles after",i,RadialProfile(0:2,i)
      enddo


c       Only do corrections if there are more than one ring
      if (nRings .ge. 2) then
        ALLOCATE(SDCorr(0:nRings-1))
c       Now try to correct the radial profile for the beam
        dR=RadialProfile(0,1)-RadialProfile(0,0)
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

        RIndx=(RCorr-RadialProfile(0,0))/dR
        if(RIndx .eq. nRings-1) then
            RIndx=RIndx-1
        endif
        SDL(1)=RadialProfile(0,RIndx)
        SDL(2)=RadialProfile(1,RIndx)
        SDH(1)=RadialProfile(0,RIndx+1)
        SDH(2)=RadialProfile(1,RIndx+1)
        call SimpleInterpolateY(SDL,SDH,RCorr,SDCorr(i))
        print*, i, RadialProfile(0,i),RCorr,RadialProfile(1,i),SDCorr(i)

      enddo
        RadialProfile(1,0:nRings-1)=SDCorr(0:nRings-1)
        DEALLOCATE(SDCorr)
      endif


c
c      print*, "Position angle checkish",EstimatedProfile(2,nRings-1)
c      print*, "sum again",sum(EstimatedProfile(2,-nRings+1:0)-VSys)
      LeadingVelProfile=sum(EstimatedProfile(2,1:nRings)-VSys)



      print*, "PA Ini",PA
      if(LeadingVelProfile .le. 0.) then
        PA=PA+Pi
        print*, "PA Adjust 1",PA
        if(PA .ge. 2.*pi) PA=PA-2.*pi    !Keep the PA between 0<PA<Phi
      endif
      print*, "Final PA Check", PA

c           TEMPORARY -- TRY AVERAGING THE RC AND SD PROFILES
c      if(nRings .gt. 3) then
c        AvgSD=sum(RadialProfile(1,3:nRings-1))/(nRings-3)
c        AvgRC=sum(RadialProfile(2,3:nRings-1))/(nRings-3)
c        AvgR=sum(RadialProfile(0,3:nRings-1))/(nRings-3)
c
c        LastPoint=0
c100     continue
c        nUse=nRings-3-LastPoint
c        print*, "N Rings Used", nUse
c        if(nUse .gt. 1) then
c            LastPoint=LastPoint+1
c            print*, "Number use", nUse
c            sumX=sum(RadialProfile(0,3:nRings-LastPoint))
c            sumY=sum(RadialProfile(2,3:nRings-LastPoint))
c            sumX2=sum(RadialProfile(0,3:nRings-LastPoint)**2.)
c            sumXY=sum(RadialProfile(0,3:nRings-LastPoint)
c     &                    *RadialProfile(2,3:nRings-LastPoint))
c            print*, sumX,sumY,sumX2,sumXY
c            b=(sumXY-sumX*sumY/nUse)/(sumX2-sumX**2./nUse)
c            a=AvgRC-b*AvgR
c            print*, "Intercept And Slope",a, b
c            if( a .lt. 0.) goto 100
c            if(a .gt. VLims(2)) goto 100
c        else
c            a=RadialProfile(2,3)
c            b=0.
c        endif
c      else
c        print*, "nrings is", nrings
c        a=RadialProfile(2,nRings-1)
c        b=0
c      endif

c      j=0
c      do i=0,2
c        if(ModelerableRingSwitch(i) .eq. 1) then
c            RadialProfile(1,j)=(AvgSD+RadialProfile(1,j))/2.
c            RadialProfile(2,j)=(AvgRC+RadialProfile(2,j))/2.
c            RadialProfile(2,j)=a+b*RadialProfile(0,j)
c            j=j+1
c        endif
c      enddo

c       Finally deallocate the temporary radial profile
      DEALLOCATE(TempRadialProfile)

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



      end module
