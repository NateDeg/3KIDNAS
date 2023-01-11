cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines necessary to
c       get the moment maps of the observed cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GetMomentMapsMod
      use MakeMomentMapsMod
      use InputUnitConversionsMod
      use CommonConsts
      use DataCubeInputMod

      implicit none
      contains

ccccc
      subroutine ConstructMomentMaps(MaskedObsDC,ObsMaps)
      implicit none
      Type(DataCube), INTENT(IN) :: MaskedObsDC
      Type(DataCube), INTENT(INOUT) :: ObsMaps

      integer i

      ObsMaps%DH=MaskedObsDC%DH
      ObsMaps%DH%nChannels=3
      call AllocateDataCube(ObsMaps)
c           Calculate the maps from the cube
      call MakeMomentMaps(MaskedObsDC, ObsMaps) !/src/MomentMaps/CalculateMomentMaps.f
 

      return
      end subroutine
cccccc


cccc
c       WALLABY cubelets have moment1 and 2 maps in freqency terms.  We will need to convert to
c           km/s for analysis.
      subroutine FreqencyMomMapsConvert(Maps)
      implicit none
      Type(DataCube) Maps

      real TempMom1, TempMom2,z
      integer i,j
      real NanMaker

      NanMaker=0.
c       Loop over all pixels
      do i=0, Maps%DH%nPixels(0)-1
        do j=0, Maps%DH%nPixels(0)-1
c           Convert the moment 0 map to Jy/beam
            Maps%Flux(i,j,0)=Maps%Flux(i,j,0)/Maps%DH%delFreq
c            Set TempMom1 to the moment 1 frequency
            TempMom1=Maps%Flux(i,j,1)
c           Set TempMom2 to the moment1 + moment2 frequency in order to eventually get dV
            TempMom2=Maps%Flux(i,j,2)+Maps%Flux(i,j,1)
c               Convert the temporary mom1 pixel to a velocity
            call RedshiftCalc(z,Maps%DH%RestFreq,TempMom1)
            TempMom1=z*Lightspeed
c               Convert the temporary mom2 pixel to a velocity
            call RedshiftCalc(z,Maps%DH%RestFreq,TempMom2)
            TempMom2=z*Lightspeed
c               Set TempMom2 to be the dV between the moment1 velocity and the one just calculated.
            TempMom2=abs(TempMom2-TempMom1)
c           Set the moment1 and moment2 map fluxes to these values
            Maps%Flux(i,j,1)=TempMom1
            Maps%Flux(i,j,2)=TempMom2
            if(Maps%Flux(i,j,0) .le. 0.) then
                Maps%Flux(i,j,1:2)=0./NanMaker
            endif
        enddo
      enddo
      return
      end subroutine
ccccccc




      end module
