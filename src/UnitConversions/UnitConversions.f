cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines to do different unit
c       conversions needed in the code
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module UnitConvertMod
      use CommonConsts

      contains
ccccccc
c           Convert a length in arcseconds to degrees
      subroutine ArcSecToDegrees(L_AS,L_Deg)
      implicit none
      real,INTENT(IN) :: L_AS
      real,INTENT(INOUT) :: L_Deg
      L_Deg=L_AS/Degree_To_AS
      return
      end subroutine
cccccccc

ccccccc
c           Convert a length in degrees to arcseconds
      subroutine DegreesToArcSec(L_Deg,L_AS)
      implicit none
      real,INTENT(INOUT) :: L_AS
      real,INTENT(IN) :: L_Deg
      L_AS=L_Deg*Degree_To_AS
      return
      end subroutine
cccccccc


ccccccc
c           Convert an angle in degrees to radians
      subroutine DegreesToRadians(L_Deg,L_Rad)
      implicit none
      real,INTENT(INOUT) :: L_Rad
      real,INTENT(IN) :: L_Deg
      L_Rad=L_Deg*Pi/180.
      return
      end subroutine
cccccccc

ccccccc
c           Convert a length in degrees to arcseconds
      subroutine RadiansToDegrees(L_Rad,L_Deg)
c
      implicit none
      real,INTENT(IN) :: L_Rad
      real,INTENT(INOUT) :: L_Deg
      L_Deg=L_Rad*180./Pi
      return
      end subroutine
cccccccc


cccccc
c       Convert a length from kilometers to meters
      subroutine KmToM(LKM,LM)
      implicit none
      real,INTENT(IN) :: LKM
      real,INTENT(INOUT) :: LM
      LM=LKM*1000.
      return
      end subroutine
cccccc

cccccc
c       Convert a length from meters to kilometers
      subroutine MToKM(LM,LKM)
      implicit none
      real,INTENT(IN) :: LM
      real,INTENT(INOUT) :: LKM
      LKM=LM/1000.
      return
      end subroutine
cccccc


ccccc
c
c       Convert a length to pixels.  It assumes that the pixel units and
c           length units are the same
      subroutine Length_PixelConversion(L,L_Pixel,PixelSize)
      implicit none
      real, INTENT(IN) :: L,PixelSize
      real, INTENT(INOUT) :: L_Pixel
c       The absolute is for 'negative' pixel sizes in some datacube headers
      L_Pixel=L/abs(PixelSize)
      return
      end subroutine
ccccccc




cccccc
c       Calculate the 'optical' redshift using a frequency and rest frequency
      subroutine RedshiftCalc(z,RestFreq,Freq)
      implicit none
      real,INTENT(IN) :: RestFreq,Freq
      real,INTENT(OUT) :: z
      real numer,denom
c       Using the optical definition
      numer=RestFreq-Freq
      denom=Freq
      z=numer/denom
      return
      end subroutine
cccccccc

cccccc
c       Calculate the conversion by Jy/arcsec^2 to M_sol/pn^2
      subroutine JyAS2_To_MSolPc2(SD_J,SD_M)
      implicit none
      real,intent(IN) :: SD_J
      real,intent(INOUT) ::SD_M
      SD_M=SD_J/JyAS_To_MsolPC
      print*, "SD ConvCheck", SD_M, SD_J,JyAS_To_MsolPC

      return
      end subroutine
cccccc

cccccc
c       Calculate the conversion by Jy/arcsec^2 to M_sol/pn^2
      subroutine MSolPc2_To_JyAS2(SD_M,SD_J)
      implicit none
      real,intent(IN) :: SD_M
      real,intent(INOUT) ::SD_J
      SD_J=SD_M*JyAS_To_MsolPC

      return
      end subroutine
cccccc


      end module
