

      module CommonConsts
      implicit none

      real, parameter :: Pi=4.*atan(1.)
      real, parameter :: HIRestFreq=1.42040575179E+09     !The rest frequency of 21 cm line in Hz
      real, parameter :: lightspeed=2.99792458e5      !The speed of light in a vacuum in km/s

      real, parameter :: Degree_To_AS=3600.         !   Conversion of degrees to arcseconds
      real, parameter :: Radian_To_AS=206265.       !   Conversion of radians to arcseconds
      real,parameter :: JyAS_To_MsolPC=1.24756e+20/(6.0574E5*1.823E18
     &                      *(2.*Pi/log(256.)))

      real, parameter :: H0=70.               !The value of the Hubble constant in units of km/s/Mpc


      end module CommonConsts
