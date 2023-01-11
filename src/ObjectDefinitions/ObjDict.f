cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definition of a SoFiA
c       catalogue
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GalaxyObjDictMode

      implicit none


      Type GalaxyFlags
        integer SizeFlag        !0 == Acceptable Size
                                !1 == Object smaller than size limit (See PipelineGlobals for size limits in beams)
                                !2 == Object large enough that simple flat disks may be insufficient.

        integer CenterFlag      !0 == Center is fine
                                !1 == Iterative center did not converge

      end Type

      Type GalaxyDictionary
        character(100) GalaxyName
        character(500) BestFitModelFile
        character(450) OutputFolder
        character(500) BestFitCubeFile
        character(500) DataCubeFile
        character(500) MaskFile

        real Distance
        real logMass
        real Size
        real RMS

        real GoodnessOfFit
        integer nCells
        real Norm_GoodnessOfFit

        Type(GalaxyFlags) Flags

      end Type


      end module
