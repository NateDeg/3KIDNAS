cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definition of a SoFiA
c       catalogue
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SoFiACatalogueMod

      implicit none

      Type CatalogueItem
        integer ObID
        real ZCent
        real RA,DEC
        real XMin,XMax,YMin,YMax,ZMin,ZMax
        real EllipseMaj,EllipseMin,EllipsePA
        real kinPA
        real RMS,w20,w50,centralFreq
        integer nChannels,n_LOS
        character(100) ObjName
        real EllipseInc

      end Type

      Type SoFiACat
        integer nObjects,SofiaSwitch
        character(500) PrepScript
        Type(CatalogueItem), ALLOCATABLE :: Objects(:)
        character(100) SourceName
      end Type





      contains


cccccc
      subroutine SoFiACatalogueAllocation(SC)
      implicit none
      Type(SoFiACat), INTENT(INOUT) :: SC
c
      ALLOCATE(SC%Objects(0:SC%nObjects-1))       !Keep indexing consistent with C
      return
      end subroutine
cccccccc

cccccc
      subroutine SoFiACatalogueDeAllocation(SC)
      implicit none
      Type(SoFiACat) SC
c
      DEALLOCATE(SC%Objects)
      return
      end subroutine
cccccccc



      end module
