cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for writing out a
c       data cube to a file.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FitOutputMod
      use DataCubeOutputsMod
      use DataCubeMod
      use BeamMod
      use SoFiACatalogueMod
      use PipelineGlobals

      use ParameterVectorToTiltedRingMod
      use TiltedRingGenerationMod
      use FillDataCubeWithTiltedRingMod
      use CubeKernelConvolutionMod

      use UnitConvertMod

      use ParameterVectorToTiltedRingMod

      PROCEDURE(GeneralOutputInterface),POINTER :: OutputFit =>null()

      ABSTRACT INTERFACE
        subroutine GeneralOutputInterface(CatItem)
            import :: CatalogueItem
            Type(CatalogueItem),INTENT(IN) :: CatItem
        END subroutine GeneralOutputInterface
      END INTERFACE


      contains

ccccc
c
      subroutine OutputBestFit_Simple(CatItem)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      character(700) OutputFolder, Script


c           Make the output folder for the file
      call MakeModelOutputFolder(OutputFolder,CatItem)
c           Make a PSF convolved cube for the best fitting model
      call OutputCube(OutputFolder,trim(CatItem%ObjName))
c           Make a residual cube and output it
      call OutputResidualCube(OutputFolder,trim(CatItem%ObjName))
c           Output the tilted ring parameters
c      call OutputTiltedRingParams(OutputFolder)

c           Output the standard model output from  the proto pipeline
      call StandardModelOutput(OutputFolder,CatItem,2)

c           Write out the model moment maps
c      call WriteModelMaps(OutputFolder,CatItem)

c           Write out the model PV maps
c      call WritePVMaps(OutputFolder,CatItem)

c           Write out a file containing the set of fitting flags
      call WriteFlagFile(OutputFolder,CatItem)

c           For testing purposes, write out a file containing the
c               first fit -- for that we need to get the TR parameters
c               from the first fit parameter vector as well as the initial estimate
c      call ParamToTiltedRing(PV_FirstFit,ModelTiltedRing
c     &          ,TR_FittingOptions)
c      call StandardModelOutput(OutputFolder,CatItem,1)

c      print*, "Output initial Param",PVIni%Param(0:PVIni%nParams-1)
      call ParamToTiltedRing(PVIni,ModelTiltedRing
     &          ,TR_FittingOptions)
      call StandardModelOutput(OutputFolder,CatItem,0)
 
c           Make a diagnostic plot
c      call GenerateDiagnosticPlot()

      call CorrectRADEC()

      return
      end subroutine
ccccc


ccccc
      subroutine CorrectRADEC()
      use PipelineGlobals
      implicit none
      integer i, j,len,indx
      character(2000) PlotCmd
      character(8) ValStr
      character(500) CodePath,PythonPath
      character(9) CodeDelim

      write(ValStr,'(I4)') Version

      call GetArg(0,CodePath)

      CodeDelim="P"
      indx = SCAN(trim(CodePath),CodeDelim,.True.)
      PythonPath = CodePath(1:indx-2)
     &          //"/src/PythonScripts/GeometryFix.py "


      PlotCmd="python3.9 "//trim(PythonPath)//" "
     &          //trim(GalaxyDict%OutputFolder)//" "
     &          //trim(GalaxyDict%GalaxyName)//" "
     &          //trim(ValStr)//" "

      PlotCmd=trim(PlotCmd)//" "
     &          //trim(GalaxyDict%DataCubeFile)//" "


      print*, trim(PlotCmd)
      call system(trim(PlotCmd))

      return
      end subroutine
ccccc

ccccc
c
      subroutine MakeModelOutputFolder(OutputFolder,CatItem)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(700),INTENT(INOUT):: OutputFolder
      character(700) Script
c           Make the output folder for the file
      OutputFolder=trim(MainOutputFolder)//trim(CatItem%ObjName)
      print*, "Moving outputs to ", trim(MainOutputFolder)
      print*, "Making an output folder", trim(OutputFolder)
     &          ,trim(CatItem%ObjName)
      Script="mkdir "//trim(OutputFolder)
      call system(Script)

      GalaxyDict%OutputFolder=OutputFolder
      return
      end subroutine
ccccccc


cccccc
c
      subroutine OutputCube(OutputFolder,Name)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(*),INTENT(IN) :: Name
      real pixelarea
      character(500) CubeFile,BaseCubeFile
      character(3) VersStr
      character(10) format_string
      real BeamArea

c           The best fitting tilted ring model should have been made in Galaxy fit, so we
c           don't need to do the conversion of a parameter vector to tilted ring parameters

c       Make the tilted ring model
      pixelarea=abs(ModelDC%DH%PixelSize(0)*ModelDC%DH%PixelSize(1))
c      BeamArea=ObservedBeam%BeamMajorAxis*ObservedBeam%BeamMinorAxis
      BeamArea=2.*Pi*abs((ObservedBeam%BeamSigmaVector(0)
     &                      *ObservedBeam%BeamSigmaVector(1)))


      call BuildTiltedRingModel(ModelTiltedRing,idum)
c       Create the point-source data cube
      call FillDataCubeWithTiltedRing(ModelDC,ModelTiltedRing)
c      print*, "Filled DC", sum(ModelDC%Flux),pixelarea
c        Convolve the cube with the beam
c           Note that it is assumed that the real beam kernel has already been calculated
      call CubeBeamConvolution(ModelDC,ObservedBeam)
c           Because the cube is in units of Jy/pixel, convert back to Jy/beam
      ModelDC%Flux=ModelDC%Flux*BeamArea
c       Write the cube to a fits file
      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,format_string) Version

      BaseCubeFile=trim(GalaxyDict%GalaxyName)//"_AverageModel_v"
     &      //trim(VersStr)//".fits"
      CubeFile=trim(OutputFolder)//"/"//trim(BaseCubeFile)

      ModelDC%DH%FType='intensity'
      ModelDC%DH%PixelSize=ModelDC%DH%PixelSize/3600.

      call WriteDataCubeToFITS(ModelDC,ObservedBeam,CubeFile
     &              ,Name)

      GalaxyDict%BestFitCubeFile=trim(BaseCubeFile)



      return
      end subroutine
ccccc


cccccc
c
      subroutine OutputResidualCube(OutputFolder,Name)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(*),INTENT(IN)::Name
      character(500) CubeFile
      Type(DataCube) ResidCube
      real BeamArea

c      The model DC and the actual cube should already be in memory
      ResidCube=ModelDC
c       The Observed Cube is still in Jy/pixel so convert back to Jy/beam
c      BeamArea=abs(ObservedBeam%BeamMajorAxis
c     &              *ObservedBeam%BeamMinorAxis)
      BeamArea=2.*Pi*abs((ObservedBeam%BeamSigmaVector(0)
     &                      *ObservedBeam%BeamSigmaVector(1)))

c      print*, "Output beam Area", BeamArea,sum(ObservedDC%Flux)
c     &              ,sum(ObservedDC%Flux*BeamArea)
      ObservedDC%Flux=ObservedDC%Flux*BeamArea

c       Calculate the difference cube
      ResidCube%Flux=ObservedDC%Flux-ModelDC%Flux
      print*, "Flux Checks", sum(ObservedDC%Flux)
     &              ,sum(ModelDC%Flux),sum(ResidCube%Flux)

c       Write the cube to a fits file
      CubeFile=trim(OutputFolder)//"/"//"DifferenceCube.fits"
      call WriteDataCubeToFITS(ResidCube,ObservedBeam,CubeFile,Name)

      return
      end subroutine
ccccc

ccccc
c
      subroutine OutputTiltedRingParams(OutputFolder)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(500) ParamsFile

      integer i
      real SDConv, BeamArea,BeamPixels
      real CentPix(0:1),PAOut
      real SOut,SDTemp

c       The surface density units of the TR fitter are natively Jy/pixel.
c           The standard TR output should be Jy km/s arcsecond^-2.
c           Set the SDConversion factor to
      SDConv=abs(ObservedDC%DH%ChannelSize)
     &          /(abs(ObservedDC%DH%PixelSize(0))
     &          *abs(ObservedDC%DH%PixelSize(1)))



      ParamsFile=trim(OutputFolder)//"/"//"TiltedRingParams.txt"
      open(10,file=trim(ParamsFile),status='replace')
      write(10,*) "#    Best fitting Tilted Ring Params"
      write(10,*) "#    Generated by the WRKP"
      write(10,*) "#    Fitting options"
      write(10,*) "#    Core Code algorithm"
      if(PFlags%CoreCodeSwitch .eq. 1) then
        write(10,*) "WRKP internal algorithm"
      elseif(PFlags%CoreCodeSwitch .eq. 2) then
        write(10,*) "3DBarolo"
      endif
      write(10,*) "#    Fitting options"
      write(10,*) "cmode= ", ModelTiltedRing%cmode
      write(10,*) "cdens= ", ModelTiltedRing%CloudBaseSurfDens
      write(10,*) "rings/beam= ", TR_FittingOptions%nRingsPerBeam

      write(10,*) "#    Best Fit"
      write(10,*) "Likelihood=", PVModel%BestLike

      write(10,*) "#    Rmid (arcsec)    Rwidth(arcsec)"
     &          //"    Xcent (pixels)"
     &          //"    Ycent (pixels)"
     &          //"    Inc (degrees)"
     &          //"    PA (degrees)   VSys (km/s)"
     &          //"    VRot (km/s)    VRad (km/s)"
     &          //"    Vvert (km/s)    VDisp (km/s)"
     &          //"    dvdz (km/s)    Sigma (Jy km/s arcsec^-2"
     &          //"    z0  ('')   zGradStart ('')   "


      do i=0, ModelTiltedRing%nRings-1
        PAOut=(ModelTiltedRing%R(i)%PositionAngle*180./Pi)
        PAOut=PAOut-90.
100     continue
        if(PAOut .lt. 0.) then
            PAOut=PAOut+360.
            goto 100
        elseif(PAOut .gt. 360.) then
            PAOut=PAOut-360.
            goto 100
        endif
        if(PFlags%Linear_Log_SDSwitch .eq. 0) then
            SDTemp=ModelTiltedRing%R(i)%SigUse
        elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
            SDTemp=10.**(ModelTiltedRing%R(i)%SigUse)
        endif


c        print*, "Final SD and VRot"
c     &          ,ModelTiltedRing%R(i)%Rmid
c     &          ,ModelTiltedRing%R(i)%VRot
c     &          ,ModelTiltedRing%R(i)%Sigma,SDConv

        write(10,*) ModelTiltedRing%R(i)%Rmid
     &                  *abs(ObservedDC%DH%PixelSize(0))
     &          , ModelTiltedRing%R(i)%Rwidth
     &                  *abs(ObservedDC%DH%PixelSize(0))
     &          , ModelTiltedRing%R(i)%CentPos
     &          , ModelTiltedRing%R(i)%Inclination*180./Pi
     &          , PAOut
     &          , ModelTiltedRing%R(i)%VSys, ModelTiltedRing%R(i)%VRot
     &          , ModelTiltedRing%R(i)%VRad, ModelTiltedRing%R(i)%Vvert
     &          , ModelTiltedRing%R(i)%VDisp, ModelTiltedRing%R(i)%dvdz
     &          , SDTemp*SDConv
     &          , ModelTiltedRing%R(i)%z0
     &          , ModelTiltedRing%R(i)%zGradiantStart


      enddo
      close(10)

      return
      end subroutine
ccccc


ccccccc
c
      subroutine StandardModelOutput(OutputFolder,CatItem
     &                  ,FitNum)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500),INTENT(IN):: OutputFolder
      integer,INTENT(IN) :: FitNum
      character(500) ParamsFile
      character(10) VersStr, format_string
      character(100) DateStr

      character(300) OutStr

      character(8) ValStr,ErrStr
      character(18) PreambleStr
      character(20) RadialProfStr(6)

      character(8) date
      character(10) time
      character(5) zone
      
      integer timeArray(8)
      integer i
      real PAOut
      real BeamPixels,SDConv1
      real SDTemp
      real CentValAS(2),RA,DEC

      call date_and_time(date,time,zone,timeArray)




c       The machine SD units are Jy/pixel and we want them in M_sol/pc^2
c           First get the conversion to Jy/arcsec^2
c           This calculation requires that the pixelsize arrays be in arcseconds
      SDConv1=abs(ObservedDC%DH%ChannelSize)
     &          /(abs(ObservedDC%DH%PixelSize(0))
     &          *abs(ObservedDC%DH%PixelSize(1)))

      print*, "SD Convert factor", SDConv1


      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,format_string) Version

c      print*, "Obj Name", trim(CatItem%ObjName)

      call NameParamFile(OutputFolder,CatItem,FitNum,ParamsFile)
c      print*,trim(ParamsFile)

      open(10,file=trim(ParamsFile),status='replace')

c       Write the Preamble
      write(10,'(a,a)') "Object:  ",trim(CatItem%ObjName)
      write(10,'(a,a)') "Source:  ",trim(SCatLocal%SourceName)
      DateStr="Date:    "//date(7:8)//"-"//date(5:6)//"-"//date(1:4)
      write(10,'(a)') trim(DateStr)
      write(10,'(a,a)') "Version: ",trim(VersStr)

c       Write out the noise measurements
c           First the RMS--which needs to be converted
      write(10,*) " "
      PreambleStr="RMS (mJy/beam)"
      write(ValStr, '(F8.2)')ObservedDC%DH%Uncertainty
     &              *ObservedBeam%BeamAreaPixels*1000.
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       Next the integrated S/N
      PreambleStr="SN_Integrated "
      write(ValStr, '(F8.2)')ObservedDC%DH%SN_Int
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the peak S/N
      PreambleStr="SN_Peak "
      write(ValStr, '(F8.2)')ObservedDC%DH%SN_Peak
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the average S/N
      PreambleStr="SN_Avg "
      write(ValStr, '(F8.2)')ObservedDC%DH%SN_Avg
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the median S/N
      PreambleStr="SN_Median "
      write(ValStr, '(F8.2)')ObservedDC%DH%SN_Median
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)


c       Write out the geometric model parameters

      write(10,*) " "
      write(10,'(a)') "Geometry Paramters"
      write(10,'(a)') "Param Name            Value    Error"

c      print*, "Reference values",ObservedDC%DH%RefVal
c     &                  ,ObservedDC%DH%RefLocation
c       Get the RA and DEC value for the central point
      do i=0,1
        CentValAS(i+1)=(ModelTiltedRing%R(0)%CentPos(i)
     &              -ObservedDC%DH%RefLocation(i))
     &                  *ObservedDC%DH%PixelSize(i)
     &              +ObservedDC%DH%RefVal(i)
      enddo
c      print*, "Central Position AS", CentValAS
      call ArcSecToDegrees(CentValAS(1),RA)
      call ArcSecToDegrees(CentValAS(2),DEC)
c      print*, "Central Position Deg", RA,DEC

      do i=0, 7
        if( i.eq. 0) then
            PreambleStr="X_kin (pixels)"
            write(ValStr, '(F8.2)')ModelTiltedRing%R(0)%CentPos(0)
            write(ErrStr, '(F5.2)')0.0
        elseif(i .eq. 1) then
            PreambleStr="Y_kin (pixels)"
            write(ValStr, '(F8.2)')ModelTiltedRing%R(0)%CentPos(1)
            write(ErrStr, '(F5.2)')0.0
        elseif(i .eq. 2) then
            PreambleStr="RA_kin (degrees)"
            write(ValStr, '(F8.2)')RA
            write(ErrStr, '(F5.2)')0.0
        elseif(i .eq. 3) then
            PreambleStr="DEC_kin (degrees)"
            write(ValStr, '(F8.2)')DEC
            write(ErrStr, '(F5.2)')0.0
        elseif(i .eq. 4) then
            PreambleStr="Inc_kin (degrees)"
            write(ValStr, '(F5.2)')ModelTiltedRing%R(0)%Inclination
     &                  *180./Pi
            write(ErrStr, '(F5.2)')0.0
        elseif(i .eq. 5) then
            PreambleStr="PA_kin (degrees)"
            PAOut=(ModelTiltedRing%R(0)%PositionAngle*180./Pi)
            PAOut=PAOut-90.
 100        continue
            if(PAOut .lt. 0.) then
                PAOut=PAOut+360.
                goto 100
            elseif(PAOut .gt. 360.) then
                PAOut=PAOut-360.
                goto 100
            endif
            write(ValStr, '(F8.3)')PAOut
            write(ErrStr, '(F8.3)')0.0
        elseif(i .eq. 6) then
            PreambleStr="VSys_kin (km/s)"
            write(ValStr, '(F8.2)')ModelTiltedRing%R(0)%VSys
            write(ErrStr, '(F8.2)')0.0
        elseif(i .eq. 7) then
            PreambleStr="VDisp_kin (km/s)"
            write(ValStr, '(F8.2)')ModelTiltedRing%R(0)%VDisp
            write(ErrStr, '(F8.2)')0.0
        endif

        OutStr=PreambleStr//"    "//ValStr//" "//ErrStr
        write(10,'(a)') trim(OutStr)
      enddo
c       Write out the radial profiles
      write(10,*)" "

      write(10,'(a)') " Radial Profiles"
      write(10,'(a)')" Rad        VROT_kin    e_VRot_kin"
     &           //"       e_VROT_kin,inc      SD_kin      e_SD_kin"
      write(10,'(a)') "('')       (km/s)         (km/s)     "
     &              //"(km/s)             (Msol/pc^2)      (Msol/pc^2)"

      do i=0,ModelTiltedRing%nRings-1
c       Get the Surface density in units of M_sol/pc^2
        if(PFlags%Linear_Log_SDSwitch .eq. 0) then
            SDTemp=ModelTiltedRing%R(i)%SigUse
        elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
            SDTemp=10.**(ModelTiltedRing%R(i)%SigUse)
        endif

        SDTemp=SDTemp*SDConv1
        call JyAS2_To_MSolPc2(SDTemp,SDTemp)

        print*, "Rmid",ModelTiltedRing%R(i)%Rmid
     &                  ,abs(ObservedDC%DH%PixelSize(0))
        write(RadialProfStr(1), '(F8.2)')ModelTiltedRing%R(i)%Rmid
     &                  *abs(ObservedDC%DH%PixelSize(0))
        write(RadialProfStr(2), '(F8.2)')ModelTiltedRing%R(i)%VRot
        write(RadialProfStr(3), '(F5.2)')0.
        write(RadialProfStr(4), '(F5.2)')0.
        write(RadialProfStr(5), '(G9.2)')SDTemp
        write(RadialProfStr(6), '(F5.2)')0.
       OutStr=RadialProfStr(1)//RadialProfStr(2)//RadialProfStr(3)
     &          //RadialProfStr(4)//RadialProfStr(5)//RadialProfStr(6)
        write(10,'(a)') trim(OutStr)

      enddo

      close(10)

      return
      end subroutine
ccccccc

ccccc
      subroutine NameParamFile(OutputFolder,CatItem
     &              ,FitNum,ParamsFile)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500),INTENT(IN):: OutputFolder
      integer,INTENT(IN) :: FitNum
      character(500),INTENT(INOUT):: ParamsFile
      
      character(8) VersStr,format_string

      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,trim(format_string)) Version
    

      if(FitNum .eq. 2) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_AvgModel_v"
     &          //trim(VersStr)//".txt"
      elseif(FitNum .eq. 1) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_FirstModel_v"
     &          //trim(VersStr)//".txt"
      elseif(FitNum .eq. 0) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_IniEstimate_v"
     &          //trim(VersStr)//".txt"
      endif
      ParamsFile=trim(OutputFolder)//"/"
     &          //trim(GalaxyDict%BestFitModelFile)
      return
      end subroutine
ccccccc


cccccccc
      subroutine GenerateDiagnosticPlot()
      use PipelineGlobals
      implicit none

      integer i, j,len
      character(2000) PlotCmd

      character(8) ValStr
      character(500) MaskTemp,BaseMaskFile
      character(500) CodePath,PythonPath
      character(9) CodeDelim
      integer indx

      write(ValStr,'(I4)') Version

c           The MaskFile string needs to be adjusted in case there are any spaces
c           First get the length of the trimmed file name
c       Initialize the new string
      if(MomentMapSwitch .eq. 0) then
        BaseMaskFile=trim(GalaxyDict%MaskFile)
        MaskTemp=GalaxyDict%MaskFile(1:1)
      elseif(MomentMapSwitch .eq. 3) then
        BaseMaskFile=trim(SecondaryMaskFile)
        MaskTemp=SecondaryMaskFile(1:1)
      endif
      len=len_trim(BaseMaskFile)
      print*, "hmmm", len
      print*, trim(BaseMaskFile)
c           Loop through the trimmed file name
      do i=2,len
c           Replace spaces with '\ ' for terminal commands
        if(BaseMaskFile(i:i) .eq. " ") then
            MaskTemp=trim(MaskTemp)//"\"    !   Only add a '\' here as the next trim would remove the space
        else
            j=len_trim(MaskTemp)
c           Now double check if the last character is a '\', which needs a space to trail it
            if( MaskTemp(j:j) .eq. "\") then
                MaskTemp=trim(MaskTemp)//" "//BaseMaskFile(i:i)
            else
                MaskTemp=trim(MaskTemp)//BaseMaskFile(i:i)
            endif
        endif
      enddo


c      call GetCWD(CodePath)
      call GetArg(0,CodePath)
      print*, "Current Code Path ", trim(CodePath)

      CodeDelim="P"
      indx = SCAN(trim(CodePath),CodeDelim,.True.)
      PythonPath = CodePath(1:indx-2)
     &          //"/src/PythonScripts/BestFitDiagnosticPlot.py "

      print*, "Trial Python Path ", indx,trim(PythonPath)

      PlotCmd="python3.9 "//trim(PythonPath)//" "
     &          //trim(GalaxyDict%OutputFolder)//" "
     &          //trim(GalaxyDict%GalaxyName)//" "
     &          //trim(ValStr)//" "

      if(MomentMapSwitch .eq. 0) then
        PlotCmd=trim(PlotCmd)//" "
     &          //trim(GalaxyDict%DataCubeFile)//" "
      elseif(MomentMapSwitch .eq. 3) then
        PlotCmd=trim(PlotCmd)//" "
     &          //trim(SecondaryCubeFile)//" "
      endif

      PlotCmd=trim(PlotCmd) //" "//trim(MaskTemp)


      write(ValStr,'(F8.2)') GalaxyDict%Distance
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)

      write(ValStr,'(F8.3)') GalaxyDict%Size
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)

      write(ValStr,'(F5.2)') GalaxyDict%logMass
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)
      print*, trim(PlotCmd)
      call system(trim(PlotCmd))

    

      return
      end subroutine
cccccccc


cccccc
c
      subroutine WriteModelMaps(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(1) MomNum
      integer i
      character(500) MomMapName
      character(8) ValStr

      write(ValStr,'(I1)') Version

      print*, "Writing Model maps"

      do i=0,2
        write(MomNum,'(I1)') i
        print*, MomNum

        MomMapName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_Mom"//trim(MomNum)
     &              //"_v"//trim(ValStr)
     &              //".fits"
        print*, trim(MomMapName)
        call WriteDCSliceToFITS(ModelMaps,ObservedBeam
     &              ,MomMapName,i,trim(CatItem%ObjName))
      enddo

      return
      end subroutine
ccccccc




cccccc
c
      subroutine WritePVMaps(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500) PVMapName
      character(20) Suffix
      integer i


      do i=0,1
        if(i .eq. 0) then
            Suffix="MajorAxisPV"
        elseif(i .eq. 1) then
            Suffix="MinorAxisPV"
        endif
        PVMapName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_"//trim(Suffix)//".fits"
      call WriteDCSliceToFITS(ModelPVMaps,ObservedBeam
     &              ,PVMapName,i,trim(CatItem%ObjName))
      enddo

      return
      end subroutine
ccccccc





cccccc
c
      subroutine WriteFlagFile(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500) FlagFileName
      character(8) ValStr

      write(ValStr,'(I1)') Version

      FlagFileName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_Flags_v"//trim(ValStr)//".txt"

      open(10,file=FlagFileName,status='replace')
      write(10,'(a)') "# Goodness Of Fit  type (1==chi, 3 == log10(chi)"
      write(10,*) PFlags%LikelihoodSwitch


      write(10,'(a)') "# Goodness Of Fit "
      write(10,*) GalaxyDict%GoodnessOfFit

      write(10,'(a)') "# RMS noise "
      write(10,*) ObservedDC%DH%Uncertainty

      write(10,'(a)') "# Total cube Cells "
      write(10,*) GalaxyDict%nCells

      write(10,'(a)') "# Normalized Goodness Of Fit "
      write(10,*) GalaxyDict%Norm_GoodnessOfFit

      write(10,'(a)') "# Size Flag (1-> <2beams, 2-> >10 beams)"
      write(10,*) GalaxyDict%Flags%SizeFlag

      write(10,'(a)') "# Pre-analysis center Flag (1-> no convergence)"
      write(10,*) GalaxyDict%Flags%CenterFlag

      close(10)

      return
      end subroutine
ccccccc

      end module
