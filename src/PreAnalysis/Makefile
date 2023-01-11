mFlags=../makeflags  
include $(mFlags)

all: $(PreAnalysisObj) install


EstimateRadialProfiles.o: 
	$(MakeFObj) EstimateRadialProfiles.f

VelProfileAnalysis.o: 
	$(MakeFObj) VelProfileAnalysis.f

EstimateShape.o: 
	$(MakeFObj) EstimateShape.f

GetMomentMaps.o: 
	$(MakeFObj) GetMomentMaps.f

ModellingInitializations.o: 
	$(MakeFObj) ModellingInitializations.f

EstimateCubeNoise.o: 
	$(MakeFObj) EstimateCubeNoise.f


InitialAnalysis_PreProcessing.o: GetMomentMaps.o VelProfileAnalysis.o EstimateCubeNoise.o
	$(MakeFObj) InitialAnalysis_PreProcessing.f

InitialAnalysis.o: EstimateCubeNoise.o EstimateShape.o EstimateRadialProfiles.o VelProfileAnalysis.o ModellingInitializations.o
	$(MakeFObj) InitialAnalysis.f

PreAnalysis.o: ModellingInitializations.o InitialAnalysis_PreProcessing.o InitialAnalysis.o
	$(MakeFObj) PreAnalysis.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
