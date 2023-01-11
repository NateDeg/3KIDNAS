mFlags=makeflags 
include $(mFlags)


all: LocalLib Programs

LocalLib: 
	cd StandardMath; make;
	cd ObjectDefinitions; make;
	cd UnitConversions; make;
	cd GeneralMinimizationRoutines; make;
	cd TiltedRingModelGeneration; make;
	cd TiltedRingToDataCube; make;
	cd ConvolveCube; make;
	cd GenerateNoiseCubes; make;
	cd ParameterToTiltedRingInterface; make;
	cd MomentMaps; make;
	cd ParticleGeometryTransformations; make;
	cd Globals; make;
	cd CompareCubes; make;
	cd GalaxyAnalysis; make;
	cd Outputs; make;
	cd Inputs; make;
	cd PreAnalysis; make;
	cd PostAnalysis; make;
	cd BootstrapSampler; make;
	cd ProgramMains; make;


Programs:
	cd $(LocBinLoc); make; make install



clean:
	cd StandardMath; make clean;
	cd ObjectDefinitions; make clean;
	cd UnitConversions; make clean;
	cd GeneralMinimizationRoutines; make clean;
	cd TiltedRingModelGeneration; make clean;
	cd TiltedRingToDataCube; make clean;
	cd ConvolveCube; make clean;
	cd GenerateNoiseCubes; make clean;
	cd ParameterToTiltedRingInterface; make clean;
	cd MomentMaps; make clean;
	cd ParticleGeometryTransformations; make clean;
	cd Inputs; make clean;
	cd Globals; make clean;
	cd PreAnalysis; make clean;
	cd PostAnalysis; make clean;
	cd CompareCubes; make clean;
	cd GalaxyAnalysis; make clean;
	cd Outputs; make clean;
	cd ProgramMains; make clean;


	cd $(LocBinLoc); rm *.o;
	cd $(LocModLoc); rm *.mod
	cd ../Programs; rm $(AllProgs)

	

