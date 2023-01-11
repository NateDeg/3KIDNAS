mFlags=../makeflags  
include $(mFlags)

all: $(CompareCubeObj) install


LikelihoodFunctions.o:
	$(MakeFObj) LikelihoodFunctions.f

MaskCube.o:
	$(MakeFObj) MaskCube.f

CubeComparison.o: LikelihoodFunctions.o
	$(MakeFObj) CubeComparison.f

FullModelComparison.o: CubeComparison.o MaskCube.o
	$(MakeFObj) FullModelComparison.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
