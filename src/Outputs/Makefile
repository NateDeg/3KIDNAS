mFlags=../makeflags  
include $(mFlags)

all: $(OutputObj) install

DataCubeOutputs.o:
	$(MakeFObj) DataCubeOutputs.f

FitOutput.o: DataCubeOutputs.o
	$(MakeFObj) FitOutput.f



clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
