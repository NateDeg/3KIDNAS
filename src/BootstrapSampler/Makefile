mFlags=../makeflags  
include $(mFlags)

all: $(BootstrapObj) install


CubeDifference.o:
	$(MakeFObj) CubeDifference.f
	
PhysCoordTransform.o:
	$(MakeFObj) PhysCoordTransform.f

GenerateBootstrap.o: CubeDifference.o PhysCoordTransform.o
	$(MakeFObj) GenerateBootstrap.f
	

FlippingBootstrap.o: CubeDifference.o PhysCoordTransform.o GenerateBootstrap.o
	$(MakeFObj) FlippingBootstrap.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
