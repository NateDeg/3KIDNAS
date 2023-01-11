mFlags=../makeflags  
include $(mFlags)

all: $(UnitConvObj) install

Particle.o:
	$(MakeFObj) Particle.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
