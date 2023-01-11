mFlags=../makeflags  
include $(mFlags)

all: $(MinimizationObj) install


DownhillSimplex.o:
	$(MakeFObj) DownhillSimplex.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
