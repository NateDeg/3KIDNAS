mFlags=../makeflags  
include $(mFlags)

all: $(MakeMomentMapsObj) install



clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
