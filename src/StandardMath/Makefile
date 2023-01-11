mFlags=../makeflags  
include $(mFlags)

all: $(StandardMathObj) install

BasicConstants.o:
	$(MakeFObj) BasicConstants.f

random.o:
	$(MakeFObj) random.f

FullCircTrig.o: BasicConstants.o
	$(MakeFObj) FullCircTrig.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
