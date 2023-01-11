mFlags=../makeflags  
include $(mFlags)

all: $(PostAnalysisObj) install


PostGalaxyAnalysis.o:
	$(MakeFObj) PostGalaxyAnalysis.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
