mFlags=../makeflags
include $(mFlags)

VPATH :=$(BinLoc)

ProgramList= $(AllProgs)

all: $(ProgramList)


BootStrapSampler : BootStrapGenerator.o
	$(F77) $(FLAGS) -o BootStrapSampler  BootStrapGenerator.o $(AllObj)  $(FitsLibLoc) $(FFTW_LIB)

SingleGalaxyFitter : SingleGalaxyFitTests.o
	$(F77) $(FLAGS) -o SingleGalaxyFitter SingleGalaxyFitTests.o $(AllObj) $(FitsLibLoc) $(FFTW_LIB) 


install:
	mv $(ProgramList) ../../Programs


clean:
	rm -f $(ProgramList);
	rm -f *.o;


