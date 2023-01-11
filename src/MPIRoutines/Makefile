mFlags=../makeflags  
include $(mFlags)

all: $(MPIObj) install

MPI_CharSplit.o:
	$(MakeFObj) MPI_CharSplit.f

MPI_IntSplit.o:
	$(MakeFObj) MPI_IntSplit.f

MPISplit.o: MPI_CharSplit.o MPI_IntSplit.o
	$(MakeFObj) MPISplit.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
