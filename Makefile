

CXX = /opt/local/bin/g++
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize

dcd2txt: dcd2txt.cpp psflib.h psflib.cpp dcdlib.h dcdlib.c stringlib.h stringlib.cpp xyzlib.h xyzlib.c
	$(CXX) -c dcd2txt.cpp stringlib.cpp xyzlib.c psflib.cpp dcdlib.c $(OPTS) 
	$(CXX) dcd2txt.o stringlib.o xyzlib.o psflib.o dcdlib.o $(OPTS) -lm -o dcd2txt.x

clean:
	rm -f *.o
