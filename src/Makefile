CC=g++
CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math 
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L${CURDIR}/lib/ -Wl,-rpath,${CURDIR}/lib/
#
INCLUDE=-I${CURDIR}/ -I${CURDIR}/include/ -I${CURDIR}/include/R/
#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm 

OBJECTS= SFAmix_c.o

all: install SFAmix.so 

SFAmix.so: $(OBJECTS)
	$(CC) -shared SFAmix_c.o -o SFAmix.so $(LINKER) $(LIBS)
SFAmix_c.o:  SFAmix_c.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  SFAmix_c.cpp
myHeader.o: myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	bash ./install.sh
clean:
	rm *.o
