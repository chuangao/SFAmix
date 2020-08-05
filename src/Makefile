CC=g++
unamestr := $(shell uname)
ifeq ($(unamestr),Linux)
	CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math -fopenmp -std=c++11
endif
ifeq ($(unamestr),Darwin)
	CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math -std=c++11
endif
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L${CURDIR}/lib/ -Wl,-rpath,${CURDIR}/lib/

INCLUDE=-isystem${CURDIR}/ -isystem${CURDIR}/include/ -isystem${CURDIR}/include/R/ -w
#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm 


all: install myHeader.o SFAmix.so clean

SFAmix.so: SFAmix_c.o
	$(CC) -shared SFAmix_c.o -o SFAmix.so $(LINKER) $(LIBS)
SFAmix_c.o:  SFAmix_c.cpp myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  SFAmix_c.cpp
myHeader.o: myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	./install.sh
clean:
	rm *.o
