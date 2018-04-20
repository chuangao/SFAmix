
Rhome = $$(/bin/bash R RHOME)

CC=g++
CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp

INCLUDE=-isystem${CURDIR}/ -isystem${CURDIR}/include/ -isystem${CURDIR}/include/R/ -w

all: install SFAmix.so clean

SFAmix.so: SFAmix_c.o
	$(CC) -shared SFAmix_c.o -o SFAmix.so
SFAmix_c.o:  SFAmix_c.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  SFAmix_c.cpp
myHeader.o: myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	mkdir -p ./include/R
	cp -r $(Rhome)/include/* ./include/R/
clean:
	rm *.o
