CC=g++
CFLAGS=-g -O3 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math

INCLUDE=-isystem${CURDIR}/ -isystem${CURDIR}/include/ -isystem${CURDIR}/include/R/ -w

OBJECTS= SFAmix_c.o

all: install SFAmix.so clean

SFAmix.so: $(OBJECTS)
	$(CC) -shared SFAmix_c.o -o SFAmix.so
SFAmix_c.o:  SFAmix_c.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  SFAmix_c.cpp
myHeader.o: myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	./install.sh
clean:
	rm *.o
