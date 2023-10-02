CC=g++
CFLAGS=-ggdb -pthread
CPPFLAGS=-O2 -pthread `pkg-config --cflags freetype2` `pkg-config --cflags libpipewire-0.3` -I esutil `sdl2-config --cflags`
LDLIBS=-ljack -lrt -lfftw3_threads -lfftw3 -lX11 -lGL -lGLEW -lm -L esutil -lutil `pkg-config --libs freetype2` `pkg-config --libs libpipewire-0.3`
LDFLAGS=-O2 -pthread

all:fftcom pwfftcom

fftcom:fftcom.o vsample.o fftcodec.o

fftcom.o:fftcom.cc

vsample.o:vsample.cc vsample.h

fftcodec.o:fftcodec.cc fftcodec.h vsample.h

cbuff.o:cbuff.cpp cbuff.h

pwfftcom:pwfftcom.o vsample.o fftcodec.o cbuff.o

pwfftcom.o:pwfftcom.cc


