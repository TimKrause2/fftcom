CC=g++
CFLAGS=-pthread -ggdb
CPPFLAGS=-ggdb -O2 `pkg-config --cflags freetype2` -I esutil -pthread -fPIC
LDLIBS=-ljack -lrt -lfftw3_threads -lfftw3 -lX11 -lGL -lGLU -lGLEW -lm -L esutil -lutil `pkg-config --libs freetype2`
LDFLAGS=-ggdb -pthread


fftcom:fftcom.o vsample.o fftcodec.o font.o

fftcom.o:fftcom.cc

vsample.o:vsample.cc vsample.h

fftcodec.o:fftcodec.cc fftcodec.h vsample.h

font.o:font.cc

ptest:ptest.cc

