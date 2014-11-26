# Note, FFTW version 3 must be installed to compile this code.
#
# It is assumed that the include file fftw3.h FFTW's static library files are in standard locations that GCC knows about.
# If this is not the case then you need to add the flag -I/wherever/fftw.h/is/located/at to CXXFLAGS
# and add the path of libfftw3.la to the environment variable LIBRARY_PATH.

C = gcc
CXX = g++
CXXFLAGS = -Wall -Wno-sign-compare -O3
LDFLAGS = -lm -fopenmp -lfftw3
OBJECTS = convolution_example.o fftw_wrappers.o

all: convolution_example
	
convolution_example: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ 

clean:
	rm -rf *.o convolution_example

