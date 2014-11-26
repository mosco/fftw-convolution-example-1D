fftw-convolution-example-1D
===========================

A simple C++ example of performing a one-dimensional discrete convolution of real vectors using the Fast Fourier Transform (FFT) as implemented in the [FFTW 3](http://www.fftw.org/) library.

This code is a simple and direct application of the well-known Convolution Theorem. It is not efficient, but meant to be easy to understand.

# Description of the algorithm
Let v1, v2 be two vectors of real numbers. The discrete (linear) convolution of these vectors can be computed by the following procedure:
 1. Zero pad both vectors to length size(v1)+size(v2)-1.
 2. Compute the discrete Fourier transforms of the padded vectors.
 3. Calculate the pointwise product of these fourier transforms. i.e. result[i] = Fourier(pad(v1))[i] * Fourier(pad(v2))[i]
 4. Transform the result by the inverse Fourier transform.

The reason why one must pad the vectors beforehand is that without padding, this procedure would compute the *cyclic* convolution of v1 and v2. Since we are interested in the linear convolution we need to add enough padding so that the wraparound will not mix different parts of the linear convolution.

For more information, see [relevant article on Wikipedia](http://en.wikipedia.org/w/index.php?title=Convolution&oldid=630841165#Circular_discrete_convolution).

# Implementation notes

* The files fftw_wrappers.hh/cc contain two classes that wrap some of the unpleasantries involved in using FFTW.
* convolution_example.cc contains a naive implementation of a linear convolution with runtime O(n^2) and the FFT-based convolution with runtime O(n^2 log n). main() loads two small vectors and prints out the (hopefully identical) results of running both convolution functions.
* The code is terribly inefficient, since it performs a lot of unnecessary memory allocation and computes an FFTW "plan" which it only uses once.

# Compilation notes

To compile the code, install FFTW version 3 on your computer and then run 'make'.

It is assumed that the include file fftw3.h and FFTW's static library files are in standard locations that GCC knows about. If this is not the case then you need to add the flag -I/wherever/fftw.h/is/located/at to CXXFLAGS (defined in Makefile) and add the path of libfftw3.la to the environment variable LIBRARY_PATH.

# Contact
Feel free to ask any questions: moscovich@gmail.com.
     
