#ifndef __FFTWWrapper_h__
#define __FFTWWrapper_h__

// Typical use of FFTW entails the following steps:
// 1. Allocate input and output buffers.
// 2. Compute a "plan" struct that. This tells FFTW what algorithms it should use when actually computing the FFT.
// 3. Execute the FFT/IFFT operation on the buffers from step 1.
//
// This file contains two classes that wrap this process nicely.
// Currently only one-dimensional transforms of real data to complex and back are supported.
//
//

#include <vector>
#include <complex.h>
#include <fftw3.h>

// Usage: (after initializing the class)
// 1. Fill input_buffer with input containing n_real_samples double numbers
//    (note, set_input_zeropadded will copy your buffer with optional zero padding)
// 2. Run execute().
// 3. Extract output by calling get_output() or directly access output_buffer[0], ..., output_buffer[output_size-1].
//    Note that the output is composed of n_real_samples/2 + 1 complex numbers.
// 
// These 3 steps can be repeated many times.
class FFTW_R2C_1D_Executor {
public:
    FFTW_R2C_1D_Executor(int n_real_samples);
    ~FFTW_R2C_1D_Executor();
    void set_input_zeropadded(const double* buffer, int size);
    void set_input_zeropadded(const std::vector<double>& vec);
    void execute();
    std::vector<double complex> get_output();

    const int input_size;
    double* const input_buffer;

    const int output_size;
    double complex* const output_buffer;

private:
    fftw_plan plan;
};

// Usage of this class is similar to that of FFTW_R2C_1D_Executor, only the input is n_real_samples/2+1 complex samples.
class FFTW_C2R_1D_Executor {
public:
    FFTW_C2R_1D_Executor(int n_real_samples);
    ~FFTW_C2R_1D_Executor();
    void set_input(const double complex* buffer, int size);
    void set_input(const std::vector<double complex>& vec);
    void execute();
    std::vector<double> get_output();

    const int input_size;
    double complex* const input_buffer;

    const int output_size;
    double* const output_buffer;

private:
    fftw_plan plan;
};


#endif

