#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"

using namespace std;

FFTW_R2C_1D_Executor::FFTW_R2C_1D_Executor(int n_real_samples) :
    input_size(n_real_samples),
    input_buffer(fftw_alloc_real(n_real_samples)),
    output_size(n_real_samples/2 + 1),
    output_buffer(fftw_alloc_complex(n_real_samples/2 + 1))
{
    plan = fftw_plan_dft_r2c_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

FFTW_R2C_1D_Executor::~FFTW_R2C_1D_Executor()
{
    fftw_destroy_plan(plan);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const double* buffer, int size)
{
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(double)*size);
    memset(&input_buffer[size], 0, sizeof(double)*(input_size - size));
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const vector<double>& vec)
{
    set_input_zeropadded(&vec[0], vec.size());
}

void FFTW_R2C_1D_Executor::execute()
{
    fftw_execute(plan);
}

vector<double complex> FFTW_R2C_1D_Executor::get_output()
{
    return vector<double complex>(output_buffer, output_buffer + output_size);
}

FFTW_C2R_1D_Executor::FFTW_C2R_1D_Executor(int n_real_samples) : 
    input_size(n_real_samples/2 + 1),
    input_buffer(fftw_alloc_complex(n_real_samples/2 + 1)),
    output_size(n_real_samples),
    output_buffer(fftw_alloc_real(n_real_samples))
{
    plan = fftw_plan_dft_c2r_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

FFTW_C2R_1D_Executor::~FFTW_C2R_1D_Executor()
{
    fftw_destroy_plan(plan);
    fftw_free(input_buffer);
    fftw_free(output_buffer);
}

void FFTW_C2R_1D_Executor::set_input(const double complex* buffer, int size)
{
    assert(size == input_size);
    memcpy(input_buffer, buffer, sizeof(double complex)*size);
    memset(&input_buffer[size], 0, sizeof(double complex)*(input_size - size));
}

void FFTW_C2R_1D_Executor::set_input(const vector<double complex>& vec)
{
    set_input(&vec[0], vec.size());
}

void FFTW_C2R_1D_Executor::execute()
{
    fftw_execute(plan);
}

vector<double> FFTW_C2R_1D_Executor::get_output()
{
    return vector<double>(output_buffer, output_buffer + output_size);
}
