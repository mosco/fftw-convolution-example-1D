#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"

using namespace std;

void print_double_array(const double* arr, int length)
{
    for (int i = 0; i < length; ++i) {
        cout << arr[i] << ", ";
    }
    cout << endl;
}
void print_complex_array(const double complex* arr, int length)
{
    for (int i = 0; i < length; ++i) {
        cout << creal(arr[i]) << "+" << cimag(arr[i]) << "i, ";
    }
    cout << endl;
}

template<class T>
void print_vector(const vector<T>& vec)
{
    for (unsigned int i = 0; i < vec.size(); ++i) {
        cout << vec[i] << ", ";
    }
    cout << endl;
}

// This function computes the discrete convolution of two arrays:
// result[i] = a[i]*b[0] + a[i-1]*b[1] + ... + a[0]*b[i]
// a and b can be vectors of different lengths, this function is careful to never
// exceed the bounds of the vectors.
vector<double> convolve(const vector<double>& a, const vector<double>& b)
{
    int n_a = a.size();
    int n_b = b.size();
    vector<double> result(n_a + n_b - 1);

    for (int i = 0; i < n_a + n_b - 1; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i; ++j) {
            sum += ((j < n_a) && (i-j < n_b)) ? a[j]*b[i-j] : 0.0;
        }
        result[i] = sum;
    }
    return result;
}

template <class T>
vector<T> vector_elementwise_multiply(const vector<T> a, const vector<T> b)
{
    assert(a.size() == b.size());
    vector<T> result(a.size());
    for (int i = 0; i < result.size(); ++i) {
        result[i] = a[i]*b[i];
    }
    return result;
}

// Convolution of real vectors using the Fast Fourier Transform and the convolution theorem.
// See http://en.wikipedia.org/w/index.php?title=Convolution&oldid=630841165#Fast_convolution_algorithms
vector<double> fftw_convolve(vector<double>& a, vector<double>& b)
{
    // Recall that element-wise
    int padded_length = a.size() + b.size() - 1;
    
    // Compute Fourier transform of vector a
    
    FFTW_R2C_1D_Executor fft_a(padded_length);
    fft_a.set_input_zeropadded(a);

    cout << "a: ";
    print_double_array(fft_a.input_buffer, fft_a.input_size);

    fft_a.execute();

    cout << "FFT(a): ";
    print_complex_array(fft_a.output_buffer, fft_a.output_size);
    cout << endl;

    // Compute Fourier transform of vector b
    
    FFTW_R2C_1D_Executor fft_b(padded_length);
    fft_b.set_input_zeropadded(b);

    cout << "b: ";
    print_double_array(fft_b.input_buffer, fft_b.input_size);

    fft_b.execute();

    cout << "FFT(b): ";
    print_complex_array(fft_b.output_buffer, fft_b.output_size);
    cout << endl;

    // Perform element-wise product of FFT(a) and FFT(b)
    // then compute inverse fourier transform.
    FFTW_C2R_1D_Executor ifft(padded_length);
    assert (ifft.input_size == fft_a.output_size);
    ifft.set_input(vector_elementwise_multiply(fft_a.get_output(), fft_b.get_output()));

    ifft.execute();

    // FFTW returns unnormalized output. To normalize it one must divide each element
    // of the result by the number of elements.
    assert(ifft.output_size == padded_length);
    vector<double> result = ifft.get_output();
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] /= padded_length;
    }

    return result;
}

int main()
{
    vector<double> a;
    a.push_back(2);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    cout << "First vector (a): ";
    print_vector(a);

    vector<double> b;
    b.push_back(1);
    b.push_back(0);
    b.push_back(7);
    cout << "Second vector (b): ";
    print_vector(b);

    cout << "==== Naive convolution ===========================================\n";

    vector<double> result_naive = convolve(a, b);
    cout << "Naive convolution result:\n";
    print_vector(result_naive);

    cout << "==== FFT convolution =============================================\n";

    vector<double> result_fft = fftw_convolve(a, b);
    cout << "FFT convolution result:\n";
    print_vector(result_fft);
}

