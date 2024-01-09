#ifndef FFT_HPP
#define FFT_HPP

#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#define _USE_MATH_DEFINES
#include "math.h"

class FFT
{
    public:

    FFT(uint16_t n); //input parameter: number of terms in fft
    ~FFT();

    uint16_t fft_n;
    double* fft_input_real;
    double* fft_input_imag;
    double* fft_output_real;
    double* fft_output_imag;

    void fft(double* fft_input_real, double* fft_input_imag, double* fft_output_real, double* fft_output_imag);
    void fft_update(double update_data);
    void ifft(double* ifft_output_real, double* ifft_output_imag, double* ifft_output_sqrt);

    private:

    uint16_t storey_sum;
    uint16_t double_size;
    uint16_t* fft_tag;
    double** w_array;

    void rader_array(uint16_t* array);
    void trigonometric_matrix(double** w_array);
};

#endif
