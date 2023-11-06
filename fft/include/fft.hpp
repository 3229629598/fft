#ifndef FFT_HPP
#define FFT_HPP

#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "windows.h"
#include <psapi.h>
#define _USE_MATH_DEFINES
#include "math.h"

class FFT
{
    public:

    FFT(uint16_t n); //input parameter: number of terms in fft
    ~FFT();

    typedef struct
    {
        double real;
        double imag;
    }complex;

    uint16_t fft_n;
    double* fft_data;
    complex* fft_result;

    void fft(double* fft_data, complex* fft_result);
    void fft_update(double update_data);

    private:

    uint16_t storey_sum;
    uint16_t double_size;
    uint16_t* fft_tag;
    double** w_array;
    complex a, b, w;

    void rader_array(uint16_t* array);
    void trigonometric_matrix(double** w_array);
};

#endif
