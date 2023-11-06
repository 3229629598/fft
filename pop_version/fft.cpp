#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#define _USE_MATH_DEFINES
#include "math.h"

#define fft_n 1024    //fft_n = 2^x

typedef struct
{
    double real;
    double imag;
}complex;

double fft_data[fft_n]={0};
complex fft_result[fft_n]={0};

void rader_array(uint16_t* array, uint16_t n)
{
    uint16_t k=n/2;
    array[n-1]=n-1;
    for(uint16_t i=1; i<n-1; ++i) 
    {
        k=n/2;
        array[i]=array[i-1];
        while (array[i]>=k)
        {
            array[i]-=k;
            k/=2;
        }
        array[i]+=k;
    }
}

void trigonometric_matrix(complex*** w_array, uint16_t n, uint16_t storey_sum)
{
    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<n; pair_site+=n_half<<1)
        {
            for(uint16_t k=0; k<n_half; ++k)
            {
                double w_ang=-M_PI*k/n_half;
                w_array[storey_cnt][pair_site][k].real=cos(w_ang);
                w_array[storey_cnt][pair_site][k].imag=sin(w_ang);
            }
        }
    }
}

void fft(double* fft_data, uint16_t* fft_tag, complex*** w_array, complex* fft_result, uint16_t n, uint16_t storey_sum)
{
    complex a, b, w;
    for(uint16_t i=0; i<n; ++i)
    {
        fft_result[i].real=fft_data[fft_tag[i]];
    }

    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<n; pair_site+=n_half<<1)
        {
            for(uint16_t k=0; k<n_half; ++k)
            {
                w.real=w_array[storey_cnt][pair_site][k].real;
                w.imag=w_array[storey_cnt][pair_site][k].imag;
                a.real=fft_result[pair_site+k].real;
                a.imag=fft_result[pair_site+k].imag;
                b.real=fft_result[pair_site+k+n_half].real;
                b.imag=fft_result[pair_site+k+n_half].imag;

                fft_result[pair_site+k].real               =a.real + w.real * b.real - w.imag * b.imag;
                fft_result[pair_site+k].imag             =a.imag + w.real * b.imag + w.imag * b.real;
                fft_result[pair_site+k+n_half].real  =a.real - w.real * b.real + w.imag * b.imag;
                fft_result[pair_site+k+n_half].imag=a.imag - w.real * b.imag - w.imag * b.real;
            }
        }
    }

}

int main()
{
    uint16_t storey_sum=log(fft_n)/log(2);
    uint16_t double_size=sizeof(double);

    uint16_t fft_tag[fft_n]={0};

    complex*** w_array=NULL;
    w_array=(complex***)malloc(storey_sum*sizeof(complex**));
    for(int i=0; i<storey_sum; ++i)
    {
        w_array[i]=(complex**)malloc(fft_n*sizeof(complex*));
        for(int j=0; j<fft_n; ++j)
            w_array[i][j]=(complex*)malloc(fft_n*sizeof(complex));
    }

    rader_array(fft_tag, fft_n);
    trigonometric_matrix(w_array, fft_n, storey_sum);

    for(int i=0; i<fft_n; ++i)
        fft_data[i]=i+1;

    time_t start=0,stop=0;
    start=clock();

    for(int i=0; i<1000; ++i)
    {
        fft(fft_data, fft_tag, w_array, fft_result, fft_n, storey_sum);
        memcpy(fft_data,fft_data+1,(fft_n-1)*double_size);
        fft_data[fft_n-1]=fft_n+i;
    }

    // fft(fft_data, fft_tag, w_array, fft_result, fft_n, storey_sum);
    // for(int i=0;i<fft_n;i++)
    //     printf("(%.4f, %.4fi) ", fft_result[i].real, fft_result[i].imag);
    
    stop=clock();
    printf("\nTime used: %f ms",(double)(stop-start));

    return 0;
}
