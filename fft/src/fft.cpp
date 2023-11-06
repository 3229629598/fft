#include <fft.hpp>

FFT::FFT(uint16_t n)
{
    fft_n=n;
    storey_sum=log(fft_n)/log(2);
    double_size=sizeof(double);
    
    fft_data=(double*)malloc(fft_n*sizeof(double));
    memset(fft_data,0,sizeof(fft_data));
    fft_tag=(uint16_t*)malloc(fft_n*sizeof(uint16_t));
    memset(fft_tag,0,sizeof(fft_tag));
    w_array=(complex***)malloc(storey_sum*sizeof(complex**));
    for(int i=0; i<storey_sum; ++i)
    {
        w_array[i]=(complex**)malloc(fft_n*sizeof(complex*));
        for(int j=0; j<fft_n; ++j)
        {
            w_array[i][j]=(complex*)malloc(fft_n*sizeof(complex));
            memset(w_array[i][j],0,sizeof(w_array[i][j]));
        }
    }
    fft_result=(complex*)malloc(fft_n*sizeof(complex));
    memset(fft_result,0,sizeof(fft_result));

    rader_array(fft_tag);
    trigonometric_matrix(w_array);
}

FFT::~FFT()
{
    free(fft_data);
	fft_data=NULL;
    free(fft_tag);
	fft_tag=NULL;
    free(w_array);
	w_array=NULL;
    free(fft_result);
	fft_result=NULL;
}

void FFT::rader_array(uint16_t* array)
{
    uint16_t k=fft_n/2;
    array[0]=0;
    array[fft_n-1]=fft_n-1;
    for(uint16_t i=1; i<fft_n-1; ++i) 
    {
        k=fft_n/2;
        array[i]=array[i-1];
        while (array[i]>=k)
        {
            array[i]-=k;
            k/=2;
        }
        array[i]+=k;
    }
}

void FFT::trigonometric_matrix(complex*** w_array)
{
    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<fft_n; pair_site+=n_half<<1)
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

void FFT::fft(double* fft_data, complex* fft_result)
{
    for(uint16_t i=0; i<fft_n; ++i)
    {
        fft_result[i].real=fft_data[fft_tag[i]];
    }

    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<fft_n; pair_site+=n_half<<1)
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

void FFT::fft_update(double update_data)
{
    memcpy(fft_data,fft_data+1,(fft_n-1)*double_size);
    fft_data[fft_n-1]=update_data;
    fft(fft_data, fft_result);
}

int main()
{
    FFT my_fft(1024);

    time_t start=0,stop=0;
    start=clock();

    for(int i=0; i<1000; ++i)
    {
        my_fft.fft_update(1024+i);
    }

    // for(int i=0; i<my_fft.fft_n; ++i)
    //     my_fft.fft_data[i]=i+1;
    // my_fft.fft(my_fft.fft_data, my_fft.fft_result);
    // for(int i=0;i<my_fft.fft_n;i++)
    //     printf("(%.4f, %.4fi) ", my_fft.fft_result[i].real, my_fft.fft_result[i].imag);
    
    stop=clock();
    printf("\nTime used: %f ms",(double)(stop-start));

    return 0;
}