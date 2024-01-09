#include <fft.hpp>

FFT::FFT(uint16_t n)
{
    storey_sum=log(n)/log(2);
    fft_n=1<<storey_sum;
    double_size=sizeof(double);
    
    do
    {
        fft_tag=(uint16_t*)malloc(fft_n*sizeof(uint16_t));
        memset(fft_tag,0,fft_n*sizeof(uint16_t));
        w_array=(double**)malloc(storey_sum*sizeof(double*));
        for(int i=0; i<storey_sum; ++i)
        {
            w_array[i]=(double*)malloc(fft_n*sizeof(double));
            memset(w_array[i],0,fft_n*sizeof(double));
        }
        fft_input_real=(double*)malloc(fft_n*sizeof(double));
        memset(fft_input_real,0,fft_n*sizeof(double));
        fft_input_imag=(double*)malloc(fft_n*sizeof(double));
        memset(fft_input_imag,0,fft_n*sizeof(double));
        fft_output_real=(double*)malloc(fft_n*sizeof(double));
        memset(fft_output_real,0,fft_n*sizeof(double));
        fft_output_imag=(double*)malloc(fft_n*sizeof(double));
        memset(fft_output_imag,0,fft_n*sizeof(double));
    } while (fft_tag==NULL||w_array==NULL||fft_input_real==NULL||fft_input_imag==NULL||fft_output_real==NULL||fft_output_imag==NULL);
    
    rader_array(fft_tag);
    trigonometric_matrix(w_array);
}

FFT::~FFT()
{
    free(fft_tag);
	fft_tag=NULL;
    free(w_array);
	w_array=NULL;
    free(fft_input_real);
	fft_input_real=NULL;
    free(fft_input_imag);
    fft_input_imag=NULL;
    free(fft_output_real);
	fft_output_real=NULL;
    free(fft_output_imag);
    fft_output_imag=NULL;
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

void FFT::trigonometric_matrix(double** w_array)
{
    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<fft_n; pair_site+=n_half<<1)
        {
            for(uint16_t k=0; k<n_half; ++k)
            {
                double w_ang=-M_PI*k/n_half;
                w_array[storey_cnt][pair_site+k]=cos(w_ang);
                w_array[storey_cnt][pair_site+k+n_half]=sin(w_ang);
            }
        }
    }
}

void FFT::fft(double* fft_input_real, double* fft_input_imag, double* fft_output_real, double* fft_output_imag)
{
    for(uint16_t i=0; i<fft_n; ++i)
    {
        fft_output_real[i]=fft_input_real[fft_tag[i]];
        fft_output_imag[i]=fft_input_imag[fft_tag[i]];
    }

    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<fft_n; pair_site+=n_half<<1)
        {
            for(uint16_t k=0; k<n_half; ++k)
            {
                w1=w_array[storey_cnt][pair_site+k];
                w2=w_array[storey_cnt][pair_site+k+n_half];
                a1=fft_output_real[pair_site+k];
                a2=fft_output_imag[pair_site+k];
                b1=fft_output_real[pair_site+k+n_half];
                b2=fft_output_imag[pair_site+k+n_half];

                double w1xb1=w1*b1;
                double w1xb2=w1*b2;
                double w2xb1=w2*b1;
                double w2xb2=w2*b2;

                fft_output_real[pair_site+k]               =a1 + w1xb1 - w2xb2;
                fft_output_imag[pair_site+k]             =a2 + w1xb2 + w2xb1;
                fft_output_real[pair_site+k+n_half]  =a1 - w1xb1 + w2xb2;
                fft_output_imag[pair_site+k+n_half]=a2 - w1xb2 - w2xb1;
            }
        }
    }
}

void FFT::fft_update(double update_data)
{
    memcpy(fft_input_real,fft_input_real+1,(fft_n-1)*double_size);
    fft_input_real[fft_n-1]=update_data;
    fft(fft_input_real,fft_input_imag,fft_output_real,fft_output_imag);
}

void FFT::ifft(double* ifft_output_real, double* ifft_output_imag, double* ifft_output_sqrt)
{
    for(uint16_t i=0; i<fft_n; ++i)
    {
        ifft_output_real[i]=fft_output_real[fft_tag[i]];
        ifft_output_imag[i]=-fft_output_imag[fft_tag[i]];
    }

    for(uint16_t storey_cnt=0, n_half=0; storey_cnt<storey_sum; ++storey_cnt)
    {
        n_half=1<<storey_cnt;
        for(uint16_t pair_site=0; pair_site<fft_n; pair_site+=n_half<<1)
        {
            for(uint16_t k=0; k<n_half; ++k)
            {
                w1=w_array[storey_cnt][pair_site+k];
                w2=w_array[storey_cnt][pair_site+k+n_half];
                a1=ifft_output_real[pair_site+k];
                a2=ifft_output_imag[pair_site+k];
                b1=ifft_output_real[pair_site+k+n_half];
                b2=ifft_output_imag[pair_site+k+n_half];

                double w1xb1=w1*b1;
                double w1xb2=w1*b2;
                double w2xb1=w2*b1;
                double w2xb2=w2*b2;

                ifft_output_real[pair_site+k]               =a1 + w1xb1 - w2xb2;
                ifft_output_imag[pair_site+k]             =a2 + w1xb2 + w2xb1;
                ifft_output_real[pair_site+k+n_half]  =a1 - w1xb1 + w2xb2;
                ifft_output_imag[pair_site+k+n_half]=a2 - w1xb2 - w2xb1;
            }
        }
    }

    for(uint16_t i=0; i<fft_n; ++i)
    {
        ifft_output_real[i]=ifft_output_real[i]/fft_n;
        ifft_output_imag[i]=-ifft_output_imag[i]/fft_n;
        ifft_output_sqrt[i]=sqrt(ifft_output_real[i]*ifft_output_real[i]+ifft_output_imag[i]*ifft_output_imag[i]);
    }
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
    //     my_fft.fft_input_real[i]=i+1;
    // my_fft.fft(my_fft.fft_input_real,my_fft.fft_input_imag,my_fft.fft_output_real,my_fft.fft_output_imag);
    // printf("\nfft_data:\n");
    // for(int i=0;i<my_fft.fft_n;i++)
    //     printf("(%.4fr, %.4fi) ", my_fft.fft_output_real[i], my_fft.fft_output_imag[i]);
    
    // double ifft_real[8]={0}, ifft_imag[8]={0}, ifft_sqrt[8]={0};
    // my_fft.ifft(ifft_real,ifft_imag,ifft_sqrt);
    // printf("\nifft_data:\n");
    // for(int i=0;i<my_fft.fft_n;i++)
    //     printf("(%.4fr, %.4fi, %.4fs) ", ifft_real[i], ifft_imag[i], ifft_sqrt[i]);
    
    stop=clock();
    printf("\nTime used: %f ms\n",(double)(stop-start));
    
    return 0;
}