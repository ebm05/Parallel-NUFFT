// Deconvolution routines for the last step of NUFFT
// Utilizes the forward FFT routine from the FFTw library
// Not implemented for speed.
// EB 240209
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <algorithm>
#include <fftw3.h>
#include "defs.h"

using namespace std;

void fftshift(complex<double> *F, complex<double> *Fs, int M1, int M2, int M3)
// Shifting of the outpt of the FFT, from 0,...,M/2,-M/2,...,-1 to
// -M/2,...M/2-1.
// Similar to the Matlab routine fftshift.
// Not implemented for speed
{
    for (int i1=0; i1<M1; i1++){
        for (int i2=0; i2<M2; i2++){
            for (int i3=0; i3<M3; i3++){
                int j1 = (i1 + M1/2) % M1;
                int j2 = (i2 + M2/2) % M2;
                int j3 = (i3 + M3/2) % M3;
                int n = j1 + j2*M1 + j3*M2*M3;
                int m = i1 + i2*M1 + i3*M2*M3;
                Fs[m++] = F[n];
            }
        }
    }
}

void deconv3d(complex<double> *F, complex<double> *ftau, int M, int Mr, double tau, bool debug)
// Note that the output of the forward transform from FFTw is non-shifted; i.e., it starts 
// with the zero wavenumber in each direction as follows: 0,1,...,M/2,-M/2,...,-1. To 
// downsample correctly we need to that into account. We want our output to be 
// M/2,...,M/2-1,-M/2,...,-1, in each direction.
{
    // Local array for the Fourier transform of ftau
    complex<double> *Ftau_unshifted = (complex<double>*)malloc(sizeof(complex<double>)*Mr*Mr*Mr);

    // Compute FFT of the convolved array on the oversampled grid
    fftw_plan plan = fftw_plan_dft_3d(Mr, Mr, Mr, reinterpret_cast<fftw_complex*>(ftau), reinterpret_cast<fftw_complex*>(Ftau_unshifted), FFTW_FORWARD, FFTW_ESTIMATE);  
    fftw_execute(plan);

    // Shift output of the FFT from wavenumber ordering 0,1,...,Mr/2,-Mr/2,...,-1, to 
    // -Mr/2,...,Mr/2-1
    complex<double> *Ftau = (complex<double>*)malloc(sizeof(complex<double>)*Mr*Mr*Mr);
    fftshift(Ftau_unshifted, Ftau, Mr, Mr, Mr);

    // Precomputed stuff for the inverse of the Fourier transform of the Gaussian
    double gfact = pow(PI/tau,1.5);
    double gexp = exp(tau);
    double sf = 1/((double)Mr*(double)Mr*(double)Mr);

    // Downsample and deconvolve Ftau with the Fourier transform of the Gaussian
    int d = Mr/2-M/2; // shift of downsampling from left
    if (debug) printf("(debug) d=%d\n",d);
    int i1 = 0;
    for (int k1=-M/2; k1<M/2; k1++)
    {
        int i2 = 0;
        for (int k2=-M/2; k2<M/2; k2++) 
        {
            int i3 = 0;
            for (int k3=-M/2; k3<M/2; k3++) 
            {
                int s1 = i1+d;
                int s2 = i2+d;
                int s3 = i3+d;                
                int n = s1 + s2*Mr + s3*Mr*Mr; // Index of the point on the oversampled grid
                int m = i1 + i2*M + i3*M*M;
                F[m] = gfact*pow(gexp,k1*k1+k2*k2+k3*k3)*Ftau[n]*sf;
                //if (debug) printf("k1,k2,k3,j1,j2,j3=%d,%d,%d,%d,%d,%d\n",k1,k2,k3,j1,j2,j3);
                //if (debug) printf("%f\n",real(Ftau[n]));
                i3++;
            }
            i2++;
        }
        i1++;
    }
    free(Ftau);
}