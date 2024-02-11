// Direct sum of the FFT
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include "defs.h"
#include "dirsum3d.h"

using namespace std;

void direct3d(complex<double> *F, complex<double> *f, double *x, double *y, double *z, int N, int M)
//
// Direct sum that computes the 3d Fourier coefficients.
// It is mainly used to validate the nufft implementations.
// 
// Input:
//  F : Complex array of size M*M*M that stores the Fourier coefficients;
//      F(k1,k2,k3) for k1,k2,k3 in -M/2,...,M/2-1 (shifted).
//  f : Complex array of size N with source values.
//  x : Real array with N non-uniform x-coordinates.
//  y : Real array with N non-uniform y-coordinates. 
//  z : Real array with N non-uniform y-coordinates.
//  N : Nr of non-uniform points.
//  M : Nr of Fourier modes in each direction (total number of modes is M*M*M)
//  
// EB 2024-01-27.
//
{

    for (int m = 0; m < M * M * M; ++m)
        F[m] = complex<double>(0, 0);

    for (int k1 = -M / 2; k1 < M / 2; ++k1)
    {
        for (int k2 = -M / 2; k2 < M / 2; ++k2)
        {
            for (int k3 = -M / 2; k3 < M / 2; ++k3)
            {
                for (int n = 0; n < N; ++n)
                {
                    int i1 = k1 + M / 2;
                    int i2 = k2 + M / 2;
                    int i3 = k3 + M / 2;
                    int i = i1 + M * i2 + M * M * i3;

                    F[i] += f[n] * exp( -1i * ((double)k1*x[n] + (double)k2*y[n] + (double)k3*z[n]));
                }
            }
        }
    }
}
