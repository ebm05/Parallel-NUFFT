//
// Main testing program for the different gridding algorithms.
//
// The results are compared to the naive gridding result. The naive nufft (gridding + deconvolution)
// is tested in test_naivenufft3d.cpp.
//
// E Boström 2024-01-28
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "omp.h"
#include "utils.h"
#include "dirsum3d.h"
#include "gridding3d.h"
#include "defs.h"

using namespace std;

int main(int argc, char *argv[])
{
    printf(
        "#*******************************************************************\n"
        "# Testing the performance of the serial FGG algorithm.\n"
        "#*******************************************************************\n");

    double t;
    int M = 100;       // Nr of uniform gridpoints, one of the three directions
    int R = 2;                   // Oversampling ratio
    int Mr = M * R;              // Nr of gridpoints in one direction, oversampled grid
    int Msp = 12;    // Nr of points in spreading in each direction
    double tau = TAU(M, R, Msp); // Spreading width parameter
    int niter = 16;

    if (argc == 4){
        M = atoi(argv[1]);
        Msp = atoi(argv[2]);
        niter = atoi(argv[3]);
    }

    printf("# M=%d, R=%d, Mr=%d, Msp=%d, tau=%f\n#\n", M,R,Mr,Msp,tau);

    printf("# N\tnaive  \tnaive FGG\tFGG  \tSpeedup ng/FGG\tSpeedup nFGG/FGG\n"
           "#--------------------------------------------------------------------------------\n");

    for (int j = 1; j <= niter; j++)
    {
        int N = pow(2, j);
        printf("%d",N);

        /* Non-unform points in [0,2*pi) and source values in (-1,1) */
        double *x = (double *)malloc(sizeof(double) * N);
        double *y = (double *)malloc(sizeof(double) * N);
        double *z = (double *)malloc(sizeof(double) * N);
        complex<double> *f = (complex<double> *)malloc(sizeof(complex<double>) * N);
        for (int j = 0; j < N; ++j)
        {
            x[j] = 2 * PI * rand01();
            y[j] = 2 * PI * rand01();
            z[j] = 2 * PI * rand01();
            f[j] = PI * randm11() + 1i * PI * randm11();
        }

        // Naive gridding
        complex<double> *ftau_naive = (complex<double> *)malloc(sizeof(complex<double>) * Mr * Mr * Mr);
        t = gettime();
        naive_gridding(ftau_naive, f, x, y, z, N, Mr, Msp, tau);
        double time_naive_griddning = gettime() - t;
        printf("\t%f", time_naive_griddning);

        // Naive Fast Gaussian gridding
        complex<double> *ftau_fgg_naive = (complex<double> *)malloc(sizeof(complex<double>) * Mr * Mr * Mr);
        t = gettime();
        fgg3d_naive(ftau_fgg_naive, f, x, y, z, N, Mr, Msp, tau);
        double time_fgg_naive = gettime() - t;
        printf("\t%f", time_fgg_naive);

        // Fast Gaussian gridding
        complex<double> *ftau_fgg = (complex<double> *)malloc(sizeof(complex<double>) * Mr * Mr * Mr);
        t = gettime();
        fgg3d_v2(ftau_fgg, f, x, y, z, N, Mr, Msp, tau);
        double time_fgg = gettime() - t;
        printf("\t%f", time_fgg);

        printf("\t%f", time_fgg_naive/time_fgg);
        printf("\t%f\n", time_naive_griddning/time_fgg);
        free(x); free(y); free(z);
        free(ftau_naive); free(ftau_fgg_naive); free(ftau_fgg);
    }
    return 0;
}