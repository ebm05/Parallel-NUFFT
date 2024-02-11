//
// Testing parallel speedup.
//
// E Boström 2024-01-28
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <complex>
#include "omp.h"
#include "utils.h"
#include "dirsum3d.h"
#include "gridding3d.h"
#include "defs.h"

using namespace std;

#define MAX_NUM_THREADS 12

int main(int argc, char *argv[])
{
    printf(
        "#*******************************************************************\n"
        "# Testing the speedup of the parallel algorithms.\n"
        "#*******************************************************************\n");

    double t;
    int N = 1e5;
    int M = 1000;                // Nr of uniform gridpoints, one of the three directions
    int R = 2;                   // Oversampling ratio
    int Msp = 12;                // Nr of points in spreading in each direction
    double tau = TAU(M, R, Msp); // Spreading width parameter
    double time_naive_parallel_1;
    double time_parallel_1;
    if (argc == 4)
    {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
        Msp = atoi(argv[3]);
    }
    int Mr = M * R;              // Nr of gridpoints in one direction, oversampled grid

    printf("# N=%d, M=%d, R=%d, Mr=%d, Msp=%d, tau=%f\n#\n", N, M, R, Mr, Msp, tau);
    printf("# num threads\t time naive \ttime parallel \tspeedup naive\tspeedup parallel\n"
           "#--------------------------------------------------------------------------------\n");

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
    
    complex<double> *ftau_naive_parallel = (complex<double> *)malloc(sizeof(complex<double>) * Mr * Mr * Mr);
    complex<double> *ftau_parallel = (complex<double> *)malloc(sizeof(complex<double>) * Mr * Mr * Mr);

    for (int np = 1; np <= MAX_NUM_THREADS; np++)
    {
        omp_set_num_threads(np);

        printf("%d", np);

        // Naive parallel
        t = gettime();
        fgg3d_parallel_naive(ftau_naive_parallel, f, x, y, z, N, Mr, Msp, tau);
        double time_naive_parallel = gettime() - t;
        printf("\t%f", time_naive_parallel);

        // Parallel
        t = gettime();
        fgg3d_parallel_v1(ftau_parallel, f, x, y, z, N, Mr, Msp, tau, 0, 0);
        double time_parallel = gettime() - t;
        printf("\t%f", time_parallel);

        if (np == 1)
        {
            time_naive_parallel_1 = time_naive_parallel;
            time_parallel_1 = time_parallel;
        }

        printf("\t%f", time_naive_parallel_1 / time_naive_parallel);
        printf("\t%f\n", time_parallel_1 / time_parallel);
    }
    free(x); free(y); free(z); free(f);
    free(ftau_naive_parallel); free(ftau_parallel);
    return 0;
}
