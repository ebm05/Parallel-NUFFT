// Gridding routines
////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include <sys/time.h>
#include "../include/defs.h"
#include "sorting.cpp"
#include <complex>
#include <vector>

using namespace std;

void FGG_expansion(vector<double> E, double x, int Mr, int Msp, double tau)
/*
 * Fast way to compute the exponentials for the fast Gaussian gridding by
 * Lindbo and Tornberg 2011.
 * The algorithm only does one exp, one pow and 2*Msp multiplications
 * instead of multiple pow operations.
 */
{
    double temp = exp(x * PI / ((double)Mr * tau));
    E[0] = pow(temp, -(double)Msp + 1);
    for (int i = 1; i < 2 * Msp; i++)
        E[i] = E[i - 1] * temp;
}

int positive_mod(int i, int n)
/*
 * Positive modulo operation used for periodic indices. Similar to the mod
 * function in MATLAB.
 */
{
    return (i % n + n) % n;
}

void naive_gridding(complex<double> *ftau, complex<double> *f, double *x, double *y, double *z, int N, int Mr, int Msp, double tau)
{

    double a = 2 * PI / (double)Mr;
    double ainv = 1 / a;
    double b = 1 / (4 * tau);

    for (int n = 0; n < N; n++)
    {
        /* Define off-grid coordinates */
        double xn = x[n];
        double yn = y[n];
        double zn = z[n];

        /* Find nearest grid point below */
        int m1 = floor(ainv * xn);
        int m2 = floor(ainv * yn);
        int m3 = floor(ainv * zn);
        double xi1 = a * m1;
        double xi2 = a * m2;
        double xi3 = a * m3;

        /* Define distances */
        double dx = xn - xi1;
        double dy = yn - xi2;
        double dz = zn - xi3;

        /* Spreading around (xi1,xi2,xi3) */
        for (int ix = -Msp + 1; ix <= Msp; ix++)
        {
            for (int iy = -Msp + 1; iy <= Msp; iy++)
            {
                for (int iz = -Msp + 1; iz <= Msp; iz++)
                {
                    int jx = positive_mod(m1 + ix, Mr);
                    int jy = positive_mod(m2 + iy, Mr);
                    int jz = positive_mod(m3 + iz, Mr);
                    int m = jx + Mr * jy + Mr * Mr * jz;
                    ftau[m] += f[n] * exp(-b * (pow(dx - a * ix, 2) + pow(dy - a * iy, 2) + pow(dz - a * iz, 2)));
                }
            }
        }
    }
}

void fgg3d_naive(complex<double> *ftau, complex<double> *f, double *x, double *y, double *z, int N, int Mr, int Msp, double tau)
/*
 * This is a naive version of the fast Gaussian gridding. For complete description
 * see the fgg3d routine below. EB 2024-01-28.
 */
{

    double a = 2 * PI / (double)Mr;
    double ainv = 1 / a;
    double b = 1 / (4 * tau);
    int P = 2 * Msp;

    //double *E2x = (double *)malloc(P * sizeof(double));
    //double *E2y = (double *)malloc(P * sizeof(double));
    //double *E2z = (double *)malloc(P * sizeof(double));
    //double *E3 = (double *)malloc(P * sizeof(double));
    //vector<double> E2x(P);
    //vector<double> E2y(P);
    //vector<double> E2z(P);
    //vector<double> E3(P);
    double E2x[12];
    double E2y[12];
    double E2z[12];
    double E3[12];

    /* Precompute E3, which is independent on (x,y,z) */
    for (int j = -Msp + 1; j <= Msp; j++)
        E3[j + Msp - 1] = exp(-pow(a * j, 2) * b);

    for (int n = 0; n < N; n++)
    {

        /* Define off-grid coordinates */
        double xn = x[n];
        double yn = y[n];
        double zn = z[n];

        /* Find nearest grid point below */
        int m1 = floor(ainv * xn);
        int m2 = floor(ainv * yn);
        int m3 = floor(ainv * zn);
        double xi1 = a * m1;
        double xi2 = a * m2;
        double xi3 = a * m3;

        /* Define distances */
        double dx = xn - xi1;
        double dy = yn - xi2;
        double dz = zn - xi3;

        /* Precompute the remaining exponentials that depend on (xn,yn,zn) */
        double E1 = exp(-b * (dx * dx + dy * dy + dz * dz));
        FGG_expansion(E2x, dx, Mr, Msp, tau);
        FGG_expansion(E2y, dy, Mr, Msp, tau);
        FGG_expansion(E2z, dz, Mr, Msp, tau);

        /* Spreading around (xi1,xi2,xi3) */
        for (int ix = 0; ix < P; ix++)
        {
            for (int iy = 0; iy < P; iy++)
            {
                for (int iz = 0; iz < P; iz++)
                {
                    int jx = positive_mod(m1 + ix - Msp + 1, Mr);
                    int jy = positive_mod(m2 + iy - Msp + 1, Mr);
                    int jz = positive_mod(m3 + iz - Msp + 1, Mr);
                    int m = jx + Mr * jy + Mr * Mr * jz;
                    ftau[m] += f[n] * E1 * E2x[ix] * E2y[iy] * E2z[iz] * E3[ix] * E3[iy] * E3[iz];
                }
            }
        }
    }

    //free(E2x);
    //free(E2y);
    //free(E2z);
    //free(E3);
}

void fgg3d(complex<double> *ftau, complex<double> *f, double *x, double *y, double *z, int N, int Mr, int Msp, double tau)
/*
 * Fast Gaussian gridding in 3d (Greengaard & Lee 2004).
 *
 * Input:
 *  ftau :  Complex array of size Mr*Mr*Mr that stores spreaded values on the
 *          oversampled grid.
 *     f :  Complex array of size N with source values f(xj,yj,zj), j=0,...,N-1.
 *     x :  Real array with N non-uniform x-coordinates xj, j=0,...,N-1.
 *     y :  Real array with N non-uniform y-coordinates yj, j=0,...,N-1.
 *     z :  Real array with N non-uniform y-coordinates zj, j=0,...,N-1.
 *     N :  Nr of non-uniform points.
 *    Mr :  Nr of point in the oversampled uniform mesh.
 *   Msp :  Nr of points used for the Gaussian spreading in each direction.
 *          The full spread is 2*Msp.
 *   tau :  Gaussian width parameter.
 *
 * EB 2024-01-28.
 */
{

    double a = 2 * PI / (double)Mr;
    double ainv = 1 / a;
    double b = 1 / (4 * tau);
    int P = 2 * Msp;

    //double *E2x = (double *)malloc(P * sizeof(double));
    //double *E2y = (double *)malloc(P * sizeof(double));
    //double *E2z = (double *)malloc(P * sizeof(double));
    //double *E3 = (double *)malloc(P * sizeof(double));
    vector<double> E2x(P);
    vector<double> E2y(P);
    vector<double> E2z(P);
    vector<double> E3(P);    

    /* Precompute E3, which is independent on (x,y,z) */
    for (int j = -Msp + 1; j <= Msp; j++)
        E3[j + Msp - 1] = exp(-pow(a * j, 2) * b);

    for (int n = 0; n < N; n++)
    {

        /* Define off-grid coordinates */
        double xn = x[n];
        double yn = y[n];
        double zn = z[n];

        /* Find nearest grid point below */
        int m1 = floor(ainv * xn);
        int m2 = floor(ainv * yn);
        int m3 = floor(ainv * zn);
        double xi1 = a * m1;
        double xi2 = a * m2;
        double xi3 = a * m3;

        /* Define distances */
        double dx = xn - xi1;
        double dy = yn - xi2;
        double dz = zn - xi3;

        /* Precompute the remaining exponentials that depend on (xn,yn,zn) */
        double E1 = exp(-b * (dx * dx + dy * dy + dz * dz));
        FGG_expansion(E2x, dx, Mr, Msp, tau);
        FGG_expansion(E2y, dy, Mr, Msp, tau);
        FGG_expansion(E2z, dz, Mr, Msp, tau);

        /* Spreading around (xi1,xi2,xi3) */
        complex<double> V0 = f[n] * E1;
        for (int ix = 0; ix < P; ix++)
        {
            complex<double> Vx = V0 * E2x[ix] * E3[ix];
            int jx = positive_mod(m1 + ix - Msp + 1, Mr);
            for (int iy = 0; iy < P; iy++)
            {
                complex<double> Vy = Vx * E2y[iy] * E3[iy];
                int jy = positive_mod(m2 + iy - Msp + 1, Mr);
                for (int iz = 0; iz < P; iz++)
                {
                    int jz = positive_mod(m3 + iz - Msp + 1, Mr);
                    int m = jx + Mr * jy + Mr * Mr * jz;
                    ftau[m] += Vy * E2z[iz] * E3[iz];
                }
            }
        }
    }
    //free(E2x);
    //free(E2y);
    //free(E2z);
    //free(E3);
}

void fgg3d_parallel_naive(complex<double> *ftau, complex<double> *f, double *x, double *y, double *z, int N, int Mr, int Msp, double tau)
/*
 * A naive parallel version of the fgg algorithm.
 * Each thread allocates a local ftau array that is spread onto. The elements
 * of the local arrays are added to the ftau array in a critical region.
 * EB 2024-01-29.
 */
{

    double a = 2 * PI / (double)Mr;
    double ainv = 1 / a;
    double b = 1 / (4 * tau);
    int P = 2 * Msp;

    // Precompute E3 (that is independent on (x,y,z))
    //double *E3 = (double *)malloc(P * sizeof(double));
    vector<double> E3(P);    
    for (int j = -Msp + 1; j <= Msp; j++)
        E3[j + Msp - 1] = exp(-pow(a * j, 2) * b);  

#pragma omp parallel
{ 
    complex<double> *ftau_local = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
    //double *E2x = (double *)malloc(P * sizeof(double));
    //double *E2y = (double *)malloc(P * sizeof(double));
    //double *E2z = (double *)malloc(P * sizeof(double));
    vector<double> E2x(P);
    vector<double> E2y(P);
    vector<double> E2z(P);

#pragma omp for
    for (int n = 0; n < N; n++)
    {
        // Define off-grid coordinates
        double xn = x[n];
        double yn = y[n];
        double zn = z[n];

        // Find nearest grid point below
        int m1 = floor(ainv * xn);
        int m2 = floor(ainv * yn);
        int m3 = floor(ainv * zn);
        double xi1 = a * m1;
        double xi2 = a * m2;
        double xi3 = a * m3;

        // Distances
        double dx = xn - xi1;
        double dy = yn - xi2;
        double dz = zn - xi3;

        // Precompute the remaining exponentials
        double E1 = exp(-(dx * dx + dy * dy + dz * dz) / (4 * tau));

        FGG_expansion(E2x, dx, Mr, Msp, tau);
        FGG_expansion(E2y, dy, Mr, Msp, tau);
        FGG_expansion(E2z, dz, Mr, Msp, tau);

        // Spreading around (xi1,xi2,xi3)
        complex<double> V0 = f[n] * E1;
        for (int ix = 0; ix < P; ix++)
        {
            complex<double> Vx = V0 * E2x[ix] * E3[ix];
            int jx = positive_mod(m1 + ix - Msp + 1, Mr);
            for (int iy = 0; iy < P; iy++)
            {
                complex<double> Vy = Vx * E2y[iy] * E3[iy];
                int jy = positive_mod(m2 + iy - Msp + 1, Mr);
                for (int iz = 0; iz < P; iz++)
                {
                    int jz = positive_mod(m3 + iz - Msp + 1, Mr);
                    int m = jx + Mr * jy + Mr * Mr * jz;
                    ftau_local[m] += Vy * E2z[iz] * E3[iz];
                }
            }
        }
    } // End of nu point loop
    // Critical region to add the local arrays to the global one.
    // (Whithout this the speedup should be close to linear.)
# pragma omp critical
    {
    for(int m = 0; m<Mr*Mr*Mr; m++)
        ftau[m] += ftau_local[m]; 
    }
    //free(E2x);
    //free(E2y);
    //free(E2z);
    free(ftau_local);
}
//free(E3);
}

void fgg3d_parallel_v1(complex<double> *ftau, complex<double> *f, double *x, double *y, double *z, int N, int Mr, int Msp, double tau, bool debug, int version)
{

    double a = 2 * PI / (double)Mr;
    double ainv = 1 / a;
    double b = 1 / (4 * tau);
    int P = 2 * Msp;
    int ns2 = (double)Msp/2;

// Precompute E3
    //double *E3 = (double *)malloc(P * sizeof(double));
    vector<double> E3(P);
    for (int j = -Msp + 1; j <= Msp; j++)
        E3[j + Msp - 1] = exp(-pow(a * j, 2) * b);
    

    // Sort the indices of the points into bins
    int64_t *sort_indices = (int64_t *)malloc(sizeof(int64_t) * N);
    double bin_size_x = 16, bin_size_y = 4, bin_size_z = 4;
    bin_sort(sort_indices, N, x, y, z, Mr, Mr, Mr, bin_size_x, bin_size_y, bin_size_z, 0);

#pragma omp parallel
    {
    //Create non-uniform index breakpoints defining the index limits for the
    //subproblems of each thread
    int nb = omp_get_num_threads();
        vector<int64_t> brk(nb + 1);
        for (int p = 0; p <= nb; ++p)
            brk[p] = (int64_t)(0.5 + N * p / (double)nb);

#pragma omp for //schedule(dynamic,1)
        for (int isub=0; isub<nb; isub++)
        {
        //int isub = omp_get_thread_num(); // Subproblem index
        int64_t N0 = brk[isub + 1] - brk[isub]; // Nr of points in subproblem

        // Copy the location and data vectors for the nonuniform points.
        // The non-uniform points are scaled to lie in [0,N1)x[0,N2)x[0,N3).
        double *kx0 = (double *)malloc(sizeof(double) * N0);
        double *ky0 = (double *)malloc(sizeof(double) * N0);
        double *kz0 = (double *)malloc(sizeof(double) * N0);
        complex<double> *f0 = (complex<double> *)malloc((sizeof(complex<double>) * N0));
        for (int64_t j = 0; j < N0; j++)
        {
            int64_t n = sort_indices[j + brk[isub]];
            kx0[j] = FOLDRESCALE(ainv*x[n], Mr, 0);
            ky0[j] = FOLDRESCALE(ainv*y[n], Mr, 0);
            kz0[j] = FOLDRESCALE(ainv*z[n], Mr, 0);
            f0[j] = f[n];
        }

        // Compute offsets of the subdomain
        int64_t os1, os2, os3, N1, N2, N3;
        get_subgrid(os1, os2, os3, N1, N2, N3, N0, kx0, ky0, kz0, Msp);
        if (debug)
            printf("%d, Index offsets = {%ld,%ld,%ld}, "
            "\tSize = %ld*%ld*%ld, "
            "\tN0 = %ld\n",isub, os1, os2, os3, N1, N2, N3, N0);

        // Allocate memory for the local cuboid of the subproblem
        complex<double> *ftau0 = (complex<double> *)calloc(sizeof(complex<double>),N1*N2*N3);

        // Allocate memory for fast Gaussian gridding arrays
        //double *E2x = (double *)malloc(P * sizeof(double));
        //double *E2y = (double *)malloc(P * sizeof(double));
        //double *E2z = (double *)malloc(P * sizeof(double));
        vector<double> E2x(P);
        vector<double> E2y(P);
        vector<double> E2z(P);

        // Loop over the subset of non-uniform points of the subproblem
        for (int n = 0; n < N0; n++)
        {

            // Local non-uniform coordinates
            double xn = a*kx0[n];
            double yn = a*ky0[n];
            double zn = a*kz0[n];

            // Find nearest grid point below
            int m1 = floor(kx0[n]);
            int m2 = floor(ky0[n]);
            int m3 = floor(kz0[n]);
            double xi1 = a * m1;
            double xi2 = a * m2;
            double xi3 = a * m3;

            // Distances
            double dx = xn - xi1;
            double dy = yn - xi2;
            double dz = zn - xi3;

            // Precompute the remaining exponentials
            double E1 = exp(-(dx * dx + dy * dy + dz * dz) / (4 * tau));
            FGG_expansion(E2x, dx, Mr, Msp, tau);
            FGG_expansion(E2y, dy, Mr, Msp, tau);
            FGG_expansion(E2z, dz, Mr, Msp, tau);

            // Spreading over the subproblem cuboid
            complex<double> V0 = f0[n] * E1;
            for (int ix = 0; ix < P; ix++)
            {
                complex<double> Vx = V0 * E2x[ix] * E3[ix];
                int jx = m1 + ix - Msp + 1 - os1;
                for (int iy = 0; iy < P; iy++)
                {
                    complex<double> Vy = Vx * E2y[iy] * E3[iy];
                    int jy = m2 + iy - Msp + 1 - os2;
                    for (int iz = 0; iz < P; iz++)
                    {
                        int jz = m3 + iz - Msp + 1 - os3;
                        int64_t m = jx + N1*jy + N1*N2*jz;
                        ftau0[m] += Vy * E2z[iz] * E3[iz];
                    }
                }
            }
        }

        // Now add the contribution to the full array.
#pragma omp critical
        add_subproblem_result(os1,os2,os3,N1,N2,N3,Mr,ftau,ftau0);

        //free(E2x);
        //free(E2y);
        //free(E2z);
        free(f0);
        free(ftau0);
        }
        //free(E3);        
    }
}



