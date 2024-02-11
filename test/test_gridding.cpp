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
#include <complex>
#include "omp.h"
#include "utils.h"
#include "dirsum3d.h"
#include "gridding3d.h"
#include "defs.h"

using namespace std;

int main(int argc, char *argv[])
{
  printf(
    "*******************************************************************\n"
    "Testing the gridding algorithms.\n"
    "\n" 
    "The result of naive gridding (tested in test_naive_nufft3d.cpp) is\n"
    "used to verify correctness.\n"
    "*******************************************************************\n\n");
  int num_threads;
  double t; 
  int N=10000;
  int M=100;
  bool debug=0;
  int R = 2;
  int Msp = 12;
  double tau = TAU(M,R,Msp);
  if (argc == 4){
    N = atoi(argv[1]);
    M = atoi(argv[2]);
    Msp = atoi(argv[3]);
  }
  int Mr = M*R; // Nr of gridpoints in one direction, oversampled grid
  
#pragma omp parallel master
  num_threads = omp_get_num_threads();

  printf("Parameters: ");
  printf("N=%d, M=%d, R=%d, Mr=%d, Msp=%d, tau=%f\n", N, M, R, Mr, Msp, tau);
  printf("Nr of Open MP threads used: %d\n", num_threads);

  /* Non-unform points in [0,2*pi) and source values in (-1,1) */
  double *x = (double *)malloc(sizeof(double) * N);
  double *y = (double *)malloc(sizeof(double) * N);
  double *z = (double *)malloc(sizeof(double) * N);
  complex<double> *f = (complex<double> *)malloc(sizeof(complex<double>) * N);
  for (int j = 0; j < N; ++j)
  {
    x[j] = 2*PI * rand01();
    y[j] = 2*PI * rand01();
    z[j] = 2*PI * rand01();
    f[j] = PI*randm11()+1i*PI*randm11();
  }

  /* Naive gridding */
  complex<double> *ftau_naive = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  printf("\nNaive gridding:\n");
  t = gettime();
  naive_gridding(ftau_naive, f, x, y, z, N, Mr, Msp, tau);
  double time_naive_griddning = gettime()-t;
  printf("\tTime: %f seconds\n", time_naive_griddning);
  free(ftau_naive);

  /* Naive Fast Gaussian gridding */
  complex<double> *ftau_fgg_naive = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  printf("\nNaive fast Gaussian gridding:\n");
  t = gettime();
  fgg3d_naive(ftau_fgg_naive, f, x, y, z, N, Mr, Msp, tau);
  double time_fgg_naive = gettime()-t;
  printf("\tTime: %f seconds\n", time_fgg_naive);
  printf("\tRelative error vs naive gridding = %e\n",norm2_rel(Mr*Mr*Mr, ftau_naive, ftau_fgg_naive));
  free(ftau_fgg_naive);

  /* Fast Gaussian gridding*/
  complex<double> *ftau_fgg = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  printf("\nFast Gaussian gridding:\n");
  t = gettime();
  fgg3d_v2(ftau_fgg, f, x, y, z, N, Mr, Msp, tau);
  double time_fgg = gettime()-t;
  printf("\tTime: %f seconds\n", time_fgg);
  free(ftau_fgg);

  /* Naive parallel Fast Gaussian gridding*/
  complex<double> *ftau_fgg_naive_parallel = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  printf("\nNaive parallel Fast Gaussian gridding:\n");
  t = gettime();
  fgg3d_parallel_naive(ftau_fgg_naive_parallel, f, x, y, z, N, Mr, Msp, tau);
  double time_fgg_naive_parallel = gettime()-t;
  printf("\tTime: %f seconds\n", time_fgg_naive_parallel);
  free(ftau_fgg_naive_parallel); 

  /* Parallel Fast Gaussian gridding v1*/
  complex<double> *ftau_fgg_parallel_v1 = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  printf("\nParallel Fast Gaussian gridding v1:\n");
  t = gettime();
  fgg3d_parallel_v1(ftau_fgg_parallel_v1, f, x, y, z, N, Mr, Msp, tau, debug,0);
  double time_fgg_parallel_v1 = gettime()-t;
  printf("\tTime: %f seconds\n", time_fgg_parallel_v1);
  free(ftau_fgg_parallel_v1);
  
  free(f); free(x); free(y); free(z);

  return 0;
}