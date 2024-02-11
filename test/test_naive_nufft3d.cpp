//
// Testing the direct sum and the naive nufft (naive gridding + deconvolution) 
// implementations.
// 
// E Boström 2024-01-27
// 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <sys/time.h>
#include "omp.h"
#include "defs.h"
#include "utils.h"
#include "gridding3d.h"
#include "dirsum3d.h"
#include "deconv3d.h"

using namespace std;

int main(int argc, char *argv[])
{
  //double tot_time, sec_naive, sec_FGG, sec_par_naive_FGG, sec_par_FGG;
  int num_threads;

  int N = 1e4; // Nr of non-uniform points
  int M = 10; // Nr of uniform gridpoints, one of the three directions
  int R = 2; // Oversampling ratio
  int Mr = R*M;
  //int Nr_timing_avg = 5; // Nr of runs to do timing to compute average
  
#pragma omp parallel master
  num_threads = omp_get_num_threads();

  printf("Parameters: ");
  printf("N=%d, M=%d, R=%d\n", N, M, R);
  printf("Nr of Open MP threads used: %d\n\n", num_threads);

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

  /* Compute the direct sum */
  printf("Direct sum: \n");
  complex<double> *F_ds = (complex<double> *)malloc(sizeof(complex<double>) * M * M * M);
  //complex<double> *F_ds_test = (complex<double> *)malloc(sizeof(complex<double>) * M * M * M);
  direct3d(F_ds, f, x, y, z, N, M);
  //dirft3d1(N, x, y, z, f, -1, M, M, M, F_ds_test);
  //printf("\trelative error=%e\t(Compared to FINUFFT sum.)\n\n",norm2_rel(M*M*M, F_ds_test, F_ds));

  /* Compute naive gridding of f */
  complex<double> *ftau_naive_msp6 = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  complex<double> *ftau_naive_msp12 = (complex<double> *)malloc(sizeof(complex<double>) * Mr*Mr*Mr);
  naive_gridding(ftau_naive_msp6, f, x, y, z, N, Mr, 6, TAU(M,R,6));
  naive_gridding(ftau_naive_msp12, f, x, y, z, N, Mr, 12, TAU(M,R,12));

  /* Compute deconvolution */
  complex<double> *F_naive_msp6 = (complex<double> *)malloc(sizeof(complex<double>) * M * M * M);
  deconv3d(F_naive_msp6, ftau_naive_msp6, M, Mr, TAU(M,R,6),0);
  complex<double> *F_naive_msp12 = (complex<double> *)malloc(sizeof(complex<double>) * M * M * M);
  deconv3d(F_naive_msp12, ftau_naive_msp12, M, Mr, TAU(M,R,12),0);

  /* Compare the naive nufft to the direct sum */
  printf("Naive nufft vs direct sum:\n");
  printf("\trelative error=%e\t(Msp=6)\n",norm2_rel(M*M*M, F_ds, F_naive_msp6));
  printf("\trelative error=%e\t(Msp=12)\n",norm2_rel(M*M*M, F_ds, F_naive_msp12));


  free(x); free(y); free(z); free(f);
  free(F_ds);
  free(ftau_naive_msp6); free(ftau_naive_msp12);
  free(F_naive_msp6); free(F_naive_msp12);
  return 0;
}