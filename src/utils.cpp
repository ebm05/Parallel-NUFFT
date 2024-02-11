// Utility functions
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <complex>
#include "defs.h"
#include "utils.h"

using namespace std;

double gettime(void)
// Timing function
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

double rand1(){
// Random value in the interval (0,1).
    return (double)rand()/(double)RAND_MAX;
}

double periodic_abs(double x1, double x2, double L) {
    double dx;
    dx = fabs(x1 - x2);
    if (dx > L/2) {
        dx = L - dx;
    }
    return dx;
}

double norm2_rel(int64_t n, complex<double>* a, complex<double>* b)
// ||a-b||_2 / ||a||_2
{
  double err = 0.0, nrm = 0.0;
  for (int64_t m=0; m<n; ++m) {
    nrm += real(conj(a[m])*a[m]);
    complex<double> diff = a[m]-b[m];
    err += real(conj(diff)*diff);
  }
  return sqrt(err/nrm);
}
double norm2_dir(int64_t n, complex<double>* a, complex<double>* b)
// ||a-b||_2
{
  double err = 0.0;   // compute error 2-norm
  for (int64_t m=0; m<n; ++m) {
    complex<double> diff = a[m]-b[m];
    err += real(conj(diff)*diff);
  }
  return sqrt(err);
}
double norm2(int64_t n, complex<double>* a)
// ||a||_2
{
  double nrm = 0.0;
  for (int64_t m=0; m<n; ++m)
    nrm += real(conj(a[m])*a[m]);
  return sqrt(nrm);
}
double norminf(int64_t n, complex<double>* a)
// ||a||_infty
{
  double nrm = 0.0;
  for (int64_t m=0; m<n; ++m) {
    double aa = real(conj(a[m])*a[m]);
    if (aa>nrm) nrm = aa;
  }
  return sqrt(nrm);
}

void print_array(complex<double> *arr, int N)
{
  for (int n=0; n<N; n++){
    printf("arr[%d]=%f+i(%f)\n",n,real(arr[n]),imag(arr[n]));
  }
}
