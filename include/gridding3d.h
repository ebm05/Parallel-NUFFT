#ifndef GRIDDING3D_H
#define GRIDDING3D_H

#include <complex>
using namespace std;

void FGG_expansion(double*, double, int, int, double);
int positive_mod(int, int);
void naive_gridding(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double);
void fgg3d_naive(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double);
void fgg3d(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double);
void fgg3d_v2(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double);
void fgg3d_parallel_naive(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double);
void fgg3d_parallel_v1(complex<double>*, complex<double>*, double*, double*, double*, int, int, int, double, bool, int);

#endif // GRIDDING3D_H
