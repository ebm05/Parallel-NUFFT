#ifndef DECONV3D_H
#define DECONV3D_H

#include <complex>
using namespace std;

void fftshift(complex<double>*, complex<double>*, int, int, int);
void deconv3d(complex<double>*, complex<double>*, int, int, double, bool);

#endif // End of the #ifndef DECONV3D_H