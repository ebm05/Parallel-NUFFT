#ifndef UTILS_H
#define UTILS_H

#include <complex>
using namespace std;

double gettime();
double rand1();
double periodic_abs(double, double, double);
double norm2_rel(int64_t, complex<double>*, complex<double>*);
double norm2_dir(int64_t, complex<double>*, complex<double>*);
double norm2(int64_t, complex<double>*);
double norminf(int64_t, complex<double>*);
void print_array(complex<double>*, int);

#endif // End of the #ifndef UTILS_H
