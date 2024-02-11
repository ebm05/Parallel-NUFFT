#ifndef SORTING_H
#define SORTING_H

#include <complex>
using namespace std;

void arrayrange(int64_t, double*, double*, double*);
void bin_sort(int64_t*, int64_t, double*, double*, double*,int64_t, int64_t, int64_t, double, double, double, int);
void get_subgrid(int64_t&, int64_t&, int64_t&, int64_t&, int64_t&, int64_t&, int64_t, double*, double*, double*, int);
void add_subproblem_result(int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, complex<double>*, complex<double>*);

#endif // End of the #ifndef SORTING_H