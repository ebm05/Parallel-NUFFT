// Sorting routines to be used in the gridding
//
#include <vector>
#include <complex>
#include "omp.h"
#include <stdlib.h>
#include "defs.h"
#include "sorting.h"

using namespace std;

void arrayrange(int64_t n, double *a, double *lo, double *hi)
// Modified version of routine from the FNUFFT project
// https://github.com/flatironinstitute/finufft
//
// With a a length-n array, writes out min(a) to lo and max(a) to hi,
// so that all a values lie in [lo,hi].
// If n==0, lo and hi are not finite.
//
{
    *lo = INFINITY;
    *hi = -INFINITY;
    for (int64_t m = 0; m < n; ++m)
    {
        if (a[m] < *lo)
            *lo = a[m];
        if (a[m] > *hi)
            *hi = a[m];
    }
}

void bin_sort(int64_t *ret, int64_t M, double *x, double *y, double *z,
                           int64_t N1, int64_t N2, int64_t N3,
                           double bin_size_x, double bin_size_y, double bin_size_z, int debug)
// Modified version of routine from the FNUFFT project
// https://github.com/flatironinstitute/finufft
{

    // Calculate the number of bins. 
    // (Here the +1 is needed to allow round-off error causing i1=N1/bin_size_x,
    // for kx near +pi, ie foldrescale gives N1 (exact arith would be 0 to N1-1).
    // Note that round-off near kx=-pi stably rounds negative to i1=0.)
    int64_t nbins1 = N1 / bin_size_x + 1;
    int64_t nbins2 = N2 / bin_size_y + 1;
    int64_t nbins3 = N3 / bin_size_z + 1;
    int64_t nbins = nbins1 * nbins2 * nbins3;
    double a1 = (double)N1 / (2 * PI);
    double a2 = (double)N2 / (2 * PI);
    double a3 = (double)N3 / (2 * PI);

    // Count how many points there are in each bin.
    std::vector<int64_t> counts(nbins, 0);
    for (int64_t i = 0; i < M; i++)
    {
        int64_t i1 = FOLDRESCALE(a1*x[i], N1, 0) / bin_size_x;
        int64_t i2 = FOLDRESCALE(a2*y[i], N2, 0) / bin_size_y;
        int64_t i3 = FOLDRESCALE(a3*z[i], N3, 0) / bin_size_z;
        int64_t bin = i1 + nbins1 * (i2 + nbins2 * i3);
        counts[bin]++;
    }
    if (debug) { printf("Number of points in each bin:\n");
        for (int i=0; i<nbins; i++) printf("counts[%d]=%ld\n", i, counts[i]);}

    // Overwrite the counts vector with the offsets
    int64_t current_offset = 0;
    for (int64_t i = 0; i < nbins; i++)
    {
        int64_t tmp = counts[i];
        counts[i] = current_offset;
        current_offset += tmp;
    }
    if (debug) { printf("Offsets:\n");
        for (int i=0; i<nbins; i++) printf("counts[%d]=%ld\n", i, counts[i]);}

    // The ret array is now updated with the indices.
    for (int64_t i = 0; i < M; i++)
    {
        // Find the bin indices
        int64_t i1 = FOLDRESCALE(a1*x[i], N1, 0) / bin_size_x;
        int64_t i2 = FOLDRESCALE(a2*y[i], N2, 0) / bin_size_y;
        int64_t i3 = FOLDRESCALE(a3*z[i], N3, 0) / bin_size_z;
        int64_t bin = i1 + nbins1 * (i2 + nbins2 * i3);
        
        // Update indices. 
        // (Note that the counts array is incremented in each iteration. Therefore moves 
        // the index one step forward inside the bin.)
        ret[counts[bin]] = i;
        ++counts[bin];
    }
}

void get_subgrid(int64_t &offset1, int64_t &offset2, int64_t &offset3, int64_t &size1, int64_t &size2, int64_t &size3, int64_t N, double *kx, double *ky, double *kz, int Msp)
// Modified version of routine from the FNUFFT project
// https://github.com/flatironinstitute/finufft
{
    double min_kx, max_kx;
    arrayrange(N, kx, &min_kx, &max_kx);
    offset1 = (int64_t)ceil(min_kx - Msp); // -Msp+1 if floor 
    size1 = (int64_t)ceil(max_kx - Msp) - offset1 + 2*Msp;

    double min_ky, max_ky;
    arrayrange(N, ky, &min_ky, &max_ky);
    offset2 = (int64_t)ceil(min_ky - Msp);
    size2 = (int64_t)ceil(max_ky - Msp) - offset2 + 2*Msp;

    double min_kz, max_kz;
    arrayrange(N, kz, &min_kz, &max_kz);
    offset3 = (int64_t)ceil(min_kz - Msp);
    size3 = (int64_t)ceil(max_kz - Msp) - offset3 + 2*Msp;

}

void add_subproblem_result(int64_t offset1, int64_t offset2, int64_t offset3, int64_t M1, int64_t M2, int64_t M3, int64_t Mr, complex<double> *ftau, complex<double> *ftau0)
// Modified version of routine from the FNUFFT project
// https://github.com/flatironinstitute/finufft
{
    // Array of offsets for subregion in y
    vector<int64_t> o2(M2);
    int64_t y = offset2;
    for (int i = 0; i < M2; ++i)
    {
        if (y < 0)
            y += Mr;
        if (y >= Mr)
            y -= Mr;
        o2[i] = y++; // Indices of y start-points
    }

    // Array of offsets for subregion in z
    vector<int64_t> o3(M3);
    int64_t z = offset3;
    for (int i = 0; i < M3; ++i)
    {
        if (z < 0)
            z += Mr;
        if (z >= Mr)
            z -= Mr;
        o3[i] = z++; // Indices of z start-points
    }

    // Offsets in x
    int64_t nlo = (offset1 < 0) ? -offset1 : 0;
    int64_t nhi = (offset1 + M1 > Mr) ? offset1 + M1 - Mr : 0; // " above in x


    for (int dz = 0; dz < M3; dz++)
    {
        int64_t oz = Mr * Mr * o3[dz];
        for (int dy = 0; dy < M2; dy++)
        {
            int64_t oy = oz + Mr * o2[dy]; // offsets in y & z

            complex<double> *out = ftau + oy; // 2 * oy;
            complex<double> *in = ftau0 + M1 * (dy + M2 * dz); 
            //2 * M1 * (dy + M2 * dz); // ptr to subgrid array

            // x-indices of ghost points low
            int64_t o = offset1+Mr;
            for (int j = 0; j < nlo; j++) 
                out[j + o] += in[j];

            // x-indices of in-domain points
            o = offset1;
            for (int j = nlo; j < (M1 - nhi); j++)
                out[j + o] += in[j];

            // x-indices of ghost points high
            o = (offset1 - Mr);
            for (int j = (M1 - nhi); j < M1; j++)
                out[j + o] += in[j];
        }
    }
}

// static inline void set_kernel_args(double *args, double x, int ns)
// // Fills vector args[] with kernel arguments x, x+1, ..., x+ns-1.
// // needed for the vectorized kernel eval of Ludvig af K.
// {
//     for (int i = 0; i < ns; i++)
//         args[i] = x + (double)i;
// }