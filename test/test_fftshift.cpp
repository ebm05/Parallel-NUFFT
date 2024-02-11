#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex>
#include "deconv3d.h"

int main(int argc, char *argv[])
{
    int M = 4;
    complex<double> *Ik1 = (complex<double>*)malloc(sizeof(complex<double>)*M*M*M);
    complex<double> *Ik2 = (complex<double>*)malloc(sizeof(complex<double>)*M*M*M);

    for (int m=0; m<M*M*M; m++){
        Ik1[m] = m+1;
    }

    fftshift(Ik1, Ik2, M, M, M);

    for (int m=0; m<M*M*M; m++){
        printf("IK1[%d]=%f --> %f\n",m,real(Ik1[m]),real(Ik2[m]));
    }

    // Printing indices of multidimensional array. 
    for (int i1=0; i1<M; i1++){
            for (int i2=0; i2<M; i2++){
                    for (int i3=0; i3<M; i3++){
                        int m = i1 + M*i2 + M*M*i3;
                        printf("IK1[%d,%d,%d]=%f\n",i1,i2,i3,real(Ik1[m]));
                    }
            }
    }
    
    return 0;
}