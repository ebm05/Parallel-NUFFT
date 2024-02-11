// Definitions
//
#define PI 3.141592653589793
#define M_1_2PI 0.159154943091895336

// affine transform to x:
// when p=true:   map [-3pi,-pi) and [-pi,pi) and [pi,3pi)    each to [0,N)
// otherwise,     map [-N,0) and [0,N) and [N,2N)             each to [0,N)
#define FOLDRESCALE(x,N,p) (p ? (x + (x>=-PI ? (x<PI ? PI : -PI) : 3*PI)) * ((double)M_1_2PI*N) : (x>=0.0 ? (x<(double)N ? x : x-(double)N) : x+(double)N))

// Pseudo random numbers
#define rand01() ((double)rand()/(double)RAND_MAX)
#define randm11() (2*rand01() - (double)1.0)
#define crandm11() (randm11() + IMA*randm11())

// Spreading parameter
#define TAU(M,R,Msp) (1 / ((double)(M * M))) * (PI * Msp) / ((double)R * ((double)R - 0.5))
