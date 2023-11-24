#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#define root 0

int Nx, Ny, Nz; 
double Lx, Ly, Lz; 
double W, U, U2; 
double L1, L2, L3, L4;  
int chiral;
double qch;
double redshift;
int geo; 
int degenerate, infinite; 
double tiltAngle;
int Np;
double Rp, Wp; 
int rand_seed;
int seed;
int conf_seed;
double iRx, iRy, iRz; 
double tmin, tmax, dt; 
double increment; 
double accuracy; 
double init_dir[3], dir1[3], dir2[3]; 
bool surfdegen, uppersurf, lowersurf, DoubleU;
int save_every, check_every, stopat, trace_checker;

bool read_param();
//external functions
extern bool norm_v(double* vec);
//external variables
extern int myid;