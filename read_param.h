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
int Np;
double Rp, Wp; 
int rand_seed;
int vseed;
double iRx, iRy, iRz; 
double tmin, tmax, dt; 
int increment; 
double accuracy; 
double init_dir[3], dir1[3], dir2[3]; 
bool surfdegen, uppersurf, lowersurf, DoubleU;
int save_every, check_every, stopat;

bool read_param();
//external functions
extern bool norm_v(double* vec);
//external variables
extern int myid;