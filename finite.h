#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#define root 0

extern int Nx, Ny, Nz; 
extern double Lx, Ly, Lz; 
extern double W, U, U2; 
extern double L1, L2, L3, L4; 
extern int chiral; 
//int chiral;
extern double qch;
extern int geo, degenerate, infinite; 
extern int Np;

extern double Rp, Wp; 

// ************************ // 
extern double iRx, iRy, iRz; //aquí está agregado el segundo radio.
// ************************ // 

//int pdegenerate, pinfinite;
extern int rand_seed;
extern int seed;

//New variables for backup and saving
int save_every, check_every;
int stopat;

//Redshift
extern double redshift;

//Volúmenes de las regiones
double dVi, dVo;

extern double tmin, tmax, dt; 
extern int increment; 
extern double accuracy; 
extern double init_dir[3], dir1[3], dir2[3]; 
extern bool uppersurf, lowersurf, DoubleU; //Variables introducidas para nuevo código.

/*SurfDegen por defecto debe ser 1. Es el flag para que comience a trabajar de manera regular. 
Cuando sea 0, significa que las superficies degenerarán de diferentes formas.
Por ejemplo, la superficie superior no degenerará y la inferior sí.
Por ahora se mantiene como estándar la superficie de arriba no degenerada.*/
extern bool surfdegen;

int bulk, surf, tot, droplet;
double S, S2;
double dV, dAdrop, dApart;
double idx, idy, idz, iddx, iddy, iddz;
double *nu;
double  *Qold, *Qo;
int *neighbor;
bool *drop, *boundary;
bool *nboundary;
int *share;
int *bulktype_MPI, *bulktype, *init_bulktype;
bool AnchNInf;

double en_ldg[3];
double en_tot, dE, el_old;
int cycle;
double en_el[5], en_surf[2];
double en_el_in[5], en_el_out[5];

int myid, numprocs;
double *q, *qn, *qo_p;
double *nu_p;
int *neigb, *sign;
int length;
MPI_Win win, win2;
MPI_Comm shmcomm; 

bool read_param();
bool read_nppos(double** pos);
bool initial();
bool initial_bulk();
bool initial_channel();
bool initial_sandwich();
bool initial_cylinder();
bool initial_coaxialcyl();
bool initial_halfcylinder();
bool initial_quartercylinder();
bool initial_drop();
bool initial_halfdrop();
bool initial_quarterdrop();
bool initial_ellip();
int peri(int node, int dir);
bool conf(double** pos);
void free_q();
bool norm_v(double *vec);
double dir2ten(double *dir, int n, double Sin);
double trqq(double *Q);
double trqqq(double *Q);
double matr_mult(double *vec);
double q_mult(double *q1, double *q2);
bool checktr(double *Q);
void free_energy();
void output();
bool scatter();
//double energy_ldg();
void energy_ldg(double* ans);
//void energy_el(double* ans);
void energy_el(double* ans, double* ans_in, double* ans_out);
void energy_surf(double* ans);
void relax_bulk();
void relax_surf();
void en_degen(double* Qin, double* loc_nu, double* Qdiff);
void relax_degen(double* Qin, double* loc_nu, double* Qdiff);