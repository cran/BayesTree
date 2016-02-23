#ifndef GUARD_Rlob
#define GUARD_Rlob

//#include <stdio.h>
#include <cstdio>
#include <iostream>
//#include <math.h>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib>

//extern "C" {
#include <R.h>
#include <Rmath.h>
//};

typedef double *dp;
typedef int *ip;

double **almat(long n, long m);
void dealmat(double **m);
int **ialmat(long n, long m);
void idealmat(int **m);

void print_mat(long n,double **m);
void print_mat(long n,long k,double **m);
void print_mat(FILE *fp,long n,double **m);
void print_mat(FILE *fp,long n,long k,double **m);
void print_vec(long n,double *vec);
void print_vec(FILE *fp,long n,double *vec);

//double ran1(long *idum);
//double gasdev(long *idum);

int Bern(double p);
int Disc(double *p);
double min(double a,double b);
double myDoubleAbs(double a);

void indtd(int k,int ind,int *d);
int dtind(int k,int *d);

int ISum(int n, int *Iv);
double max(double a, double b);

void choldc(double **a, int n, double p[]);
void sym_chol_inv(int n,double **a,double **li);
double sym_inv_det(int n,double **a, double **ai);
void mul_ltl(int n,double **l,double **a);
void solve_rtxb(int n,double **r,double *x,double *b);
void solve_rxb(int n, double **r,double *x,double *b);
double gammln(double xx);

int compare( const void *arg1, const void *arg2 );
void stanAndSortForCart(int n, int k, double **raw, double **stan, int *numUnique, double **uniqueVals, double* meanV, double* scaleV);

#endif



