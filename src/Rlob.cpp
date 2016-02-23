//#include <stdio.h>
#include <cstdio>
#include <iostream>
//#include <math.h>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib>

#include "Rlob.h"

//extern "C" {
#include <R.h>
#include <Rmath.h>
//};

typedef double *dp;
typedef int *ip;

double **almat(long n, long m)
{
	double *temp,**p;
	temp = new double [n*(m+1)];
	if(!temp) Rprintf("allocation error for temp in almat\n");
	p = new dp [n+1];
	if(!p) Rprintf("allocation error for p in almat\n");
	for(int i=1; i<=n; i++) *(p+i) = temp + (m+1)*(i-1);
	return p;
}

void dealmat(double **m)
{
	delete [] m[1];
	delete [] m;
}

int **ialmat(long n, long m)
{
	int *temp,**p;
	temp = new int [n*(m+1)];
	if(!temp) Rprintf("allocation error for temp in ialmat\n");
	p = new ip [n+1];
	if(!p) Rprintf("allocation error for p in almat\n");
	for(int i=1; i<=n; i++) *(p+i) = temp + (m+1)*(i-1);
	return p;
}

void idealmat(int **m)
{
	delete [] m[1];
	delete [] m;
}


void print_mat(long n,double **m)
{
	int i,j;
	Rprintf("\n");
	for(i=1;i<=n;i++) {
		for(j=1;j<=n;j++) {
			Rprintf("%f ",m[i][j]);
		}
		Rprintf("\n");
	}
	Rprintf("\n");
}

void print_mat(FILE *fp,long n,double **m)
{
	int i,j;
	//fprintf(fp,"\n");
	for(i=1;i<=n;i++) {
		for(j=1;j<=n;j++) {
			fprintf(fp,"%f ",m[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

void print_mat(long n,long k,double **m)
{
	int i,j;
	Rprintf("\n");
	for(i=1;i<=n;i++) {
		for(j=1;j<=k;j++) {
			Rprintf("%f ",m[i][j]);
		}
		Rprintf("\n");
	}
	Rprintf("\n");
}

void print_mat(FILE *fp,long n,long k,double **m)
{
	int i,j;
	//Rprintf("\n");
	for(i=1;i<=n;i++) {
		for(j=1;j<=k;j++) {
			fprintf(fp,"%f ",m[i][j]);
		}
		fprintf(fp,"\n");
	}
	//Rprintf("\n");
}

void print_vec(long n,double *vec)
{
	int i;
	for(i=1;i<=n;i++) Rprintf("%f\n",vec[i]);
}
void print_vec(FILE *fp,long n,double *vec)
{
	int i;
	for(i=1;i<=n;i++) fprintf(fp,"%f ",vec[i]);
	fprintf(fp,"\n");
}

int Bern(double p)
//Bernouill random number generator
{
        if(unif_rand() < p) {
                return 1;
        } else {
                return 0;
        }
}

int Disc(double *p)
// draw from discrete distributin given by p, return index
{
   double sum;
   double u = unif_rand();

   int i=1;
   sum=p[1];
   while(sum<u) {
      i += 1;
      sum += p[i];
   }
   return i;
}

double min(double a,double b)
{
	if(a<b) {
		return a;
	} else {
		return b;
	}
}

double  max(double a, double b)
{
	if(a>b) {
		return a;
	} else {
		return b;
	}
}

double myDoubleAbs(double a)
{
	if(a<0) return -a;
	else return a;
}

void indtd(int k,int ind,int *d)
{
int j,i,nind;

nind=ind;

for(i=0;i<k;i++) {
	j=k-i-1;
	d[j+1] = (int)(((double) nind)/(pow(2.0,(double)j)));
	nind=nind-d[j+1]*(int)pow(2.0,(double)j);
}
}

int dtind(int k,int *d)
{
	int ind=0;
	int i;
	for(i=1;i<=k;i++) ind += d[i] * (int)pow(2.0,(double)(i-1));
	return ind;
}

int ISum(int n, int *Iv)
{
	int sum=0;
	int i;
	for(i=1;i<=n;i++) sum += Iv[i];
	return sum;
}

// added for hierarchical -----------------------------------------------------------------------
void mul_ltl(int n,double **l,double **a)
{
	int i,j,k;
	
	for(i=1;i<=n;i++) {
		a[i][i] = 0.0;
		for(k=i;k<=n;k++) a[i][i] += l[k][i]*l[k][i];
		for(j=(i+1);j<=n;j++) {
			a[i][j] = 0.0;
			for(k=i;k<=n;k++) a[i][j] += l[k][i]*l[k][j];
			a[j][i] = a[i][j];
		}
	}
}

void choldc(double **a, int n, double p[])
{
        
        int i,j,k;
        double sum;
 
        for (i=1;i<=n;i++) {
                for (j=i;j<=n;j++) {
                        for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
                        if (i == j) {
                                if (sum <= 0.0)
                                        Rprintf("choldc failed\n");
                                p[i]=sqrt(sum);
                        } else a[j][i]=sum/p[i];
                }
        }
}

void sym_chol_inv(int n,double **a,double **li)
//
// put inverse of the choleski root of a into li
//
{
	
	int i,j,k;
	
	// copy a into li
	for(i=1;i<=n;i++) {
		li[i][i] = a[i][i];
		for(j=(i+1);j<=n;j++) li[i][j] = a[i][j];
	}
	
	double *p = new double [n+1];
	
	// get lower triangular chol root of ai=a
	choldc(li,n,p);
	
	
	double sum;
	
	//get inverse of li = chol(a)
	for(i=1;i<=n;i++) {
		li[i][i] = 1.0/p[i];
		for(j=(i+1);j<=n;j++) {
			sum = 0.0;
			for(k=i;k<j;k++) sum -= li[j][k]*li[k][i];
			li[j][i] = sum/p[j];
			li[i][j] = 0.0;
		}
	}
	
	delete [] p;

	
}

double sym_inv_det(int n,double **a, double **ai)
//
// invert a symmetric matrix
// and return determinant of the inverse
//
{
	double **li = almat(n,n);
	sym_chol_inv(n,a,li);
	mul_ltl(n,li,ai);

	double det=1.0;
	int i;
	for(i=1;i<=n;i++) det *= li[i][i];
	
	dealmat(li);

	return det*det;
	
}

void solve_rtxb(int n,double **r,double *x,double *b)
/* solves r'x=b for x where r is n by n upper triangular */
{
        int k,i,j;
        double d;
        
        k=n;
        
        for(i=1;i<=k;i++)
        {
                d=b[i];
                for(j=1;j<=(i-1);j++)
                {
                        d=d-r[j][i]*x[j];
                }
                x[i]=d/r[i][i];
        }
        return;
}

void solve_rxb(int n, double **r,double *x,double *b)
/* solves rx=b for x where r is n by n and upper triangular */
{
        
        int i,j,ii,k;
        double d;
        
        k=n;
        
        for(i=1;i<=k;i++)
        {
                ii=k-i+1;
                d=b[ii];
                for(j=(ii+1);j<=k;j++)
                {
                        d=d-r[ii][j]*x[j];
                }
                x[ii]=d/r[ii][ii];
        }
        return;
}

double gammln(double xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}

int compare( const void *arg1, const void *arg2 )
{
	double d1,d2;
	d1 = *(double *)arg1;
	d2 = *(double *)arg2;

	if(d1<d2) return -1;
	if(d1>d2)
	{
		return 1;
	} 
	else {
		return 0;
	}
}


void stanAndSortForCart(int n, int k, double **raw, double **stan, int *numUnique, double **uniqueVals, double* meanV, double* scaleV)
// standardize each column and get the sorted values:
{
	

	int i,j;
	double **sortStan=almat(n,k);

	// first you have to sort the orginal values (also get the mean values)

	double *temp = new double[n]; // need to copy into this vec which indexes from 0
	
	for(i=1;i<=k;i++) meanV[i]=0.0;
	for(i=1;i<=k;i++)
	{
		for(j=1;j<=n;j++) temp[j-1] = raw[j][i];
		qsort((void *)temp, n, sizeof(double), compare);
		for(j=1;j<=n;j++) sortStan[j][i] = temp[j-1];

		for(j=1;j<=n;j++) meanV[i] += raw[j][i];
		meanV[i] /= n;
	}
	delete [] temp;

	// now standardize both the original values and the sorted values
	
	for(i=1;i<=k;i++)
	{
		scaleV[i] = (sortStan[n][i]-sortStan[1][i]);
		for(j=1;j<=n;j++)
		{
			
			sortStan[j][i] = (sortStan[j][i]-meanV[i])/scaleV[i];
			stan[j][i] = (raw[j][i]-meanV[i])/scaleV[i];
		}
	}

	

	// now get the unique values-------------------------------------------
	int countU;
	for(i=1;i<=k;i++)
	{
		countU=1;
		uniqueVals[1][i] =  sortStan[1][i];
		for(j=2;j<=n;j++)
		{
			if(sortStan[j][i]!=sortStan[j-1][i])
			{
				countU++;
				uniqueVals[countU][i]= sortStan[j][i];
			}
		}
		numUnique[i]=countU;
		
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	dealmat(sortStan);

}




//------------------------------------------------------------------------------------------------
