#ifndef GUARD_Lib
#define GUARD_Lib

#include <iostream>
#include <string>
#include <vector>
#include <cstring>

//extern "C" {
#include <R.h>
#include <Rmath.h>
//};

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat>

typedef double *dp;
typedef int *ip;
typedef std::vector<double> Vec;
typedef std::vector<int> IVec;

class Lib
{
public:
        //constructors, destructor-----------------------------------
        Lib() {}
        ~Lib() {}
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // static public methods-----------------------------------------
        static double** almat(long n, long m);
        static void dealmat(double** mat);
        static int** ialmat(long n, long m);
        static void idealmat(int** mat);

        static double median(Vec x);
        static void sortedUnique(int n, double *x, Vec& uv);
        static double mean(const Vec& x);
        static double sdev(const Vec& x);
        static void acov(Vec& x,int nl, Vec& acov, bool cor=true);
        static double tssd(Vec& x,int n, int nl);
        static void batchMeans(Vec& x, int bsize, Vec& means);

        static void printVec(Vec& v);
        static void printVec(Vec& v, char* fname);
        static void printVec(std::vector<int>& iv);

        static double nextDouble() { return unif_rand();}
        static double nextGaussian() { return norm_rand();}
        static double genChi2(int nu) {return rchisq((double)nu);}
        static int Disc(Vec& p);
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

};


#endif // GUARD_Lib
