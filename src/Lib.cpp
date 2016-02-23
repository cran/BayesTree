#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <fstream>

#include "Lib.h"

//extern "C" {
#include <R.h>
#include <Rmath.h>
//};


//public methods-----------------------------------------------------
void Lib::printVec(std::vector<int>& iv)
{
      Rprintf("defunct Lib::printVec");
        //for(std::vector<int>::const_iterator it = iv.begin(); it !=iv.end(); it++) std::cout << *it << std::endl;
}
void Lib::printVec(Vec& v)
{
      Rprintf("defunct Lib::printVec");
        //for(Vec::const_iterator it = v.begin(); it !=v.end(); it++) std::cout << *it << std::endl;
}
void Lib::printVec(Vec& v , char* fname)
{
        std::ofstream of(fname);
        for(Vec::const_iterator it = v.begin(); it != v.end(); it++) of << *it << std::endl;
}
double** Lib::almat(long n, long m)
{
        double *temp,**p;
        temp = new double [n*(m+1)];
        if(!temp) Rprintf("allocation error for temp in almat\n");
        p = new dp [n+1];
        if(!p) Rprintf("allocation error for p in almat\n");
        for(int i=1; i<=n; i++) *(p+i) = temp + (m+1)*(i-1);
        return p;
}
void Lib::dealmat(double **m)
{
        delete [] m[1];
        delete [] m;
}
int** Lib::ialmat(long n, long m)
{
        int *temp,**p;
        temp = new int [n*(m+1)];
        if(!temp) Rprintf("allocation error for temp in ialmat\n");
        p = new ip [n+1];
        if(!p) Rprintf("allocation error for p in almat\n");
        for(int i=1; i<=n; i++) *(p+i) = temp + (m+1)*(i-1);
        return p;
}
void Lib::idealmat(int **m)
{
        delete [] m[1];
        delete [] m;
}

double Lib::median(Vec x)
{
        typedef std::vector<double>::size_type vec_sz;

        vec_sz size = x.size();
        if(size==0)
                throw std::domain_error("median of an empty vector");

        std::sort(x.begin(),x.end());
        vec_sz mid = size/2;
        return size%2 == 0 ? (x[mid]+x[mid-1])/2 : x[mid];
}
void Lib::sortedUnique(int n, double *x, Vec& uv)
{
   uv.clear();
   if(n==0) return;
   typedef std::vector<double>::size_type vec_sz;

   //copy x into vector xv
   Vec xv(n);
   vec_sz nv = (vec_sz)n;
   for(vec_sz i=0;i<nv;i++) xv[i] = *(x+i);

   //sort xv
   std::sort(xv.begin(),xv.end());

   //get unique values of xv into uv
   uv.push_back(xv[0]);
   vec_sz nu = 1;
   double ov = uv[0];
   for(vec_sz i=1;i<nv;i++) {
      if(xv[i] != ov) {
         ov = xv[i];
         uv.push_back(ov);
      }
   }

}
double Lib::mean(const Vec& x)
{
        return std::accumulate(x.begin(),x.end(),0.0)/x.size();
}
double Lib::sdev(const Vec& x)
{
        double d;
        double m = mean(x);
        double s = 0.0;
        int n = (int)x.size();
        for(int i=0;i<n;i++) { d = x[i]-m; s += d*d;}
        return sqrt(s/n);
}
void Lib::acov(Vec& x,int nl, Vec& acov, bool cor)
{
        typedef Vec::size_type ST;
        double c;
        int n = (int)x.size();
        double m = mean(x);
        acov.clear();
        for(int i=0;i<=nl;i++)
        {
                c=0.0;
                for(int j=0;j<(n-i);j++) c += (x[j]-m)*(x[j+i]-m);
                acov.push_back(c);
        }
        if(cor)
        {
                double c0 = acov[0];
                for(int i=0;i<=nl;i++)
                {
                        acov[i] = acov[i]/c0;
                }
        }
        else
                for(int i=0;i<=nl;i++) acov[i] /= n;

}
double Lib::tssd(Vec& x, int n, int nl)
{
        Vec gamma;
        acov(x,nl,gamma,false);
        double v=gamma[0];
        for(int i=1;i<=nl;i++) v += 2*(1.0-((double)i/n))*gamma[i];
        return sqrt(v/n);
}
void Lib::batchMeans(Vec& x, int bsize, Vec& means)
{
        means.clear();
        int n = (int)x.size();
        bool someleft = (n>=bsize);
        int nbatch=0;
        while(someleft)
        {
                double sum = 0.0;
                int off = nbatch*bsize;
                for(int i=0;i<bsize;i++) sum += x[off+i];
                means.push_back(sum/bsize);
                nbatch += 1;
                someleft = (n>=(bsize)*(nbatch+1));
        }
}
int Lib::Disc(Vec& p)
// draw from discrete distributin given by p, return index
{
    double sum;
    double u = unif_rand();

    int i=0;
    sum=p[0];
    while(sum<u) {
        i += 1;
        sum += p[i];
    }
    return i;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
