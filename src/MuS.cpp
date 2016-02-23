#include <iostream>
#include <cmath>
#include "MuS.h"
#include "Lib.h"

//extern "C" {
#include <R.h>
//};

// public methods -------------------------------------------------
void MuS::drawPost()
{
   mu = post_m + post_s*norm_rand();
}
double MuS::getLogILik()
{
   //double pi=3.14159265359;  
   double res = .5*log(a/(a+b));
   //std::cout << "p1: " << res << std::endl;
   //res -= (nob/2.0)*std::log(2*pi); //kill later
   //std::cout << "p2: " << res << std::endl;
   //res -= (nob/2.0)*log(sigma2); //kill later
   //std::cout << "p3: " << res << std::endl;
   res -= s2/(2.0*sigma2);
   //std::cout << "p4: " << res << std::endl;
   res -= .5*(a*b*ybar*ybar)/(a+b);
   //std::cout << "p5: " << res << std::endl;
   return res;  
}
void MuS::updatepost()
{
   int i;
   double d;
   if(nob) {
      ybar=0.0;
      for(i=1;i<=nob;i++) ybar += y[indices[i]];
      ybar /= nob;
      s2=0.0;
      for(i=1;i<=nob;i++) {d=y[indices[i]]-ybar; s2 += d*d;}
      b = nob/sigma2;
      post_m = (b*ybar)/(a+b);
      post_s = 1.0/std::sqrt(a+b);
   }
   else {
      post_m=0.0;
      post_s= 1.0/std::sqrt(a);
      b=0.0;
   }
}
void MuS::setData(int nob, double **x, double *y,
                  int *indices, double *w)
{
   this->nob=nob;
   this->y=y;
   this->indices=indices;
   updatepost();
}
double* MuS::getParameterEstimate()
{
   double *rv = new double[2];
   rv[1] = post_m;
   return rv;
}
double* MuS::getFits(int np, double **xpred, int* indpred)
{
   double *rv = new double[np+1];
   for(int i=1;i<=np;i++) rv[i]=post_m;
   return rv;
}
void MuS::toScreen() const
{
   Rprintf("defunct MuS::toScreen called");
   /*
   std::cout << "MuS EndNodeModel:\n";
   std::cout << "  Normal mean given the standard deviation\n";
   std::cout << "mu: " << mu << " ,sigma: " << std::sqrt(sigma2) << "\n";
   std::cout << "prior sd: " << 1.0/std::sqrt(a) << "\n";
   std::cout << "nob: " << nob << std::endl;
   if(nob) { 
      std::cout << "posterior mean and sd: " <<
          post_m << " , " << post_s << "\n";
      std::cout << "nob, ybar,s2: " << nob << ", " <<  ybar << ", " << s2 << std::endl;
   }
   */
}
