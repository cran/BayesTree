#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <search.h>
#include <vector>
#include <valarray>
#include <algorithm>

extern "C" {
#include <R.h>
#include <Rmath.h>
};

#include "global.h"
#include "Node.h"
#include "Funs.h"
#include "Prior.h"
#include "MuS.h"
#include "Sdev.h"

extern "C" {
void F77_NAME(dcopy)(const int *n, const double *dx, const int *incx, double *dy, const int *incy);
void F77_NAME(daxpy)(const int *n, const double *alpha, const double *dx, const int *incx,
                double *dy, const int *incy);
double F77_NAME(ddot)(const int *n, const double *dx,
                  const int *incx, const double *dy, const int *incy);
};

using namespace std;

typedef double *dp;
typedef int *ip;
typedef Node *NodeP;
typedef void *voidP;
typedef EndNodeModel* EndNodeModelP;

//define global variables -----------------------------------------------------------
int NumX; // number of x variables
int NumY; //number of y variables
int *VarType; // for each variable tell what kind, CAT or ORD
int NumObs; // number of observations
double **XDat; // x data, note: cats are double
double **YDat; // y data, note: cats are double
double* YDat1=0;	// y data with just one y
double **XDatR;	// x data, used in regression model
int NumXR;		// number of columns in XDatR
double* weights;

int *RuleNum; // integer vec of length NumX, ith is number of split
				// points for ORD var and number of CATs for CAT var
double **RuleMat; // ragged array, ith row has RuleNum values, split values
					// for ORD, cat values for CAT

int *Ivec;

double pBD;
double pSwap;
double pChange;


CPriParams PriParams;
EndNodeModel* endNodeModel=0;

extern "C" {
void mbart(int *iNumObs, int *iNumX, int *inrowTest,
           double *iXDat, double *iYDat,
	   double *iXTest,
	   double *isigma, int *isigdf, double *isigquant,
	   double *ikfac,
	   double *ipower, double *ibase,
	   double *ibinary_offset,
	   int *iNTree, int *indPost, 
	   int *iprintevery, int *ikeepevery, int *ikeeptrainfits,
	   int *inumcut, int *iusequants, int *iprintcutoffs,
	   int *verbose,
	   double *sdraw, double *trdraw, double *tedraw, int *vcdraw)
{
   GetRNGstate();

   bool binary = (*ibinary_offset > -1000.0);
   double binary_offset = *ibinary_offset;

   if(*verbose) {
      if(binary)
         Rprintf("\n\nRunning BART with binary y\n\n");
      else
         Rprintf("\n\nRunning BART with numeric y\n\n");
   }

   NumObs = *iNumObs;
   NumX = *iNumX;
   int nrowTest = *inrowTest;

   //XDat, YDat, and XTest are copied into other storage below

   double sigma;
   if(binary)
      sigma = 1.0;
   else
      sigma = *isigma;

   int sigdf = *isigdf;
   double sigquant = *isigquant;

   double kfac = *ikfac;	

   int NTree = *iNTree;
   int ndPost = *indPost;
   
   int printevery= *iprintevery;
   int keepevery = *ikeepevery;
   bool keeptrainfits = true;
   if(!(*ikeeptrainfits)) keeptrainfits=false;
   bool usequants = true;
   if(!(*iusequants)) usequants=false;
   int printcutoffs = *iprintcutoffs;

   //note: the meaning of kfac is different for binary y 
   //  for numeric y, it is standardized so E(Y) is probably in (-.5,.5)
   //  for binary y, it is in (-3,3)
   double musig;
   if(binary)
      musig = 3.0/(kfac*sqrt((double)NTree)); 
   else
      musig = .5/(kfac*sqrt((double)NTree)); 

   PriParams.base = *ibase;
   PriParams.power = *ipower;

   // Pr(sigma < k) = sigquant <=>
   //        lambda = qchisq(1-sigquant,df)*k^2/nu, nu=sigdf
   double lambda;
   double dsigdf = (double)sigdf;
   double qchi = qchisq(1.0-sigquant,dsigdf,1,0);
   lambda = (sigma*sigma*qchi)/dsigdf;

   if(*verbose) {
   Rprintf("number of trees: %d\n",NTree);
   Rprintf("Prior:\n");
   Rprintf("\tk: %lf\n",kfac);
   if(binary) {
      Rprintf("\tbinary offset is: %lf\n",binary_offset);
   } else {
      Rprintf("\tdegrees of freedom in sigma prior: %d\n",sigdf);
      Rprintf("\tquantile in sigma prior: %lf\n",sigquant);
   }
   Rprintf("\tpower and base for tree prior: %lf %lf\n",PriParams.power,PriParams.base);
   Rprintf("\tuse quantiles for rule cut points: %d\n",usequants);
   Rprintf("data:\n");
   Rprintf("\tnumber of training observations: %d\n",NumObs);
   Rprintf("\tnumber of test observations: %d\n",nrowTest);
   Rprintf("\tnumber of explanatory variables: %d\n\n",NumX);
   }

   Ivec = new int[NumObs+1];
   for(int i=1;i<=NumObs;i++) Ivec[i]=i;
   NumY=1;

   XDat = Lib::almat(NumObs,NumX);
   YDat = Lib::almat(NumObs,NumY); //note: for binary y this the 0-1 y (never changes), for continuous never used
   YDat1 = new double[NumObs+1];//used for resids in backfitting
   double *Y = new double[NumObs+1];//used for y in numeric, latent - binary_offset for binary y

   int tcnt = 0;
   for(int j=1;j<=NumX;j++) {
      for(int i=1;i<=NumObs;i++) {
         XDat[i][j] = *(iXDat+tcnt); 
              tcnt++;
      }
   }
   if(binary) {
      for(int i=0;i<NumObs;i++) {
         YDat[i+1][1] = iYDat[i];
	 Y[i+1] = 2*(iYDat[i]) - 1.0 - binary_offset; //guess at good starting for latents
      }
   } else {
      for(int i=0;i<NumObs;i++) Y[i+1] = iYDat[i];
   }


   double** XTest = Lib::almat(nrowTest,NumX);
   tcnt = 0;
   for(int j=1;j<=NumX;j++) {
      for(int i=1;i<=nrowTest;i++) {
         XTest[i][j] = *(iXTest+tcnt);
         tcnt++;
      }
   } 

   VarType = new int [NumX+1];
   for(int i=1;i<=NumX;i++) VarType[i]=ORD;
   weights = new double[NumObs+1];
   for(int i=1;i<=NumObs;i++) weights[i]=1.0;
   RuleNum = new int[NumX+1];
   RuleMat = new dp [NumX+1];

   int cnq,cfac,cnc,cind,coffset;
   double maxx,minx,xinc;
   double* xcol;
   if(usequants) {
      Vec xv;
      for(int i=1;i<=NumX;i++) {
         Lib::sortedUnique(NumObs,iXDat+((i-1)*NumObs),xv);
         cnq = (int)xv.size();
         if(cnq<=(inumcut[i-1]+1)) {
            cfac = 1;
            cnc = cnq-1;
            coffset=0;
         } else {
            cnc = inumcut[i-1];
            cfac = cnq/cnc;
            coffset = (cfac/2);
         }
         RuleNum[i] = cnc; 
         RuleMat[i] = new double[cnc+1];
         for(int j=0;j<cnc;j++) {
            cind = std::min(j*cfac+coffset,cnq-2);
            RuleMat[i][j+1] = (xv[cind]+xv[cind+1])/2.0;
         }
      }
   } else {
      for(int i=1;i<=NumX;i++) {
         xcol = iXDat + ((i-1)*NumObs);
         maxx = *(std::max_element(xcol,xcol+NumObs));
         minx = *(std::min_element(xcol,xcol+NumObs));
         RuleNum[i] = inumcut[i-1];
         xinc = (maxx-minx)/(RuleNum[i]+1);
         RuleMat[i] = new double[RuleNum[i]+1];
         for(int j=1;j<=RuleNum[i];j++) RuleMat[i][j] = minx + j*xinc;
      }
   }
   if(*verbose) {
   Rprintf("\nCutoff rules c in x<=c vs x>c\n");
   Rprintf("Number of cutoffs: (var: number of possible c):\n");
   for(int i=1;i<=NumX;i++) {
      Rprintf("(%d: %d) ",i,RuleNum[i]);
      if(i%5 ==0) Rprintf("\n");
   }
   Rprintf("\n");
   if(printcutoffs>0) {
      Rprintf("cutoffs:\n");
      for(int i=1;i<=NumX;i++) {
         Rprintf("x(%d) cutoffs: ",i);
         int j=1;
         while((j<RuleNum[i]) && (j<printcutoffs) ) {
            Rprintf("%lf ",RuleMat[i][j]);
            if(j%5 ==0) Rprintf("\n\t");
            j+=1;
         }
         if((j>2) && (j!=RuleNum[i])) Rprintf("...");
         Rprintf("%lf",RuleMat[i][RuleNum[i]]);
         Rprintf("\n");
      }
   }
   Rprintf("\n\n");
   }

   //temporary: print out rules:
   /*
   for(int i=1;i<=NumX;i++) {
      std::cout << "\n\n";
      std::cout << "for var: " << i << ", numrules: " << RuleNum[i] << std::endl;
      std::cout << "cuts:\n";
      for(int j=1;j<=RuleNum[i];j++) {
         std::cout << RuleMat[i][j] << "  ";
         if(j%10 ==0) std::cout << std::endl;
      }
   }
   */

   MuS mu;
   mu.setSigma(sigma);
   mu.setPriorS(musig);
   endNodeModel = &mu;

   Sdev sd;
   sd.setPrior(sigdf,lambda);

   pBD = .5;
   pChange = .4;
   pSwap = .1;

   std::vector<Node*> theTrees(NTree+1);
   typedef std::vector<Node*>::size_type nvs;
   for(nvs i=1;i<theTrees.size();i++)
      {theTrees[i] = new Node; theTrees[i]->SetData();}

   double** mtrainFits = Lib::almat(NTree,NumObs);
   for(int i=1;i<=NumObs;i++)
   {
      for(int j=1;j<=NTree;j++) mtrainFits[j][i] =0.0;
   }
   double** mtestFits = Lib::almat(nrowTest,NTree);
   //double** mtestFits = Lib::almat(NTree,nrowTest);

   int Done=0;
   double alpha=0.0;
   int step=0;

   double* eps = new double[NumObs+1];
   double* mtotalfit = new double[NumObs+1];
   // invariant: the rows sums of mtrainFits must be this mtotalfit
   for(int i=1;i<=NumObs;i++) mtotalfit[i]=0.0;

   double **mfits;
   mfits = new dp[3];
   mfits[1] = new double[NumObs+1];
   mfits[2] = new double[nrowTest+1];

   int scnt=0;    // count draws of sigma
   int trcnt=0;   //count draws of train fits
   int tecnt=0;   //count draws of test  fits
   std::vector<int> varcnt; //store var counts, each draw
   typedef std::vector<int>::size_type ivs;
   int vcnt=0;    // count draws of var counts
   int inc=1;
   double mone=-1.0;
   double pone=1.0;
   double *onev = new double[NTree+1];
   for(int k=1;k<=NTree;k++) onev[k]=1.0;
   
   time_t tp;
   int time1 = time(&tp);

   if(*verbose) Rprintf("Running mcmc loop:\n");
   for (int k=1;k<=ndPost;k++) {
      //if(k%printevery== 0) std::cout << "iteration: " << k << " (of " << ndPost << ")" << std::endl;
      if(*verbose && (k%printevery== 0)) Rprintf("iteration: %d (of %d)\n",k,ndPost);
      for(nvs i=1;i<theTrees.size();i++) {
         //for(int j=1;j<=NumObs;j++) {
            //YDat1[j] = Y[j]-mtotalfit[j]+mtrainFits[i][j];
         //}
         F77_CALL(dcopy)(&NumObs,Y+1,&inc,YDat1+1,&inc); //copy Y into YDat1
         F77_CALL(daxpy)(&NumObs,&mone,mtotalfit+1,&inc,YDat1+1,&inc); //subtract mtotalfit from YDat1
         F77_CALL(daxpy)(&NumObs,&pone,mtrainFits[i]+1,&inc,YDat1+1,&inc);//add mtrainFits[i]
	 alpha =  Metrop(&theTrees[i],&Done,&step);

	 if(k%keepevery==0)
            theTrees[i]->currentFits(&mu,NumObs,XDat,YDat1,nrowTest,XTest,weights,mfits);
	 else
            theTrees[i]->currentFits(&mu,NumObs,XDat,YDat1,0,XTest,weights,mfits);

	 //for(int j=1;j<=NumObs;j++) mtotalfit[j] += (mfits[1][j]-mtrainFits[i][j]);
         F77_CALL(daxpy)(&NumObs,&mone,mtrainFits[i]+1,&inc,mtotalfit+1,&inc);//sub old fits
         F77_CALL(daxpy)(&NumObs,&pone,mfits[1]+1,&inc,mtotalfit+1,&inc); //add new fits
         //for(int j=1;j<=NumObs;j++) mtrainFits[i][j] = mfits[1][j];
         F77_CALL(dcopy)(&NumObs,mfits[1]+1,&inc,mtrainFits[i]+1,&inc);
         for(int j=1;j<=nrowTest;j++) mtestFits[j][i] = mfits[2][j];
         //F77_CALL(dcopy)(&NumObs,mfits[2]+1,&inc,mtestFits[i]+1,&inc);
      }
      //for(int m=1;m<=NumObs;m++) 
      //   eps[m]=YDat[m][1]-mtotalfit[m];
      if(!binary) {
         F77_CALL(dcopy)(&NumObs,Y+1,&inc,eps+1,&inc);
         F77_CALL(daxpy)(&NumObs,&mone,mtotalfit+1,&inc,eps+1,&inc);
         sd.setData(NumObs,eps);
         sd.drawPost();
         mu.setSigma(sd.getS());
      }
      if(binary) {
         double u,Z;
         for(int i=1;i<=NumObs;i++) {
	    u = unif_rand();
	    if(YDat[i][1] > 0) {
	       Z = qnorm((1.0-u)*pnorm(-mtotalfit[i]-binary_offset,0.0,1.0,1,0) + u,0.0,1.0,1,0);
	    } else {
	       Z = -qnorm((1.0-u)*pnorm(mtotalfit[i]+binary_offset,0.0,1.0,1,0) + u,0.0,1.0,1,0);
	    }
	    Y[i] = mtotalfit[i] + Z;
	 }
      }
      if(k%keepevery==0) {
         //double sum;
         if(keeptrainfits) {
            for(int i=1;i<=NumObs;i++) {trdraw[trcnt] = mtotalfit[i]; trcnt++;}
         }
         for(int i=1;i<=nrowTest;i++) {
            //sum=0.0;
            //for(int j=1;j<=NTree;j++) sum += mtestFits[i][j];
            //tedraw[tecnt] = sum; tecnt++;
            tedraw[tecnt++] = F77_CALL(ddot)(&NTree,onev+1,&inc,mtestFits[i]+1,&inc);
         }
         if(!binary) sdraw[scnt] = sd.getS(); scnt++;
         countVarUsage(theTrees,varcnt);
         for(int i=1;i<=NumX;i++) { vcdraw[vcnt] = varcnt[i]; vcnt++;} 
      }
   }
   int time2 = time(&tp);

   if(*verbose) {
   Rprintf("time for loop: %d\n",time2-time1);

   Rprintf("\nTree sizes, last iteration:\n");
   for(nvs i=1;i<theTrees.size();i++) {
      Rprintf("%d ",theTrees[i]->NumBotNodes());
      if(i%20 == 0) Rprintf("\n");
   }
   std::cout << "\n\n";

   Rprintf("Variable Usage, last iteration (var:count):\n");
   countVarUsage(theTrees,varcnt);
   for(ivs i=1;i<varcnt.size();i++) {  
      Rprintf("(%d: %d) ",i,varcnt[i]);
      if(i%5 == 0) Rprintf("\n");
   }

   Rprintf("\nDONE BART\n\n");
   }

   //delete 
   if(nrowTest) {
      Lib::dealmat(mtestFits);
      Lib::dealmat(XTest);
   }
   Lib::dealmat(mtrainFits);
   Lib::dealmat(XDat);
   Lib::dealmat(YDat);
   delete [] YDat1;
   delete [] Y;
   delete [] Ivec;
   delete [] VarType;
   delete [] weights;
   delete [] RuleNum;
   for(int i=1;i<=NumX;i++) delete [] RuleMat[i];
   delete [] RuleMat;
   for(nvs i=1;i<theTrees.size();i++)
      theTrees[i]->deall();
   delete [] eps;
   delete [] mtotalfit;
   delete [] mfits[1];
   delete [] mfits[2];
   delete [] mfits;

   PutRNGstate();
}
};


