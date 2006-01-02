#include <fstream>
#include <cstring>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <search.h>
#include <vector>

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
	   int *iNTree, int *indPost, 
	   int *iprintevery, int *ikeepevery, int *ikeeptrainfits,
	   int *inumcut,
	   double *sdraw, double *trdraw, double *tedraw, int *vcdraw)
{
   GetRNGstate();
   Rprintf("\n\nRunning BART\n\n");

   NumObs = *iNumObs;
   NumX = *iNumX;
   int nrowTest = *inrowTest;

   //XDat, YDat, and XTest are copied into other storage below

   double sigma = *isigma;
   int sigdf = *isigdf;
   double sigquant = *isigquant;

   double kfac = *ikfac;	

   int NTree = *iNTree;
   int ndPost = *indPost;
   
   int printevery= *iprintevery;
   int keepevery = *ikeepevery;
   bool keeptrainfits = true;
   if(!(*ikeeptrainfits)) keeptrainfits=false;

   int numcut= *inumcut;

   double musig = .5/(kfac*sqrt((double)NTree)); 
   if(NTree==1) PriParams.power = 1.0;
   else PriParams.power=2.0;
   PriParams.base=.95;

   // Pr(sigma < k) = sigquant <=>
   //        lambda = qchisq(1-sigquant,df)*k^2/nu, nu=sigdf
   double lambda;
   double dsigdf = (double)sigdf;
   double qchi = qchisq(1.0-sigquant,dsigdf,1,0);
   lambda = (sigma*sigma*qchi)/dsigdf;

   Rprintf("number of trees: %d\n",NTree);
   Rprintf("Prior:\n");
   Rprintf("\tk: %lf\n",kfac);
   Rprintf("\tdegrees of freedom in sigma prior: %d\n",sigdf);
   Rprintf("\tquantile in sigma prior: %lf\n",sigquant);
   Rprintf("data:\n");
   Rprintf("\tnumber of training observations: %d\n",NumObs);
   Rprintf("\tnumber of test observations: %d\n",nrowTest);
   Rprintf("\tnumber of explanatory variables: %d\n\n",NumX);

   Ivec = new int[NumObs+1];
   for(int i=1;i<=NumObs;i++) Ivec[i]=i;
   NumY=1;

   XDat = Lib::almat(NumObs,NumX);
   YDat = Lib::almat(NumObs,NumY); 
   YDat1 = new double[NumObs+1];
   double *Y = new double[NumObs+1];

   int tcnt = 0;
   for(int j=1;j<=NumX;j++) {
      for(int i=1;i<=NumObs;i++) {
         XDat[i][j] = *(iXDat+tcnt); 
              tcnt++;
      }
   }
   for(int i=0;i<NumObs;i++) {
      YDat[i+1][1] = iYDat[i];
      YDat1[i+1] = YDat[i+1][1];
      Y[i+1] = YDat1[i+1];
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
   for(int i=1;i<=NumX;i++) RuleNum[i] = numcut;
   RuleMat = new dp [NumX+1];
   for(int i=1;i<=NumX;i++)
   {
      int nr = RuleNum[i];
      RuleMat[i] = new double[nr+1];
      for(int j=1;j<=nr;j++) RuleMat[i][j] = -.5 + (j/(1.0+nr));
   }
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
   double *onev = new double[nrowTest+1];
   for(int k=1;k<=nrowTest;k++) onev[k]=1.0;
   
   time_t tp;
   int time1 = time(&tp);

   for (int k=1;k<=ndPost;k++) {
      //if(k%printevery== 0) std::cout << "iteration: " << k << " (of " << ndPost << ")" << std::endl;
      if(k%printevery== 0) Rprintf("iteration: %d (of %d)\n",k,ndPost);
      for(nvs i=1;i<theTrees.size();i++) {
         //for(int j=1;j<=NumObs;j++) {
            //YDat1[j] = Y[j]-mtotalfit[j]+mtrainFits[i][j];
         //}
         F77_CALL(dcopy)(&NumObs,Y+1,&inc,YDat1+1,&inc); //copy Y into YDat1
         F77_CALL(daxpy)(&NumObs,&mone,mtotalfit+1,&inc,YDat1+1,&inc); //subtract mtotalfit from YDat1
         F77_CALL(daxpy)(&NumObs,&pone,mtrainFits[i]+1,&inc,YDat1+1,&inc);//add mtrainFits[i]
	 alpha =  Metrop(&theTrees[i],&Done,&step);
         theTrees[i]->currentFits(&mu,NumObs,XDat,YDat1,nrowTest,XTest,weights,mfits);
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
      F77_CALL(dcopy)(&NumObs,Y+1,&inc,eps+1,&inc);
      F77_CALL(daxpy)(&NumObs,&mone,mtotalfit+1,&inc,eps+1,&inc);
      sd.setData(NumObs,eps);
      sd.drawPost();
      mu.setSigma(sd.getS());
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
         sdraw[scnt] = sd.getS(); scnt++;
         countVarUsage(theTrees,varcnt);
         for(int i=1;i<=NumX;i++) { vcdraw[vcnt] = varcnt[i]; vcnt++;} 
      }
   }
   int time2 = time(&tp);
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


