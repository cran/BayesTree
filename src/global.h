#define CAT 1  //  1 is flag for categorical type of variable
#define ORD 2 //   2 is flag for ordered variable

#define BIRTH 1
#define DEATH 2
#define SWAP 3
#define CHANGE 4

#define NORMAL_REG		  1
#define LOGIT_REG		  2
#define POISSON_REG		  3

extern double pBD;
extern double pSwap;
extern double pChange;

extern int NumX;
extern int *VarType;
extern int NumY;

extern int NumObs;
extern double **XDat;
extern double **YDat;
extern double* YDat1;
extern double **XDatR;
extern int NumXR;
extern double* weights;

extern int *RuleNum; 
extern double **RuleMat; 

extern int *Ivec; //global vector defined to be 1,2,3,...NumObs in main

#include "CPriParams.h"
#include "EndNodeModel.h"

extern CPriParams PriParams;
extern EndNodeModel* endNodeModel;


