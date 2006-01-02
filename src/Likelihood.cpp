#include "global.h"
#include "List.h"
#include "Node.h"
#include "Rlob.h"
#include "Funs.h"
#include "Likelihood.h"
#include "Prior.h"
#include "BirthDeath.h"
#include "ChangeRule.h"
#include "Swap.h"

double LogLNode(Node *n)
//returns of integrated likelihood of data at node n
{	
	int nob;
	int *divec;
	double lp=0.0;
	
	MakeIntVec(&(n->DataList),&divec,&nob);

	if(nob==0) lp = -10000000.0;
	else
	{

		endNodeModel->setData(nob,XDatR,YDat1,divec,weights);
		lp = endNodeModel->getLogILik();
	}

	delete [] divec;

	return lp;
}


double LogLT(Node *branch,Node *top)
{
	int nbot;
	NodeP *botvec;
	MakeBotVec(branch,&botvec,&nbot);

	int i;
	double LP=0.0;

	for(i=1;i<=nbot;i++) {

		LP += LogLNode(botvec[i]);	
	}

	delete [] botvec;
	return LP;
}
