//extern "C" {
#include <R.h>
#include <Rmath.h>
//};


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

typedef Node *NodeP;


void MakeBotVec(Node *top,NodeP **botvec,int *NBot)
// allocates and defines array of bot nodes for tree top, sets NBot
{

	int i;

	List *bots;
	top->GetBotList(&bots);
	*NBot = bots->length;
	*botvec = new NodeP [*NBot+1];
	
		
	Cell *cell = bots->first;
	(*botvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*NBot);i++) {
		cell = cell->after;
		(*botvec)[i]=(Node *)cell->contents;	
	}
	
	bots->deall();
	delete bots;
}

void MakeNogVec(Node *top,NodeP **nogvec,int *NNog)
// allocates and defines list of nog nodes for tree top, sets NNog
{

	int i;

	List *bots;
	top->GetNogList(&bots);
	*NNog = bots->length;
	*nogvec = new NodeP [*NNog+1];
	
	if (*NNog){
	Cell *cell = bots->first;
	(*nogvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*NNog);i++) {
		cell = cell->after;
		(*nogvec)[i]=(Node *)cell->contents;	
	}
	}
	bots->deall();
	delete bots;
}

void MakeSwapVec(Node *top,NodeP **swapvec,int *Nswap)
// 
{

	int i;

	List *swaps;
	top->GetSwapsList(&swaps);
	*Nswap = swaps->length;
	*swapvec = new NodeP [*Nswap+1];
	
	if (*Nswap){
	Cell *cell = swaps->first;
	(*swapvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*Nswap);i++) {
		cell = cell->after;
		(*swapvec)[i]=(Node *)cell->contents;	
	}
	}
	swaps->deall();
	delete swaps;
}

void MakeNotBotVec(Node *top,NodeP **notbotvec,int *Nnotbot)
// allocates and defines list of nog nodes for tree top, sets NNog
{

	int i;

	List *bots;
	top->GetNotBotList(&bots);
	*Nnotbot = bots->length;
	*notbotvec = new NodeP [*Nnotbot+1];
	
	if (*Nnotbot){
	Cell *cell = bots->first;
	(*notbotvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*Nnotbot);i++) {
		cell = cell->after;
		(*notbotvec)[i]=(Node *)cell->contents;	
	}
	}
	bots->deall();
	delete bots;
}


void MakeIntVec(List *intlist, int **ivec, int *n)
//allocate and define vec of integers corresponding to list of int pointers
{
	int i;

	*n = intlist->length;
	*ivec = new int [*n + 1];

	Cell *cell = intlist->first;
	if(*n>0) (*ivec)[1]=*((int *)(cell->contents));
	
	for(i=2;i<=(*n);i++) {
		cell = cell->after;
		(*ivec)[i]=*((int *)(cell->contents));
	}
}

void GetDataInd(Node *top,int *ind)
// sets ind[i] to be index of nodes corresponding to obs i
// ind is already allocated of length NumObs
{

	int i,j;

	int nbot;
	NodeP *botvec;
	MakeBotVec(top,&botvec,&nbot);

	Node *nn;

	for(i=1;i<=NumObs;i++) {
	
		top->FindNode(XDat[i],&nn);
		for(j=1;nn!=botvec[j];j++);
		ind[i] = j;
		
	}

	delete [] botvec;
}
void GetDataInd(Node *top,int *ind, int NumObsPred, double** data)
// sets ind[i] to be index of nodes corresponding to obs i
// ind is already allocated of length NumObsPred
{

	int i,j;

	int nbot;
	NodeP *botvec;
	MakeBotVec(top,&botvec,&nbot);

	Node *nn;

	for(i=1;i<=NumObsPred;i++) {
	
		top->FindNode(data[i],&nn);
		for(j=1;nn!=botvec[j];j++);
		ind[i] = j;
		
	}

	delete [] botvec;
}

int getMaxDepth(Node *top)
{
	int nbot;
	NodeP *botvec;
	MakeBotVec(top,&botvec,&nbot);

	int i;
	int max=0;
	int tDepth;
	for(i=1;i<=nbot;i++)
	{
		tDepth = Depth(botvec[i]);
		if(tDepth>max) max = tDepth;
	}

	delete [] botvec;
	return max;
}




Rule *GetRulePointer(int index, int curindex, int depth, int curdepth, Node* n)
// gets the pointer to the rule of the node having index value index
// index: index of node whose rule we want
// curindex: index of node n
// depth: number of levels below the top
// curdepth: depth of node n
//
// when you call at the top you should have
// curindex = 1
// depth = top->DepthBelow()
// curdepth = 0
// node = top

{

	//return &(n->rule);

	int nn;


	// if you are a bottom node then you have no rule so return 0 for the pointer value
	if(n->Bot) return 0;
	
	
	// if you are at the right node return the address of the rule
	if(index==curindex)
	{
		return &(n->rule);
	}
	// let d= depth-curdepth = # of levels below n
	// the number of nodes in each branch below n is 2^d -1
	// the index of the left child is the curindex + 1
	// the index of the left is current + #below to the left + 1 = curindex + 2^d
	// either go left or right
	else
	{

	    nn = (int)pow(2.0,(double)(depth-curdepth));

		if(index >= (curindex + nn))
		{
			return GetRulePointer(index,curindex+nn,depth,curdepth+1,n->RightC);
		}
		else
		{
			return GetRulePointer(index,curindex+1,depth,curdepth+1,n->LeftC);
		}
	}
}

int RulesDifferent(Rule *r1,Rule *r2)
{
	if((r1==0) && (r2==0)) return 0;

	if((r1==0) && (r2!=0)) return 1;

	if((r1!=0) && (r2==0)) return 1;

	int v1 = r1->Var;
	int v2 = r2->Var;

	if(v1!=v2)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

		
int ShannonBanksMetric(Node *top1,Node *top2)
{
	int d1 = top1->DepthBelow();
	int d2 = top2->DepthBelow();

	int depth = (int)max(d1,d2);

	int numnodes = (int)pow(2.0,depth+1.0) - 1;

	int index;
	Rule *r1;
	Rule *r2;

	int numdiff=0;

	for(index=1;index<=numnodes;index++)
	{
		r1 = GetRulePointer(index, 1, depth, 0, top1);
		r2 = GetRulePointer(index, 1, depth, 0, top2);

		numdiff += RulesDifferent(r1,r2);
	}

	return numdiff;
}



/*double FitMetric(Node* top1,Node* top2)
//gets the fits for the two trees and sums the squared differences between the fitted values
{
	double* f1 = Fits(top1);
	double* f2 = Fits(top2);

	int i;
	double met=0;
	for(i=1;i<=NumObs;i++)
	{
		met += (f1[i]-f2[i])*(f1[i]-f2[i]);
	}

	delete [] f1;
	delete [] f2;

	return met;
}

int MisCMetric(Node* top1,Node* top2)
//gets the fits for the two trees and sums the squared differences between the fitted values
{
	double* f1 = Fits(top1);
	double* f2 = Fits(top2);

	int i;
	int met=0;
	int r1,r2;
	for(i=1;i<=NumObs;i++)
	{
		r1=0;
		r2=0;
		if(f1[i]>.5) r1=1;
		if(f2[i]>.5) r2=1;
		if(!(r1==r2)) met += 1;
		
	}

	delete [] f1;
	delete [] f2;

	return met;
}*/

int AndrewsMetric(Node *top1,Node *top2)
// for each pair of observations count the number of times they are not "partioned the same"
// in the two trees.
// they are partioned the same if whether or not the two observations are in the same
// partition is the same for both trees
{
	int *ind1 = new int [NumObs+1];
	int *ind2 = new int [NumObs+1];

	GetDataInd(top1,ind1);
	GetDataInd(top2,ind2);

	int count=0;
	int Isame1,Isame2;

	int i,j;

	for(i=1;i<NumObs;i++)
	{
		for(j=(i+1);j<=NumObs;j++)
		{
			Isame1 = (ind1[i]==ind1[j])?1: 0;
			Isame2 = (ind2[i]==ind2[j])?1: 0;

			count += std::abs(Isame1-Isame2);

		}
	}

	delete [] ind1;
	delete [] ind2;

	return count;
}


void AddDatChildren(Node *n)
{
	if(!(n->rule).Var) 
		Rprintf("error in AddDatChildren: rule not set\n");
	if(((n->LeftC)->DataList.length!=0) || ((n->RightC)->DataList.length!=0))
		Rprintf("error in AddDatChildren: data already set\n");

	int *divec;
	int nob;
	MakeIntVec(&(n->DataList),&divec,&nob);

	int i;
	for(i=1;i<=nob;i++) {
		if((n->rule).Right(XDat[divec[i]])) {
			(n->RightC)->SetData(divec[i]);
			
		}
		else {
			(n->LeftC)->SetData(divec[i]);
		}
	}
		
	
	delete [] divec;
}

void FixDataBelow(Node *cnode)
{
	int *divec;
	int nobs;
	int i;
	
	(cnode->LeftC)->ClearData();
	(cnode->RightC)->ClearData();

		MakeIntVec(&(cnode->DataList), &divec, &nobs);
	for(i=1;i<=nobs;i++) {
		if ((cnode->rule).Right(XDat[divec[i]])) {
			(cnode->RightC)->SetData(divec[i]);
		} else {
			(cnode->LeftC)->SetData(divec[i]);
		}
	}
	delete [] divec;

}



void UpDateOrdVarAvail(Node *n, int VarI, int left, int right)
{
	int numsplit = right-left+1;
	if(numsplit<1) {
		(n->VarAvail)[VarI]=0;
	} else {
		(n->VarAvail)[VarI]=1;
	}

	if(!(n->Bot)) {
		int lleft,lright,rleft,rright;
		lleft=left;
		rleft=left;
		lright=right;
		rright=right;

		if(((n->rule).Var)==VarI) {
			lright = (n->rule).OrdRule - 1;
			rleft = (n->rule).OrdRule + 1;
		}

		UpDateOrdVarAvail(n->LeftC,VarI,lleft,lright);
		UpDateOrdVarAvail(n->RightC,VarI,rleft,rright);
	}
}


void UpDateCatVarAvail(Node *n, int VarI, int *cats)
{

	int i;
	int RN = RuleNum[VarI];

	if(ISum(RN,cats)<2) {
		(n->VarAvail)[VarI]=0;
	} else {
		(n->VarAvail)[VarI]=1;
	}

	if(!(n->Bot)) {
	
		int *catsl = new int [RN+1];
		int *catsr = new int [RN+1];
		for(i=1;i<=RN;i++) {
			catsl[i]=cats[i];
			catsr[i]=cats[i];
		}

		if(((n->rule).Var)==VarI) {
			for(i=1;i<=RN;i++) {
				if(cats[i]) {
					if((n->rule).CatRule[i]) {
						catsl[i]=0;
					} else {
						catsr[i]=0;
					}
				}
			}
		}

		UpDateCatVarAvail(n->LeftC,VarI,catsl);
		UpDateCatVarAvail(n->RightC,VarI,catsr);
	}

	delete [] cats;
}


void UpDateVarAvail(Node *n,int VarI)
{

	if(VarType[VarI]==CAT) {
		int *cats = new int [RuleNum[VarI]+1];
		GetSetCats(n,VarI,cats);
		UpDateCatVarAvail(n,VarI,cats);// note cats is deleted in here
	} else {
		int LeftI,RightI;
		GetSplitInterval(&LeftI,&RightI, n,VarI);
		UpDateOrdVarAvail(n,VarI,LeftI,RightI);
	}


}

void CheckTree(Node *top)
{
}

void RestrictSize(Node **n,int min)

// this isn't done yet, but I want to modify so it will trim back a tree

// if the tree has fewer than min obs in a node

{

/*	if ( (!(*n)->Bot) & (((*n)->Datalist).length>min)){

		RestrictSize(&((*n)->LeftC),min);

		RestrictSize(&((*n)->RightC),min);

	} 

	

		Rprintf("%d\n",(*n)->DataList.length);

	else {

		RestrictSize(&((*n)->LeftC),min);

		RestrictSize(&((*n)->RightC),min);

		Rprintf("%d\n",(*n)->DataList.length);

	}

	*/

}

double Metrop(Node **top,int *Done,int *step)
{
	double alpha;
	int BD;
	//disable backout----------------
	////Node *copy = new Node;
	//(*top)->CopyTree(copy);
	//++++++++++++++++++++++++++++++++

	//double u = ran1(&idum);
	double u = unif_rand();
	if(u<pBD) {
		alpha = BirthDeath(*top,&BD,Done);
		if(BD) {
			*step = BIRTH;
		} else {
			*step = DEATH;
		}
	} else if(u<pBD+pSwap) { 
		alpha = SwapRule(*top,Done);
		*step=SWAP;
	} else {
		alpha = ChangeRule(*top,Done);
		*step = CHANGE;
	}

	//disable backout
	
	// back out if any nodes of the resultant tree have less than 5 obs---------------------------
	/*
	NodeP *botvec;
	int NBot;
	MakeBotVec(*top,&botvec,&NBot);
	int i;
	int minN = NumObs;
	for (i=1;i<=NBot;i++){
		if ((((botvec[i])->DataList).length) < minN) 
			minN = (((botvec[i])->DataList).length);
	}
	delete [] botvec;

	if (minN < 10){
		(*top)->deall();
		delete (*top);
		(*top) = copy;
		alpha = -1;
	} else {
		copy->deall();
		delete copy;
	}
	*/
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	return alpha;
}
void countVarUsage(std::vector<Node*>& trees, std::vector<int>& cnt)
{
   std::vector<VarUsage> vu;
   cnt.clear(); cnt.resize(NumX+1);
   for(std::vector<Node*>::size_type i=1;i<trees.size();i++) {
      vu.clear();
      getVarUsage(trees[i],0,0,vu);
      for(std::vector<VarUsage>::size_type j=0;j<vu.size();j++) cnt[vu[j].varIndex]++; 
   }   
}
