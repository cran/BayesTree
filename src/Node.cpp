#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "MuS.h"

#include "global.h"
#include "List.h"
//#include "EndNodeModel.h"

#include "Rlob.h"
#include "Lib.h"

#include "Node.h"
#include "Funs.h"

#include "Prior.h"

typedef double *dp;

double Rule::SplitVal()
{
	if(Var==0) return -1000.0;
	if(VarType[Var]!=ORD) return -2000.0;
	return RuleMat[Var][OrdRule];
}


Rule::Rule(Rule &rule) {
	Var=rule.Var;
	if(VarType[Var]==ORD) {
		OrdRule = rule.OrdRule;
	} else {
		CatRule = new int [RuleNum[Var]+1];
		int i;
		for(i=1;i<=RuleNum[Var];i++) CatRule[i] = rule.CatRule[i];
	}
}

Rule::Rule()
{
	Var=0;
	OrdRule=0;
	CatRule=0;
	
}

Rule::~Rule()
{
	
	if(Var && (VarType[Var] == CAT)) {
		
		delete [] CatRule;
		
	
	}
}

void Rule::deall()
{
	if(Var && (VarType[Var] == CAT)) {
		
		delete [] CatRule;
	
	}

	Var=0;
	OrdRule=0;
	CatRule=0;
	
}

int Rule::Right(double *x)
// returns 1 if vector x "goes" to right node, 0 else
{
	int i;
	
	if(VarType[Var]==CAT) {

		for(i=1;i<=RuleNum[Var];i++) {
			// find which of the values x[Var] is then return CatRule[i]
			if((x[Var] == RuleMat[Var][i])) {
				if(CatRule[i]) {
					return 1;
				} else {
					return 0;
				}
			}
		}
		printf("error in Rule::Right, no match found for cat varialble\n");
		return 0;
		
	
	} else {
		if(x[Var]>RuleMat[Var][OrdRule]) {
			return 1;
		} else {
			return 0;
		}
	}
}


void CopyRule(Rule *r1,Rule *r2)
//rule of r1 copied to r2
{
	int i;
	int RN;

	if(r2->Var) r2->deall();
	if(r1->Var) {
		r2->Var = r1->Var;
		if(VarType[r1->Var]==ORD) {
			r2->OrdRule = r1->OrdRule;
		} else {
			RN = RuleNum[r1->Var];
			r2->CatRule = new int[RN + 1];
			for(i=1;i<=RN;i++) (r2->CatRule)[i] = (r1->CatRule)[i];
		}
	}
}



void Node::CopyTree(Node *copy)
{

	int i;
	
	copy->Top = Top;
	copy->Bot = Bot;
	copy->Nog = Nog;

	for(i=1;i<=NumX;i++) copy->VarAvail[i] = VarAvail[i];

	if(!Bot) {
		CopyRule(&rule,&(copy->rule));
		Node *Left,*Right;
		Left = new Node;
		Right = new Node;
		copy->LeftC = Left;
		copy->RightC = Right;
		LeftC->CopyTree(Left);
		RightC->CopyTree(Right);
		Left->Parent = copy;
		Right->Parent = copy;
	}
	if (Top) {
		copy->SetData();
	}
}


void Node::deall() 
{
	if(!Bot) {
		LeftC->deall();
		RightC->deall();
		delete LeftC;
		delete RightC;
		rule.deall();
		Bot=1;
		Nog=0;
		if(!Top) {
			Node *brother = Brother(this);
			if(brother->Bot) (Parent->Nog)=1;
		}

	}
	if(Top) {
		Bot=1;
		Nog=0;
		if(DataList.length) DataList.deall();
		rule.deall();
		int i;
		for (i=1;i<=NumX;i++)
				VarAvail[i]=1;
	}
		
}



void Node::ClearData()
{
	DataList.deall();
	if(!Bot) {
		LeftC->ClearData();
		RightC->ClearData();
	}
}

Node::Node()
{
	Top = 1;
	Bot = 1;
	Nog = 0;

	VarAvail = new int [NumX+1];
	int i;
	for(i=1;i<=NumX;i++) VarAvail[i]=1;
}

Node::~Node()
{


	delete [] VarAvail;
	if(DataList.length) DataList.deall();

	
}



void Node::PrintTree()
{
	int i;
	int d;

	d=Depth(this);

	for(i=1;i<=d;i++) printf(" ");

	//printf("node:%X",this);
	std::cout << "node:" << this;
	printf(" n:%d",DataList.length);
	printf(" TBN: %d%d%d",Top,Bot,Nog);
	printf(" Avail:");
	for(i=1;i<=NumX;i++) printf("%d",VarAvail[i]);
	if(!Top) std::cout << "parent:" << Parent << " ";
	if(!Bot) {
		printf(" Var:%d ",rule.Var);
		if(VarType[rule.Var]==CAT) {
			printf("CATRule: ");
			for(i=1;i<=RuleNum[rule.Var];i++) printf(" %d",rule.CatRule[i]);
		} else {
			printf("ORDRule:");
			printf("(%d)=%f",rule.OrdRule, rule.SplitVal());
		}
	}
	printf("\n");

	if(!Bot) {
		LeftC->PrintTree();
		RightC->PrintTree();
	}

}

void Node::PrintTree(FILE *out)
{
	int i;
	int d;

	d=Depth(this);

	for(i=1;i<=d;i++) fprintf(out," ");


	std::cout << "node:" << this;
	fprintf(out,"node:%X",1);
	fprintf(out," n:%d",DataList.length);
	fprintf(out," TBN: %d%d%d",Top,Bot,Nog);
	fprintf(out," Avail:");
	for(i=1;i<=NumX;i++) fprintf(out,"%d",VarAvail[i]);
	if(!Top) fprintf(out,"parent:%X ",2);
	if(!Bot) {
		fprintf(out," Var:%d ",rule.Var);
		if(VarType[rule.Var]==CAT) {
			fprintf(out,"CATRule: ");
			for(i=1;i<=RuleNum[rule.Var];i++) fprintf(out," %d",rule.CatRule[i]);
		} else {
			fprintf(out,"ORDRule:");
			fprintf(out,"(%d)=%f",rule.OrdRule,rule.SplitVal());
		}
	}
	fprintf(out,"\n");

	if(!Bot) {
		LeftC->PrintTree(out);
		RightC->PrintTree(out);
	}

}

void Node::PrintBernTree(FILE *out)
{
int i;
	int d;

	d=Depth(this);

	for(i=1;i<=d;i++) fprintf(out," ");

	fprintf(out,"node:%X",1);
	fprintf(out,"/%d",DataList.length);
	fprintf(out," TBN:%d%d%d",Top,Bot,Nog);
	fprintf(out," Avail:");
	for(i=1;i<=NumX;i++) fprintf(out,"%d",VarAvail[i]);
	if(!Top) fprintf(out," parent:%X",2);
	if(!Bot) {
		fprintf(out," Var:%d ",rule.Var);
		if(VarType[rule.Var]==CAT) {
			fprintf(out,"CATRule: ");
			for(i=1;i<=RuleNum[rule.Var];i++) fprintf(out,"%d",rule.CatRule[i]);
		} else {
			fprintf(out,"ORDRule:");
			fprintf(out,"(%d)=%f",rule.OrdRule,rule.SplitVal());
		}
	}
	fprintf(out,"\n");

	if(!Bot) {
		LeftC->PrintBernTree(out);
		RightC->PrintBernTree(out);
	}

}



int Node::NumBotNodes()
// returns number of bottom nodes
{
	if(Bot) {
		return 1;
	} else {
		return (LeftC->NumBotNodes() + RightC->NumBotNodes()); 
	}
}

int Node::DepthBelow()
{
	if(Nog) return 1;
	if(Bot) return 0;
	return (1+(int)max(LeftC->DepthBelow(),RightC->DepthBelow()));
}

int Node::NumNogNodes()
{
	if(Bot) return 0;
	if(Nog) return 1;
	return (LeftC->NumNogNodes() + RightC->NumNogNodes());
}

void Node::GetBotList(List **list)
// gets list of pointers to bottom nodes
{
	if(Bot) {	
		*list=new List;
		(**list).length=1;
		Cell *tempcell;
		tempcell=new Cell;
		(*tempcell).contents=this;
		(*tempcell).Beg=1;
		(*tempcell).End=1;
		(**list).first=tempcell;
		(**list).last=tempcell;
	} else {
		List *llist,*rlist;
		LeftC->GetBotList(&llist);
		RightC->GetBotList(&rlist);
		CombineLists(llist,rlist,list);
	}
}

void Node::GetNogList(List **list)
// gets list of pointers to nog nodes
{
	if(Bot) {
		*list=new List;
		(**list).length=0;
	} else {

		if(Nog) {	
			*list=new List;
			(**list).length=1;
			Cell *tempcell;
			tempcell=new Cell;
			(*tempcell).contents=this;
			(*tempcell).Beg=1;
			(*tempcell).End=1;
			(**list).first=tempcell;
			(**list).last=tempcell;
		} else {
			List *llist,*rlist;
			LeftC->GetNogList(&llist);
			RightC->GetNogList(&rlist);
			CombineLists(llist,rlist,list);
		}
	}
}



void Node::GetNotBotList(List **list)
// gets list of pointers to nog nodes
{

	if(Bot) {
		*list=new List;
		(**list).length=0;
	} else if(Nog) {
			
		*list=new List;
		(*list)->length=0;
		AddCellToEnd(*list,(void *)this);

	} else {
		List *llist,*rlist;
		LeftC->GetNotBotList(&llist);
		RightC->GetNotBotList(&rlist);
		CombineLists(llist,rlist,list);
		AddCellToEnd(*list,(void *)this);

	}
	
}

void Node::GetSwapsList(List **list)
{
	if(Nog || Bot) {
		*list=new List;
		(**list).length=0;
	} else if((LeftC->Bot || LeftC->Nog) && (RightC->Bot || RightC->Nog)) {
		*list=new List;
		(*list)->length=0;
		AddCellToEnd(*list,(void *)this);
	} else {

		List *llist,*rlist;
		LeftC->GetSwapsList(&llist);
		RightC->GetSwapsList(&rlist);
		CombineLists(llist,rlist,list);
		AddCellToEnd(*list,(void *)this);
	}
}




void Node::FindNode(double *x,Node **n)
// gets pointer to bottom node for x
{
	
	if(Bot) {
		*n=this;
	}
	else {
		if(rule.Right(x)) {
			RightC->FindNode(x,n);
			
		}
		else {
			LeftC->FindNode(x,n);
		}
	}

}

voidP* Node::GetBotArray()
{
	voidP* botArray=0;

	int nbot = NumBotNodes();

	botArray = new voidP[nbot+1];

	int i;

	List *bots;
	GetBotList(&bots);
	//*NBot = bots->length;
	//*botvec = new NodeP [*NBot+1];
	
		
	Cell *cell = bots->first;
	//(*botvec)[1]=(Node *)cell->contents;
	botArray[1] = (void*)cell->contents;
	
	for(i=2;i<=nbot;i++) {
		cell = cell->after;
		//(*botvec)[i]=(Node *)cell->contents;	
		botArray[i] = (void*)cell->contents;
	}
	
	bots->deall();
	delete bots;

	return botArray;
}

int* Node::GetIndPart(int numObs, double** xx)
{
	int* indPart= new int[numObs+1];
	int i,j;
	//int nbot = NumBotNodes();
	voidP* botvec = GetBotArray();
	Node* nn;
	for(i=1;i<=numObs;i++)
	{
		FindNode(xx[i],&nn);
		for(j=1;((void*)nn) != botvec[j];j++);
		indPart[i] = j;	
	}
	delete [] botvec;
	return indPart;
}

double** Node::GetFits(void* model,int  nTrain,double** xTrain, double** xTrainR, double* yTrain, int nTest, double** xTest, double** xTestR, double* w)
{
	int i,j;
	EndNodeModel* mod = (EndNodeModel*)model;
	double** fits=0;

	//int edim = mod->getEstimateDim();

	int* indPartTrain = GetIndPart(nTrain,xTrain);
	int* indPartTest = GetIndPart(nTest,xTest);

	//for(i=1;i<=nTest;i++) printf("Test: %d %d\n",i,indPartTest[i]);
	//for(i=1;i<=nTrain;i++) printf("Train: %d %d\n",i,indPartTrain[i]);

	fits = new dp[3];
	fits[1] = new double[nTrain+1];
	fits[2] = new double[nTest+1];

	int nbot = NumBotNodes();
	
	int* indObsTrain=0;
	int* indObsTest=0;
	int nobTrain=0,nobTest=0;
	double* tempFits=0;
	int count;

	//double* tempCoef=0;

	for(i=1;i<=nbot;i++)
	{
		nobTrain=0;
		for(j=1;j<=nTrain;j++) {if(indPartTrain[j]==i) nobTrain += 1;}

		nobTest=0;
		for(j=1;j<=nTest;j++) {if(indPartTest[j]==i) nobTest +=1;}

		indObsTrain = new int[nobTrain+1];
		count=0;
		for(j=1;j<=nTrain;j++) {if(indPartTrain[j]==i) {count +=1; indObsTrain[count]=j;}}

		indObsTest = new int[nobTest+1];
		count=0;
		for(j=1;j<=nTest;j++) {if(indPartTest[j]==i) {count +=1; indObsTest[count]=j;}}

		mod->setData(nobTrain,xTrainR,yTrain,indObsTrain,w);

		tempFits = mod->getFits(nobTrain,xTrainR,indObsTrain);
		for(j=1;j<=nobTrain;j++) fits[1][indObsTrain[j]] = tempFits[j];
		if(tempFits!=0) delete [] tempFits;

		tempFits = mod->getFits(nobTest,xTestR,indObsTest);
		for(j=1;j<=nobTest;j++) fits[2][indObsTest[j]] = tempFits[j];
		if(tempFits!=0) delete [] tempFits;

		delete [] indObsTrain;
		delete [] indObsTest;
		

	}

	delete [] indPartTrain;
	delete [] indPartTest;

	return fits;
}


void Node::currentFits(MuS* mod,int nTrain,double** xTrain,double* yTrain,int nTest,double** xTest,double* w, double **fits)
{
        double ybar,postmu,postsd,b,a; //posterior of mu in a bottom node
        double nodeMu; //draw of mu, for a bottom node

        voidP* botvec = GetBotArray(); //bottom nodes
	int* indPartTest;
	if(nTest) indPartTest = GetIndPart(nTest,xTest); //partition of test x re bottom nodes

	int nbot = NumBotNodes();
	int nobTrain=0;
        int *itr;

	for(int i=1;i<=nbot;i++) { // loop over bottom nodes-------------
                //data is list of indices of train obs in the bottom node
                List& data = ((Node *)botvec[i])->DataList;
                nobTrain = data.length;
                itr = new int[nobTrain+1]; //copy list contents to itr

                Cell *cell = data.first;
                if(nobTrain>0) itr[1]=*((int *)(cell->contents));
                ybar = yTrain[itr[1]];
                for(int j=2;j<=nobTrain;j++) {
                   cell = cell->after;
                   itr[j]=*((int *)(cell->contents));
                   ybar += yTrain[itr[j]];
                }
                ybar /= nobTrain;

                b=nobTrain/mod->getSigma2();a=mod->getA();
                postmu = b*ybar/(a+b); postsd = 1.0/sqrt(a+b);
                nodeMu = postmu + postsd*norm_rand();

		for(int j=1;j<=nTest;j++) {if(indPartTest[j]==i) fits[2][j]=nodeMu; }
		for(int j=1;j<=nobTrain;j++) fits[1][itr[j]] = nodeMu;

                delete [] itr;
	} //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if(nTest) delete [] indPartTest;
        delete [] botvec;
}
double** Node::GetEstimates(void* model,int  nTrain,double** xTrain, double** xTrainR, double* yTrain, double* w)
{
	EndNodeModel* mod = (EndNodeModel*)model;
	
	int edim = mod->getEstimateDim();
	int nbot = NumBotNodes();
	
	double** estimates = Lib::almat(nbot,edim);

	int* indPartTrain = GetIndPart(nTrain,xTrain);
	int i,j;
	int* indObsTrain=0;
	
	int nobTrain=0;
	int count;
	double* tempCoef=0;

	for(i=1;i<=nbot;i++)
	{
		nobTrain=0;
		for(j=1;j<=nTrain;j++) if(indPartTrain[j]==i) nobTrain += 1;

		indObsTrain = new int[nobTrain+1];
		count=0;
		for(j=1;j<=nTrain;j++) if(indPartTrain[j]==i) {count +=1; indObsTrain[count]=j;}

		mod->setData(nobTrain,xTrainR,yTrain,indObsTrain,w);

		tempCoef = mod->getParameterEstimate();
		for(j=1;j<=edim;j++) estimates[i][j] = tempCoef[j];
		if(tempCoef!=0) delete [] tempCoef;

		delete [] indObsTrain;
	}

	delete [] indPartTrain;

	return estimates;
}


void Node::SetData()
{
	int i;
	for(i=1;i<=NumObs;i++) SetData(i);
}

void Node::SetData(int i)
{

	
	Cell *cell;
	
	cell=new Cell;
	cell->contents= (&Ivec[i]);
	cell->End=1;
	
	

	if(DataList.length==0) {
	
		DataList.first=cell;
		DataList.last=cell;
					
		cell->Beg=1;
		
			
		
	} else {
			
			
			DataList.last->End=0;
			DataList.last->after=cell;
			cell->before=DataList.last;
			DataList.last=cell;
			
			cell->Beg=0;
			
	}
	
	DataList.length +=1;
			
	if (!Bot){

		if(rule.Right(XDat[i])) {
			RightC->SetData(i);
			
		}
		else {
			LeftC->SetData(i);
		}
		
	
	}
	
	
}
	




int Depth(Node *nn)
{
	int d=0;
	while(!(nn->Top)) {
		d += 1;
		nn = nn->Parent;
	}

	return d;
}

Node *Brother(Node *n)
{
	if(n->Top) return 0;

	if(n==((n->Parent)->LeftC)) {
		return (n->Parent)->RightC;
	} else {
		return (n->Parent)->LeftC;
	}
}







int AddChildsInd(Node *n,int var,int cut)
{
	// Add children to bottom node n
	// use rule for ordinal rule var, with index of rule cut
	// (in global array RuleMat[NumX+1][RuleNum[i]+1]
	// only works for ordinal rules right now

	if(VarType[var]!=ORD)
	{
		printf("AddChilds: error, VarType!=ORD");
		return 0;
	}
	if((n->Bot)!=1)
	{
		printf("AddChilds: error, node not bottom");
		return 0;
	}

	int LeftI,RightI;
	GetSplitInterval(&LeftI,&RightI,n,var);
	if(cut<LeftI) 
	{
		printf("AddChilds: error, cut<LeftI");
		return 0;
	}
	if(cut>RightI)
	{
		printf("AddChilds: error, cut>LeftI");
		return 0;
	}
	
	int LeftEx=0;
	int RightEx=0;

	n->Bot = 0;
	n->Nog = 1;
	if(!(n->Top)) (n->Parent)->Nog=0;

	(n->rule).Var=var;
	(n->rule).OrdRule = cut;

	
	if((n->rule).OrdRule==LeftI) LeftEx=1;
	if((n->rule).OrdRule==RightI) RightEx=1;

	SpawnChildren(n,LeftEx,RightEx);

	return 1;

}

int AddChildsVal(Node *n, int var, double cutVal)
{
	// Add children to bottom node n
	// use rule for ordinal rule var, with value of rule cutVal
	// only works for ordinal rules right now

	// first find index of rule in RuleMat closest to cuVal
	// then call the other AddChilds that uses the index

	int i;
	int index=1;
	double minV=myDoubleAbs(cutVal-RuleMat[var][1]);
	for(i=2;i<=RuleNum[var];i++) 
	{
		if(myDoubleAbs(cutVal-RuleMat[var][i])<minV) 
		{
			minV = myDoubleAbs(cutVal-RuleMat[var][i]);
			index = i;
		}
	}

	int retVal = AddChildsInd(n,var,index);
	return retVal;
}

void getVarUsage(Node* node, int depth, int nodeIndex, std::vector<VarUsage>& vu)
{
   if(!node->Bot) {
      VarUsage temp_vu; temp_vu.depth = depth; temp_vu.nodeIndex = nodeIndex; temp_vu.varIndex = node->rule.Var;
      vu.push_back(temp_vu);
      if(!node->Nog) {
         getVarUsage(node->LeftC,depth+1,2*nodeIndex+1,vu);
         getVarUsage(node->RightC,depth+1,2*nodeIndex+2,vu);
      } 
   }
}

void printVarUsageVector(const std::vector<VarUsage>& vs)
{
   std::cout << "Printing Variable Usage Vector:\n";
   std::cout << "Depth\tNode\tVariable\n";
   for(std::vector<VarUsage>::size_type i=0;i!=vs.size();i++) 
      std::cout << vs[i].depth << "\t" << vs[i].nodeIndex << "\t" << vs[i].varIndex << "\n";
}


