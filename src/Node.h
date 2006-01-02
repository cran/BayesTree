#include <stdio.h>
#include <vector>
#include <cmath>

#include "List.h"
#include "MuS.h"

typedef void *voidP;

class Rule {
public:

	Rule();
	~Rule();
	Rule(Rule &rule);

	int Var; // index of variable associated with rule, 0 on creation

	int OrdRule; // if ordered this will be the index of the split value

	int *CatRule;

	int Right(double *x); // returns 1 if vector x "goes" to right node, 0 else
	void deall();

	double SplitVal();
};

void CopyRule(Rule *r1,Rule *r2);
class Node {
public:

	Node();
	~Node();

	// next three are indicators for 
	// whether the node is a top, bottom, or nogrand
	int Top;
	int Bot;
	int Nog;
	
	// pointers for tree structure
	Node *Parent;
	Node *LeftC;
	Node *RightC;

	//the rule for splitting the observations
	Rule rule;
	int *VarAvail; // ith is 1 if variable i has rules left, 0 else

	
	//list of observations corresponding to node
	List DataList;

	//functions
	int NumBotNodes();  // returns number of bottom nodes
	int NumNogNodes(); // returns number of Nog nodes
	void GetBotList(List **list);  // gets list of pointers to bottom nodes
	void GetNogList(List **list);
	void GetNotBotList(List **list);
	void GetSwapsList(List **list);
	void FindNode(double *x,Node **n); // gets pointer to bottom node for x

	voidP* GetBotArray();
	int* GetIndPart(int numObs, double** x);
	double** GetFits(void* model,int  nTrain,double** dTrain, double** dTrainR,double* yTrain, int nTest, double** dTest, double** dTestR, double* w);
	double** GetEstimates(void* model,int  nTrain,double** dTrain, double** dTrainR,double* yTrain, double* w);
        //void  currentFits(MuS* mod,int nTrain,double** xTrain,double* yTrain,int nTest,double** xTest,double* w, Lib& gen, double** fits);
        void  currentFits(MuS* mod,int nTrain,double** xTrain,double* yTrain,int nTest,double** xTest,double* w, double** fits);

	void PrintTree();
	void PrintTree(FILE *out);
	void PrintBernTree(FILE *out);
	void SetData(int i); //goes thru the tree and adds i to list of obs indices for each appropriate node
	void SetData();
	void ClearData(); // deallocate all the DataLists for all nodes in tree
	void deall();
	void CopyTree(Node *copy);

	int DepthBelow();
};

int Depth(Node *n);
Node *Brother(Node *n);

int AddChildsInd(Node *n,int var,int cut);
int AddChildsVal(Node *n,int var,double cutVal);

class VarUsage {
public:
   int depth;
   int nodeIndex;
   int varIndex;
};

void getVarUsage(Node* node, int depth, int nodeIndex, std::vector<VarUsage>& vu);
void printVarUsageVector(const std::vector<VarUsage>& vs);
