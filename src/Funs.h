typedef Node *NodeP;

void MakeBotVec(Node *top,NodeP **botvec,int *NBot);
void MakeNogVec(Node *top,NodeP **nogvec,int *NNog);
void MakeNotBotVec(Node *top,NodeP **notbotvec,int *Nnotbot);
void MakeSwapVec(Node *top,NodeP **swapvec,int *Nswap);


void MakeIntVec(List *intlist, int **ivec, int *n);


void GetDataInd(Node *top,int *ind);
void GetDataInd(Node *top,int *ind, int NumObsPred, double** data);

int getMaxDepth(Node *top);



//double FitMetric(Node* top1,Node* top2);
//int MisCMetric(Node* top1,Node* top2);
int AndrewsMetric(Node *top1,Node *top2);

Rule *GetRulePointer(int index, int curindex, int depth, int curdepth, Node* n);
int RulesDifferent(Rule *r1,Rule *r2);
int ShannonBanksMetric(Node *top1,Node *top2);



void AddDatChildren(Node *n);
void FixDataBelow(Node *cnode);
void CopyRule(Rule *r1,Rule *r2);

void UpDateOrdVarAvail(Node *n, int VarI, int left, int right);
void UpDateCatVarAvail(Node *n, int VarI, int *cats);
void UpDateVarAvail(Node *n,int VarI);

void CheckTree(Node *top);



void RestrictSize(Node **top,int *min);

double Metrop(Node **top,int *Done,int *step);

void countVarUsage(std::vector<Node*>& trees, std::vector<int>& cnt);





