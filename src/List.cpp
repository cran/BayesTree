//#include <stdio.h>
#include <cstdio>
#include <iostream>
//#include <math.h>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib>
//#include <time.h>
#include <ctime>

typedef void *voidP;

// extern "C" {
#include <R.h>
#include <Rmath.h>
// };



class Cell {
public:
	int Beg; //flag for is it the first cell(1=yes 0=no)
	int End; //flag for is it the last cell (1=yes 0=no)
	Cell *before;
	Cell *after;
	void *contents;
};


class List {
public:
	Cell *first;
	Cell *last;
	int length;
	void deall();
	List();
	~List();
};

List::~List()
{
	//deall();
}

List::List()
{
	length=0;
}

void List::deall()
{
	if(length>0) {
		Cell *kill=0,*next=0;
		kill=first;
		for(int i=1;i<=length;i++) {
			if(i<length) next=kill->after;
			delete kill;
			kill=next;
		}
		length=0;
	}
	
}


void CombineLists(List *list1,List *list2, List **list)
{
	int n1,n2;

	n1=(*list1).length;
	n2=(*list2).length;
	
	if(n1==0) {
		(*list)=list2;
		delete list1;
		return;
	}
	if(n2==0) {
		(*list)=list1;
		delete list2;
		return;
	}
	if((n1>0) && (n2>0)) {			
	
		(*list)=new List;
	
	
		(**list).length=n1+n2;
		(**list).first=(*list1).first;
		(**list).last=(*list2).last;
	
		((*list1).last)->after=(*list2).first;
		((*list1).last)->End=0;
		((*list2).first)->before=(*list1).last;
		((*list2).first)->Beg=0;
	
		delete list1;
		delete list2;
	}
	
}


void PrintList(List *list)
{

	Cell *cell;
	int n=(*list).length;
	Rprintf("\n the length of the list is %d\n",n);
	
	if(n>0) {
	Rprintf("the pointer contents and status values are:\n");
	cell =(*list).first;
	Rprintf(" %p %d %d\n",cell->contents,cell->Beg,cell->End);
	
	for(int i=1;i<n;i++) {
		cell=(*cell).after;
		Rprintf(" %p %d %d\n",cell->contents,cell->Beg,cell->End);
		
	}
	//std::cout << "\n";
        Rprintf("\n");
	}
}	

void DelCell(List *list,Cell *cell)
{
	(list->length) -=1;
	if((list->length)>0) {
		if(cell->Beg) {
			(cell->after)->Beg=1;
			list->first = (cell->after);
		}
		else if(cell->End) {
			(cell->before)->End=1;
			list->last = cell->before;
		} else {
			(cell->before)->after = cell->after;
			(cell->after)->before = cell->before;
		}
	}

	delete cell;
}

void AddCellToEnd(List *list, void *p)
{
	int len = list->length;

	Cell *cell;
	cell = new Cell;
	cell->contents = p;
	cell->End = 1;
	if(len) {
		(list->last)->End=0;
		(list->last)->after = cell;

		cell->before = list->last;
		cell->Beg=0;
	} else {
		list->first=cell;
		cell->Beg=1;
	}
	list->last = cell;
	list->length +=1;
}

void AddCellToBeg(List *list, void *p)
{
	int len = list->length;

	Cell *cell;
	cell = new Cell;
	cell->contents = p;
	cell->Beg = 1;
	if(len) {
		(list->first)->Beg=0;
		(list->first)->before = cell;

		cell->after = list->first;
		cell->End=0;
	} else {
		list->last=cell;
		cell->End=1;
	}
	list->first = cell;
	list->length +=1;
}

void AddCellAfter(List *list,Cell *oldcell,void *p)
{
	if(oldcell->End) {
		AddCellToEnd(list,p);
	} else {
		Cell *cell;
		cell = new Cell;

		cell->contents = p;
		cell->Beg=0;
		cell->End=0;
		cell->before = oldcell;
		cell->after = oldcell->after;

		oldcell->after = cell;
		(cell->after)->before = cell;

		list->length +=1;
	}
}

void AddCellBefore(List *list,Cell *oldcell,void *p)
{
	if(oldcell->Beg) {
		AddCellToBeg(list,p);
	} else {
		Cell *cell;
		cell = new Cell;

		cell->contents = p;
		cell->Beg=0;
		cell->End=0;
		cell->before = oldcell->before;
		cell->after = oldcell;

		oldcell->before = cell;
		(cell->before)->after = cell;

		list->length +=1;
	}
}

void ListToVector(List *list,voidP **p,int *n)
// allocates and defines array of void pointers corresponding to the list
{

	int i;

	
	*n = list->length;
	*p = new voidP [*n+1];

	if((*n)) {
	
		
		Cell *cell = list->first;
		(*p)[1]=cell->contents;
	
		for(i=2;i<=(*n);i++) {
			cell = cell->after;
			(*p)[i]=cell->contents;	
		}
	}
	
	
}













