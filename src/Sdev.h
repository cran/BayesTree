#ifndef GUARD_Sdev
#define GUARD_Sdev

#include "Lib.h"

class Sdev
{
public:

	// cons dest
	Sdev(): s(1.0),nu(3),lambda(1.0),nob(0),e(0) {}
	~Sdev() {}
	// get set
	double getS() { return this->s;}
	void setS(double s) { this->s = std::abs(s);}
	// public methods
	void setPrior(int nu, double lambda) {this->nu = nu; this->lambda = lambda;}
	void setData(int nob, double* e) { this->nob = nob; this->e = e;}
	void drawPost()
	{
		double ss = 0.0;
		for(int i=1;i<=nob;i++) ss += e[i]*e[i];
		int nupost(nu+nob);
		double nlpost(nu*lambda + ss);
		//s = sqrt(nlpost/(Lib::genChi2(nupost)));
		s = sqrt(nlpost/rchisq((double)nupost));
	}
private:
	//parameter
	double s;
	// prior
	int nu;
	double lambda;
	// data
	int nob;
	double* e;
};
#endif


