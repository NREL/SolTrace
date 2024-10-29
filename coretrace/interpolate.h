#ifndef __interpolate_h
#define __interpolate_h

#include <vector>
typedef std::vector<double> VectDoub;
typedef std::vector<VectDoub >  MatDoub;


struct LUdcmp
{
	/*
	LU Decomposition (solution to matrix equation A . x = b)
	*/

	int n;

	MatDoub lu;
	MatDoub aref;

	std::vector<int> indx;
	double d;

	LUdcmp(MatDoub& a);
	void solve(VectDoub& b, VectDoub& x);
	void solve(MatDoub& b, MatDoub& x);
	void inverse(MatDoub& ainv);
	double det();
	void mprove(VectDoub& b, VectDoub& x);
};

class GaussMarkov 
{
protected:
	double pvg_dist(const double r);
	//Powvargram vgram;
	int ndim, npt;
	double lastval, lasterr;
	VectDoub dstar, vstar, yvi;
	MatDoub v;
	LUdcmp* vi;
	double SQR(const double a);
	double alph, bet, nugsq;

public:
	MatDoub x;	//
	VectDoub y;

	GaussMarkov();
	~GaussMarkov();
	
	void setup(const double beta = 1.99, const double nug = 0., const double* err = NULL);

	double interp(VectDoub& xstar);

	double interp(VectDoub& xstar, double& esterr);

	double rdist(VectDoub* x1, VectDoub* x2);
};





#endif
