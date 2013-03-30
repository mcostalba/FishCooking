#include "stat.h"

void proba_elo(double elo, double draw_elo, double *pwin, double *ploss)
{
	*pwin  = 1 / (1 + pow(10, (-elo + draw_elo) / 400));
	*ploss = 1 / (1 + pow(10, (elo + draw_elo) / 400));
}

double elo_to_score(double elo)
{
	return 2/(1+pow(10,-elo/400)) - 1;
}

double score_to_elo(double score)
{
	assert(0 < score && score < 1);
	return -400*log10(2/(1+score) - 1);
}

uint64_t rotate(uint64_t x, uint64_t k)
{ return (x << k) | (x >> (64 - k)); }

// RKISS 64-bit PRNG by Bob Jenkins
uint64_t rand64()
{
	static uint64_t		// seeds by Heinz Van Saanen (dieharder-proof)
		a = 0x46dd577ff603b540,
		b = 0xc4077bddfacf987b,
		c = 0xbbf4d93b7200e858,
		d = 0xd3e075cfd449bb1e;
		
	const uint64_t e = a - rotate(b,  7);
	a = b ^ rotate(c, 13);
	b = c + rotate(d, 37);
	c = d + e;
	
	return d = e + a;
}

static double erf_inv(double x)
{
	static const double a = 8*(PI-3)/(3*PI*(4-PI));
	const double y = log(1-x*x), z = 2/(PI*a) + y/2;
	return (x > 0 ? 1 : -1) * sqrt(sqrt(z*z -y/a) - z);
}

double Phi(double x)	 { return 0.5 + 0.5*erf(x/sqrt(2)); }
double Phi_inv(double x) { return sqrt(2)*erf_inv(2*x-1); }

double uniform() { return (double)rand64() / 0xffffffffffffffffULL; }
double gauss()	 { return Phi_inv(uniform());}