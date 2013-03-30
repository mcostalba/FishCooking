#pragma once
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#define PI 3.14159265358979323846	// GCC doesn't define this

// Bayeselo formula to get P(win) and P(Loss) as a function of elo, for a given drawelo
extern void proba_elo(double elo, double draw_elo, double *pwin, double *ploss);

// CDF and quantile function of the N(0,1) law
extern double Phi(double x);
extern double Phi_inv(double x);

// Pseudo Random Number Generators
extern uint64_t rand64();	// draw integers (uniformly) between 0 and 2^64-1
extern double uniform();	// draw U(0,1)
extern double gauss();		// draw N(0,1)

// Scores are between -1 and 1 (loss=-1, draw=0, win=1)
// So the elo<->score formulas are not the usual ones
double elo_to_score(double elo);
double score_to_elo(double score);