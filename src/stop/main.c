#include <stdio.h>
#include <stdlib.h>
#include "stat.h"

typedef struct {
	double t, threshold;
} Stop;

void foo(double elo, unsigned nb_simu)
{
	const unsigned T = 16000;
	
	double pwin, ploss;
	proba_elo(elo, 200.0, &pwin, &ploss);
	
	// Mean and Variance of Xt
	double mu = pwin-ploss, v = pwin+ploss - mu*mu;
	
	// Final stop at LOS>95%
	double final_stop = sqrt(v*T) * Phi_inv(0.95);
	
	// Stopping rules
	#define NB_STEP	2
	const Stop S[1+NB_STEP] = {
		{0, 0},						// initial time point (threshold discarded)
		{T/2, 0},		// elo < 0 after T/2 games
		{T, final_stop}				// LOS < 95% stop after T games
	};
	
	// Counter for type I and type II error
	unsigned typeI = 0, typeII = 0;
	
	// Simulation loop
	for (unsigned simu = 0; simu < nb_simu; ++simu) {
		// We simulate the random walk at the stopping points directly. This approximation is OK if
		// you look at the random walk from "far enough". It allows to reduce the computation time
		// massively (compared to generating the random walk step by step, generating one game at a
		// time)
		double W = 0;				// random walk Wt = X1 + ... + Xt (where Xi are game results)
		double rejected = false;	// set to true when a stop rejects the patch
		
		for (unsigned i = 1; i <= NB_STEP; ++i) {
			const double dt = S[i].t - S[i-1].t;	// dt = t'-t
			W += mu*dt + sqrt(v*dt)*gauss();		// use the law of Wt'-Wt ~= N(mu, v.dt)
			
			if (W < S[i].threshold) {				// apply the i-th stopping rule
				rejected = true;
				break;
			}
		}
		
		typeI += (mu <= 0) && !rejected;	// type I error = false positive (the riskiest one)
		typeII += (mu > 0) && rejected;		// typeII error = false non positive
	}
	
	printf("elo: %1.2f\ttypeI: %1.4f\ttypeII: %1.4f\n",
		elo, (double)typeI / nb_simu, (double)typeII/ nb_simu);
}

int main(int argc, char **argv)
{
	if (argc != 5) {
		puts("4 parameters requires: elo_min elo_max elo_step nb_simu\n");
		exit(EXIT_FAILURE);
	}
	
	const double elo_min = atof(argv[1]), elo_max = atof(argv[2]), elo_step = atof(argv[3]);
	const unsigned nb_simu = atoll(argv[4]);
	
	for (double elo = elo_min; elo <= elo_max; elo += elo_step)
		foo(elo, nb_simu);
}