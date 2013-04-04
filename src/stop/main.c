#include <stdio.h>
#include <stdlib.h>
#include "stat.h"

/* Use this value to cap the number of games (standard SPRT has no cap) */
#define T	0xffffffff

/* DRAW_ELO controls the proportion of draws. See proba_elo() function. To estimate this value, use
 * draw_elo = 200.log10[(1-w)/w.(1-l)/l]
 * where (w,l) are the win and loss ratio */
#define DRAW_ELO	240

/* Parametrization of the SPRT is here */
const double elo0 = 0, elo1 = 7;	// expressed in BayesELO units
const double alpha = 0.05;		// alpha = max type I error (reached on elo = elo0)
const double beta = 0.05;			// beta = max type II error for elo >= elo2 (reached on elo = elo2)

void SPRT_stop(double pwin, double ploss, unsigned nb_simu,
	double *typeI, double *typeII, unsigned *avg_stop)
{	
	// LLR bounds
	const double lower_bound = log(beta / (1-alpha));
	const double upper_bound = log((1-beta) / alpha);
	
	// Calculate the probability laws under H0 and H1
	double pwin0, ploss0, pdraw0;
	double pwin1, ploss1, pdraw1;
	proba_elo(elo0, DRAW_ELO, &pwin0, &ploss0); pdraw0 = 1-pwin0-ploss0;
	proba_elo(elo1, DRAW_ELO, &pwin1, &ploss1); pdraw1 = 1-pwin1-ploss1;
	
	// Calculate the log-likelyhood ratio increment for each game result Xt
	const double llr_inc[3] = {
		log(ploss1 / ploss0),
		log(pdraw1 / pdraw0),
		log(pwin1 / pwin0)
	};

	// Calculate the true values of E(Xt) and V(Xt), and elo
	const double mu = pwin-ploss, v = pwin+ploss - mu*mu;
	const double elo = 200*log10(pwin/ploss*(1-ploss)/(1-pwin));
	
	// Collect the risk and reward statistics along the way
	unsigned typeI_cnt = 0, typeII_cnt = 0;	// Counter for type I and type II error
	uint64_t sum_stop = 0;					// sum of stopping times (to compute average)

	for (unsigned simu = 0; simu < nb_simu; ++simu) {
		/* Simulate one trajectory of T games
		 * - Calculate the LLR random walk along the way
		 * - early stop when LLR crosses a bound */
		bool accepted;					// patch acceptation
		double LLR = 0;					// log-likelyhood ratio
		unsigned count[3] = {0,0,0};	// counts the number of: LOSS, DRAW, WIN (in this order)
		
		unsigned t;
		for (t = 0; t < T; ++t) {
			const int X = game_result(pwin, ploss);
			LLR += llr_inc[X+1];
			++count[X+1];
			
			if (LLR < lower_bound) {
				// patch rejected by early stopping
				accepted = false;
				sum_stop += t;
				break;
			} else if (LLR > upper_bound) {
				// patch accepted by early stopping
				accepted = true;
				sum_stop += t;
				break;
			}
		}
		
		if (t == T) {
			// patch was not early stopped
			
			// Calculate the empirical mean and variance E_hat(Xi) and V_hat(Xi)
			// "mu_hat" and "v_hat" are estimators of the unknown "mu" and "v"
			double mu_hat = (double)(count[WIN+1] - count[LOSS+1]) / T;
			double v_hat = (double)(count[WIN+1] + count[LOSS+1])/T - mu_hat * mu_hat;			
			
			// Apply the empirical p-value > 95% condition
			accepted = mu_hat > sqrt(v_hat/T) * Phi_inv(0.95);
			sum_stop += T;
		}

		typeI_cnt += (elo <= elo0) && accepted;		// type I error
		typeII_cnt += (elo > elo0) && !accepted;	// typeII error
	}
	
	*typeI = (double)typeI_cnt / nb_simu;
	*typeII = (double)typeII_cnt / nb_simu;
	*avg_stop = sum_stop / nb_simu;
}

int main(int argc, char **argv)
{
	if (argc != 5) {
		puts("4 parameters requires: elo_min elo_max elo_step nb_simu\n");
		exit(EXIT_FAILURE);
	}
	
	const double elo_min = atof(argv[1]), elo_max = atof(argv[2]), elo_step = atof(argv[3]);
	const unsigned nb_simu = atoll(argv[4]);
	
	// Print header
	puts("BayesELO,P(win),P(loss),ELO,typeI,typeII,avg(stop)");

	for (double elo = elo_min; elo <= elo_max; elo += elo_step) {
		double pwin, ploss, typeI, typeII;
		unsigned avg_stop;
		
		proba_elo(elo, DRAW_ELO, &pwin, &ploss);
		double score = 0.5 + (pwin - ploss) / 2;
		double ELO = -400 * log10(1/score - 1);
		
		SPRT_stop(pwin, ploss, nb_simu, &typeI, &typeII, &avg_stop);
		printf("%.2f,%.4f,%.4f,%.2f,%.4f,%.4f,%u\n", elo, pwin, ploss, ELO, typeI, typeII, avg_stop);
	}
}
