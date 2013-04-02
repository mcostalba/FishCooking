#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "stat.h"

/* number of games for each simulation */
#define T	16000

/* DRAW_ELO controls the proportion of draws. See proba_elo() function. To estimate this value, use
 * draw_elo = 200.log10[(1-w)/w.(1-l)/l]
 * where (w,l) are the win and loss ratio */
#define DRAW_ELO	240

/* Parametrization of the SPRT is here */
const double elo0 = 0.0, elo1 = 5.0; 	     // expressed in BayesELO units
const double alpha = 0.05; //0.05;			// alpha = max type I error when elo < elo1
const double beta = 0.05; //1.0  			// alpha = max type II error when elo > elo2

typedef struct {
	unsigned t;			// time (in nb of games)
	double threshold;	// stop when the random walk is below that value at time t
} Stop;

void discrete_stop(double pwin, double ploss, unsigned nb_simu,
	double *typeI, double *typeII, unsigned *avg_stop)
{
	// Mean and Variance of Xt
	double mu = pwin-ploss, v = pwin+ploss - mu*mu;
	
	// Final stop at LOS < 95%
	double final_stop = sqrt(v*T) * Phi_inv(0.95);
	
	/* Stopping rule defined here */
    /*
	#define NB_STEP	1
	const Stop S[1+NB_STEP] = {
		{0, 0},				// initial time point (threshold discarded)
		{T, final_stop}		// LOS < 95% stop after T games
	};
    */

	#define NB_STEP	4
	const Stop S[1+NB_STEP] = {
		{0, 0},				// initial time point (threshold discarded)
        {4000, 4000 * elo_to_score(-5)},
        {8000, 8000 * elo_to_score(0)},
        {12000, 12000 * elo_to_score(0)},
		{T, final_stop}		// LOS < 95% stop after T games
	};
	
		
	// Counter for type I and type II error
	unsigned typeI_cnt = 0, typeII_cnt = 0;
	
	// Sum of stopping times
	uint64_t sum_stop = 0;
	
	// Simulation loop
	for (unsigned simu = 0; simu < nb_simu; ++simu) {
		// We simulate the random walk at the stopping points directly. This approximation is OK if
		// you look at the random walk from "far enough". It allows to reduce the computation time
		// massively (compared to generating the random walk step by step, generating one game at a
		// time)
		double W = 0;				// random walk Wt = X1 + ... + Xt (where Xi are game results)
		double rejected = false;	// set to true when a stop rejects the patch
		unsigned i;
		
		for (i = 1; i <= NB_STEP; ++i) {
			const unsigned dt = S[i].t - S[i-1].t;	// dt = t'-t
			W += mu*dt + sqrt(v*dt)*gauss();		// use the law of Wt'-Wt ~= N(mu, v.dt)
			
			if (W < S[i].threshold) {	// apply the i-th stopping rule
				rejected = true;
				sum_stop += S[i].t;
				break;
			} else if (i == NB_STEP)
				sum_stop += S[i].t;
		}
		
		typeI_cnt += (mu <= 0) && !rejected;	// type I error = false positive (the riskiest one)
		typeII_cnt += (mu > 0) && rejected;		// typeII error = false non positive
	}
	
	*typeI = (double)typeI_cnt / nb_simu;
	*typeII = (double)typeII_cnt / nb_simu;
	*avg_stop = sum_stop / nb_simu;
}

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

	// Calculate the true values of E(Xt) and V(Xt)
	const double mu = pwin-ploss, v = pwin+ploss - mu*mu;
	
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
		for (t = 0; /*t < T*/; ++t) {
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
		
		/*if (t == T) {
			// patch was not early stopped
			
			// Calculate the empirical mean and variance E_hat(Xi) and V_hat(Xi)
			// "mu_hat" and "v_hat" are estimators of the unknown "mu" and "v"
			double mu_hat = (double)(count[WIN+1] - count[LOSS+1]) / T;
			double v_hat = (double)(count[WIN+1] + count[LOSS+1])/T - mu_hat * mu_hat;			
			
			// Apply the empirical LOS > 95% condition
			accepted = mu_hat > sqrt(v_hat/T) * Phi_inv(0.95);
			sum_stop += T;
		}*/

		typeI_cnt += (mu <= 0) && accepted;		// type I error = false positive
		typeII_cnt += (mu > 0) && !accepted;	// typeII error = false non positive
	}
	
	*typeI = (double)typeI_cnt / nb_simu;
	*typeII = (double)typeII_cnt / nb_simu;
	*avg_stop = sum_stop / nb_simu;
}

void real_data_test(char *filename)
{
    float values[32000];
    int valueCount = 0;

    // Read data
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;    

    fp = fopen(filename, "r");
    if (!fp)
    {
        puts("Failed to open data file!");
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, fp)) != -1) {
        if (sscanf(line, "%f", &values[valueCount++]) <= 0)
        {
            printf("Parse file for data file on line: %d", valueCount);
            exit(EXIT_FAILURE);
        }

    }

    // Loop through samples, collecting statistics:
    unsigned avg_stop_tot = 0;
    double type1_elo_lost_tot = 0.0;
    double type2_elo_lost_tot = 0.0;

    for (int i = 0; i < valueCount; i++)
    {
        double RealElo = values[i];

        // FIXME: Conversion to bayes-elo
        // HACK : START FORKING.
        double min = -1000.0, max = 1000.0;

        while (max - min > 0.001)
        {
            double tryElo = (max + min) / 2;
            double tryWin, tryLoss;

            proba_elo(tryElo, DRAW_ELO, &tryWin, &tryLoss); 
            double tryScore = 0.5 + (tryWin - tryLoss) / 2;
            double ELO = -400 * log10(1/tryScore - 1);          
 
            if (ELO > RealElo)
                max = tryElo;
            else
                min = tryElo;
        }

        double elo = (max + min) / 2;

        // HACK: END FORKING 
 
        double pwin, ploss, typeI, typeII;
        unsigned avg_stop;

        proba_elo(elo, DRAW_ELO, &pwin, &ploss);
        double score = 0.5 + (pwin - ploss) / 2;
        double ELO = -400 * log10(1/score - 1);

        //discrete_stop(pwin, ploss, 1000, &typeI, &typeII, &avg_stop);
        SPRT_stop(pwin, ploss, 1000, &typeI, &typeII, &avg_stop);

        printf("%.2f,%.4f,%.4f,%.2f,%.4f,%.4f,%u\n", elo, pwin, ploss, ELO, typeI, typeII, avg_stop); 
   
        type1_elo_lost_tot += elo * typeI;
        type2_elo_lost_tot += elo * typeII; 
        avg_stop_tot += avg_stop;
    }

    printf("Average stop: %d\n", avg_stop_tot / valueCount);
    printf("Total type1 elo loss: %1.4f\n", type1_elo_lost_tot);
    printf("Total type2 elo loss: %1.4f\n", type2_elo_lost_tot);
}

int main(int argc, char **argv)
{
    real_data_test("./framework_tests.csv");
    exit(EXIT_SUCCESS);

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
		
		discrete_stop(pwin, ploss, nb_simu, &typeI, &typeII, &avg_stop);
		//SPRT_stop(pwin, ploss, nb_simu, &typeI, &typeII, &avg_stop);
		
		printf("%.2f,%.4f,%.4f,%.2f,%.4f,%.4f,%u\n", elo, pwin, ploss, ELO, typeI, typeII, avg_stop);
	}
}
