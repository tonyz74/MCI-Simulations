#include "common.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng;
double SIGMAS[DIMS] = { 1.2, 1.2 };
double NC = 0.744903010757;

int rej = 0;

// Performing rejection sampling like this is quite inefficient
// in regions where the valid values for `r` are in a low-probability
// region. However, it does provide accurate results.
//
// Extensions to this EE might consider how Monte Carlo methods
// (random simulations) can also sample from arbitrary distributions,
// a common method is the Metropolis-Hastings algorithm.
point sampler() {
    point res;
    for (int k = 0; k < DIMS; k++) {
        double r;
        rej--;
        do {
            r = gsl_ran_gaussian(rng, SIGMAS[k]);
            rej++;
        } while (r < LOWER_BOUND.x[k] || r > UPPER_BOUND.x[k]);

        res.x[k] = r;
    }
    
    return res;
}

double pdf(point p) {
    double ans = 1.0;

    for (int k = 0; k < DIMS; k++) {
        ans *= gsl_ran_gaussian_pdf(p.x[k], SIGMAS[k]);
    }
    
    return ans / NC; 
}

double imp_integrate() {
    double ans = 0;
    int N = N_TOTAL_SAMPS;

    for (int i = 0; i < N; i++) {
        point samp = sampler(); 
        double p = pdf(samp);
        double fval = f(samp);

        ans += (1.0 / N) * (fval / p);
    }

    return ans;
}



int main(void) {
    srand(77);
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 77);

    double trials[N_TRIALS] = {0};
    for (int i = 0; i < N_TRIALS; i++) {
        int cnt = i + 1;
        if (cnt % 200 == 0) {
            printf("%.2lf%%\n", 100.0 * (double) cnt / N_TRIALS);
        }

        trials[i] = imp_integrate();
    }

    print_trials_summary(trials);
    write_trials("output/imp_trials.txt", trials);

    printf("acceptance rate = %lf\n", (double) N_TOTAL_SAMPS / ( rej + N_TOTAL_SAMPS));
}

