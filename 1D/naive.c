#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double naive_integrate(integrand fn, double lb, double ub, int n_samps) {
    double len = ub - lb;

    double estimate = 0;
    for (int i = 0; i < n_samps; i++) {
        double x = randf_range(lb, ub);
        estimate += (len * fn(x)) / n_samps;
    }

    return estimate;
}

int main(void) {
    srand(time(NULL));
    
    double trials[N_TRIALS];
    for (int i = 0; i < N_TRIALS; i++) {
        trials[i] = naive_integrate(
            f, LOWER_BOUND, UPPER_BOUND,
            N_TOTAL_SAMPS
        );
    }

    write_trials("output/naive_trials.txt", trials);
}
