#include "common.h"
#include <time.h>
#include <stdio.h>

double naive_integrate() {
    point ub = UPPER_BOUND;
    point lb = LOWER_BOUND;
    int N = N_TOTAL_SAMPS;
    integrand fn = f;

    double hypervol = 1.0;
    for (int k = 0; k < DIMS; k++) {
        hypervol *= (ub.x[k] - lb.x[k]);
    }

    double sum_ord1 = 0.0;
    double sum_ord2 = 0.0;

    double m = 0.0;
    for (int i = 0; i < N; i++) {
        point p;
        for (int k = 0; k < DIMS; k++) {
            p.x[k] = randf_range(lb.x[k], ub.x[k]);
        }
        
        double samp = hypervol * fn(p);
        sum_ord1 += samp;
        sum_ord2 += samp * samp;
        m += (1.0 / N) * samp;
    }

    double var = (sum_ord2 - 2.0 * sum_ord1 * m + N * m * m) / (N - 1);
    (void) var;
    // printf("var: %lf\n", var);

    return m;
}

int main(void) {
    srand(77);
    
    double trials[N_TRIALS] = {0};
    for (int i = 0; i < N_TRIALS; i++) {
        int cnt = i + 1;
        if (cnt % 200 == 0) {
            printf("%.2lf%%\n", 100.0 * (double) cnt / N_TRIALS);
        }

        trials[i] = naive_integrate();
    }

    print_trials_summary(trials);
    write_trials("output/naive_trials.txt", trials);
}
