#include "common.h"
#include <stdio.h>
#include <time.h>

#define STRATA 100

typedef struct ss_state {
    integrand func;
    int n_each;
    point pmin, pmax;
} ss_state;

ss_state ss_init(void) {
    ss_state st;
    st.pmin = LOWER_BOUND; st.pmax = UPPER_BOUND;
    st.func = f;

    int regs = 1;
    for (int i = 0; i < DIMS; i++) {
        regs *= STRATA;
    }
    st.n_each = N_TOTAL_SAMPS / regs;
   
    return st;
}

double _ss_integrate_region(
    ss_state *const st,
    point regmin,
    point regmax
) {
    int N = st->n_each; double ans = 0.0;
    double hypervol = 1.0;
    for (int k = 0; k < DIMS; k++) {
        hypervol *= regmax.x[k] - regmin.x[k];
    }

    for (int i = 0; i < N; i++) {
        point x;
        for (int k = 0; k < DIMS; k++) {
            x.x[k] = randf_range(regmin.x[k], regmax.x[k]);
        }

        ans += (1.0 / N) * st->func(x) * hypervol;
    }

    return ans;
}

// Currently enumerating over the k-th dimension
double _ss_integrate_recurse(
    ss_state *const st,
    point regmin, 
    point regmax,
    int k
) {
    if (k == DIMS) {
        return _ss_integrate_region(st, regmin, regmax);
    }

    double step = (st->pmax.x[k] - st->pmin.x[k]) / STRATA;
    double accum = 0;
    for (int h = 0; h < STRATA; h++) {
        point newmin = regmin;
        point newmax = regmax;

        newmin.x[k] = h * step + st->pmin.x[k];
        newmax.x[k] = (h + 1) * step + st->pmin.x[k];

        accum += _ss_integrate_recurse(st, newmin, newmax, k + 1);
    }

    return accum;
}

double ss_integrate(ss_state *const st) {
    return _ss_integrate_recurse(st, (point) { 0 }, (point) { 0 }, 0);
}

int main(void) {
    srand(77);
    
    double trials[N_TRIALS] = {0};
    for (int i = 0; i < N_TRIALS; i++) {
        int cnt = i + 1;
        if (cnt % 200 == 0) {
            printf("%.2lf%%\n", 100.0 * (double) cnt / N_TRIALS);
        }

        ss_state ss = ss_init();
        trials[i] = ss_integrate(&ss);
    }

    print_trials_summary(trials);
    write_trials("output/ss_trials.txt", trials);
}
