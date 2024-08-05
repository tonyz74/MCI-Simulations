#include "common.h"
#include <stdio.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include <stdbool.h>

#define BINS 5000
#define BIN_LEN(st, k, b) \
    ((st->lattice[k][(b) + 1]) - (st->lattice[k][(b)]))

typedef struct ais_estimate {
    double mean;
    double var;
} ais_estimate;

typedef struct ais_state {
    integrand func;
    point pmin, pmax;
    int n_adapting_samps;
    int n_adapted_samps;
    double lattice[DIMS][BINS + 1]; 
} ais_state;

ais_state ais_init(void) {
    ais_state st;
    st.func = f;
    st.pmin = LOWER_BOUND; st.pmax = UPPER_BOUND;

    assert(N_TOTAL_SAMPS >= DIMS * BINS);
    st.n_adapting_samps = 0.3 * N_TOTAL_SAMPS;
    if (DIMS * BINS > st.n_adapting_samps) {
        st.n_adapting_samps = DIMS * BINS;
    }
    st.n_adapted_samps = N_TOTAL_SAMPS - st.n_adapting_samps;

    for (int k = 0; k < DIMS; k++) {
        double len = st.pmax.x[k] - st.pmin.x[k];
        double step_size = len / BINS;

        for (int b = 0; b < BINS; b++) {
            st.lattice[k][b] = st.pmin.x[k] + b * step_size;
        }
        st.lattice[k][BINS] = st.pmax.x[k];
    }
     
    return st;
}


ais_estimate ais_integrate_iter(ais_state *st, int N, bool adapt) {
    ais_estimate ans = { .var = 0.0, .mean = 0.0 };
    assert(N >= BINS * DIMS);

    int bcnt[DIMS][BINS] = {0};
    double bsum[DIMS][BINS] = {0};
    double sum_ord1 = 0.0, sum_ord2 = 0.0;

    for (int i = 0; i < N; i++) {
        point rand_pt = { .x = {0} };
        int rand_bins[DIMS] = {0};

        // Assuming separable, so multiply by all individual
        // PDF approximations in each dimension.
        double inv_pdf_x = 1.0;

        for (int k = 0; k < DIMS; k++) {
            int b = rand() % BINS;
            rand_bins[k] = b;
            double len = BIN_LEN(st, k, b);

            rand_pt.x[k] = randf_range(st->lattice[k][b], st->lattice[k][b + 1]);
            inv_pdf_x *= (BINS * len);
        }

        double samp = st->func(rand_pt) * inv_pdf_x;
        for (int k = 0; k < DIMS; k++) {
            bcnt[k][rand_bins[k]]++;
            bsum[k][rand_bins[k]] += samp;
        }

        sum_ord1 += samp;
        sum_ord2 += samp * samp;
    }

    ans.mean = sum_ord1 / N;
    {
        double m = ans.mean;
        double s1 = sum_ord1, s2 = sum_ord2;
        ans.var = (s2 - 2.0 * s1 * m + N * m * m) / (N - 1);
    }

    if (adapt) {
        // Resolve each dimension separately
        for (int k = 0; k < DIMS; k++) {
            double tot_mass = 0.0;

            double wt[BINS] = {0};
            double len[BINS] = {0};
            double vel[BINS] = {0};

            for (int b = 0; b < BINS; b++) {
                assert(bcnt[k][b] > 0);
                double mass_b = bsum[k][b] / bcnt[k][b];
                tot_mass += mass_b;
            }

            for (int b = 0; b < BINS; b++) {
                double mass_b = bsum[k][b] / bcnt[k][b];
                wt[b] = mass_b / tot_mass;
                len[b] = BIN_LEN(st, k, b);
                vel[b] = len[b] / wt[b];
            }

            double wt_now = 0.0;
            double pos_now = st->pmin.x[k];
            double ideal_wt = 1.0 / BINS;
            double new_pos[BINS + 1] = {0};

            int j = 0;
            for (int b = 1; b < BINS; b++) {
                while (j < BINS && wt_now + wt[j] <= ideal_wt) {
                    pos_now += len[j];
                    wt_now += wt[j];
                    wt[j] = 0;
                    len[j] = 0;
                    j++;
                }

                if (j < BINS) {
                    double wt_need = ideal_wt - wt_now;
                    double dist = wt_need * vel[j];
                    pos_now += dist;

                    wt[j] -= wt_need;
                    len[j] -= dist;

                    wt_now = 0.0;
                    new_pos[b] = pos_now;
                }
            }

            new_pos[0] = st->pmin.x[k];
            new_pos[BINS] = st->pmax.x[k];
            
            memcpy(&st->lattice[k], new_pos, sizeof(new_pos));
        }
    }

    return ans;
}

void ais_plot_pdf(ais_state *st, int k, int steps, double out[]) {
    double pos = st->pmin.x[k];
    double step = (st->pmax.x[k] - st->pmin.x[k]) / steps;

    int bin = 0;
    for (int i = 0; i < steps; i++) {
        while (bin < BINS && st->lattice[k][bin + 1] <= pos) {
            bin++;
        }

        out[i] = 1.0 / (BIN_LEN(st, k, bin) * BINS);
        pos += step;
    }
}


double ais_integrate(ais_state *st) {
    ais_estimate res_init = ais_integrate_iter(
        st, st->n_adapting_samps, true
    );
    double var_e1 = res_init.var / st->n_adapting_samps;
    
    if (st->n_adapted_samps == 0) {
        return res_init.mean;
    }

    ais_estimate res_fin = ais_integrate_iter(
        st, st->n_adapted_samps, false
    );
    double var_e2 = res_fin.var / st->n_adapted_samps;
 
    double vwt_sum = (res_init.mean / var_e1) + (res_fin.mean / var_e2);
    double tot_vwt = (1.0 / var_e1) + (1.0 / var_e2);

    return vwt_sum / tot_vwt;
}

int main(void) {
    srand(77);

    // Run multiple trials
    if (false) {
        double trials[N_TRIALS];
        for (int i = 0; i < N_TRIALS; i++) {
            int cnt = i + 1;
            if (cnt % 200 == 0) {
                printf("%.2lf%%\n", 100.0 * (double) cnt / N_TRIALS);
            }

            ais_state ig = ais_init();
            trials[i] = ais_integrate(&ig);
        }

        write_trials("output/ais_trials.txt", trials);
        print_trials_summary(trials);
    }

    // Plot a PDF (1D only)
    if (true) {
        ais_state ig = ais_init();
        ais_integrate(&ig);
        double pdf_x0[N_PDF_STEPS];
        
        FILE *out = fopen("output/ais_pdf.txt", "w");
        fprintf(out, "%d\n1\n", N_PDF_STEPS);
        fprintf(out, "%lf %lf\n", LOWER_BOUND.x[0], UPPER_BOUND.x[0]);
        
        ais_plot_pdf(&ig, 0, N_PDF_STEPS, pdf_x0);
        for (int i = 0; i < N_PDF_STEPS; i++) {
            fprintf(out, "%lf\n", pdf_x0[i]);
        }
        
        fclose(out);
    }
}
