#include "common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#define BINS 6
#define BIN_LEN(st, b) ((st)->lattice[(b) + 1] - (st)->lattice[b])

typedef struct ais1d_estimate {
    double mean;
    double var;
} ais1d_estimate;

typedef struct ais1d_state {
    integrand func;
    double low, high;

    int n_adapting_samps;
    int n_adapted_samps;
    double lattice[BINS + 1];
} ais1d_state;

ais1d_state ais1d_init(
    integrand fn,
    double L, double H,
    int n_samps
) {
    ais1d_state st;
    st.low = L;
    st.high = H;
    st.func = fn;
    st.n_adapting_samps = 1.0 * n_samps;
    st.n_adapted_samps = n_samps - st.n_adapting_samps;

    double tot_domain = (H - L);
    double each_size = tot_domain / BINS;
    for (int i = 0; i <= BINS; i++) {
        st.lattice[i] = (double) i * each_size + st.low;
    }

    return st;
}

void ais1d_sample_in_bin(
    ais1d_state *st, int b,
    int bcnt[], double bsum[],
    double *sum_ord1, double *sum_ord2
) {
    double x = randf_range(st->lattice[b], st->lattice[b + 1]);
    double inv_pdf_x = (BINS * BIN_LEN(st, b));

    double samp = st->func(x) * inv_pdf_x;
    bcnt[b]++; bsum[b] += samp;

    *sum_ord1 += samp;
    *sum_ord2 += samp * samp;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Calculate an integral estimate and the variance using   *
 * importance sampling with the current approximation of   *
 * the PDF. Optionally refine the PDF approximation using  *
 * these new samples.                                      *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
ais1d_estimate ais1d_integrate_iter(
    ais1d_state *st,
    int N, bool adapt
) {
    ais1d_estimate res = (ais1d_estimate) {
        .mean = 0.0, .var = 0.0
    };

    double sum_ord1 = 0.0;
    double sum_ord2 = 0.0;

    int bcnt[BINS] = {0};
    double bsum[BINS] = {0};

    // If it's going to adapt, each bin needs at least one 
    // sample to ensure there is no division by 0, since
    // we must have bcnt[b] != 0 for all b.
    if (adapt) {
        assert(N >= BINS);

        for (int b = 0; b < BINS; b++) {
            ais1d_sample_in_bin(
                st, b,
                bcnt, bsum,
                &sum_ord1, &sum_ord2
            );
        }

        N -= BINS;
    }

    for (int i = 0; i < N; i++) {
        int b = rand() % BINS;
        ais1d_sample_in_bin(
            st, b,
            bcnt, bsum,
            &sum_ord1, &sum_ord2
        );
    }

    res.mean = sum_ord1 / N;
    {
        double m = res.mean;
        double s1 = sum_ord1;
        double s2 = sum_ord2;
        res.var = (s2 - 2.0 * m * s1 + m * m * N) / (N - 1);
    }

    if (adapt) {
        double wt[BINS] = {0};
        double len[BINS] = {0};
        double vel[BINS] = {0};

        // Computes the weights and additional information for each 
        // interval. We assume weight is equally distributed within the bin
        // so `vel` is the velocity, where distance is the length of the
        // interval and time is the weight.

        double tot_mass = 0;
        for (int b = 0; b < BINS; b++) {
            double b_mass = bsum[b] / bcnt[b];
            tot_mass += b_mass;
        }
        printf("tot mass %lf\n", tot_mass);

        for (int b = 0; b < BINS; b++) {
            double b_mass = bsum[b] / bcnt[b];
            wt[b] = b_mass / tot_mass;
            printf("in bin b %d: %lf, mass %lf\n", b, wt[b], b_mass);
            len[b] = BIN_LEN(st, b);
            vel[b] = len[b] / wt[b];
        }

        // Refines the different bin positions to ensure that weights
        // are close to identical in each interval.

        double wt_now = 0.0;
        double pos_now = st->low;
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

        new_pos[0] = st->low;
        new_pos[BINS] = st->high;
        
        memcpy(st->lattice, new_pos, sizeof(new_pos));
    }
    
    return res;
}

double ais1d_integrate(ais1d_state *st) {

    // Take `n_adapting_samps` uniform samples. We start off
    // our approximation as a uniform distribution.
    ais1d_estimate res_init = ais1d_integrate_iter(
        st, st->n_adapting_samps, true
    );

    // Now we have learned the function, we try integrating
    // again, this time doing importance sampling using our
    // improved approximation of the PDF.
    
    // If there are no samples allocated to this, skip.
    if (st->n_adapted_samps == 0) {
        return res_init.mean;
    }

    ais1d_estimate res_fin = ais1d_integrate_iter(
        st, st->n_adapted_samps, false
    );

    // Our final estimate for the integral combines these two
    // results, weighted by the inverse variances.

    double vwt_sum = (res_init.mean / res_init.var)
        + (res_fin.mean / res_fin.var);
    double tot_vwt = (1.0 / res_init.var)
        + (1.0 / res_fin.var);

    return vwt_sum / tot_vwt;
}

void ais1d_plot_pdf(ais1d_state *st, int steps, double out[]) {
    double pos = st->low;
    double step = (st->high - st->low) / steps;

    int bin = 0;
    for (int i = 0; i < steps; i++) {
        while (bin < BINS && st->lattice[bin + 1] <= pos) {
            bin++;
        }

        out[i] = 1.0 / (BIN_LEN(st, bin) * BINS);
        pos += step;
    }
}



enum output_mode {
    PDF,
    TRIALS,
    ONE,
    SAVE
};

int main(void) {
    srand(time(NULL));

    enum output_mode mode = PDF;

    if (mode == ONE) {
        ais1d_state ig = ais1d_init(
            f, LOWER_BOUND, UPPER_BOUND, N_TOTAL_SAMPS
        );
        ais1d_integrate(&ig);
    }

    if (mode == TRIALS) {
        double trials[N_TRIALS];
        for (int i = 0; i < N_TRIALS; i++) {
            ais1d_state ig = ais1d_init(
                f, LOWER_BOUND, UPPER_BOUND, N_TOTAL_SAMPS
            );
            trials[i] = ais1d_integrate(&ig);
        }

        write_trials("output/ais1d_trials.txt", trials);
    }

    if (mode == PDF) {
        ais1d_state ig = ais1d_init(
            f, LOWER_BOUND, UPPER_BOUND, N_TOTAL_SAMPS
        );
        ais1d_integrate(&ig);

        double pdf[N_PDF_STEPS];
        ais1d_plot_pdf(&ig, N_PDF_STEPS, pdf);

        write_pdf("output/ais1d_pdf.txt", pdf);
    }

    if (mode == SAVE) {
        FILE *sav = fopen("output/ais1d_save.txt", "w");
        ais1d_state ig = ais1d_init(
            f, LOWER_BOUND, UPPER_BOUND, N_TOTAL_SAMPS
        );
        ais1d_integrate(&ig);
        
        for (int i = 0; i <= BINS; i++) {
            fprintf(sav, "%lf\n", ig.lattice[i]);
        }

        fclose(sav);
    }
}
