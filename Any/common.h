#ifndef COMMON_H
#define COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define DIMS 1

#define LOWER_BOUND ((point) { { -2 } })
#define UPPER_BOUND ((point) { { 2 } })
#define N_TOTAL_SAMPS 1000000

#define N_TRIALS 5000
#define N_PDF_STEPS 5000

/* * * * * * * * * * * * * * * * * * *
 * UTILITY FUNCTIONS AND STRUCTURES  *
 * * * * * * * * * * * * * * * * * * */

typedef struct point {
    double x[DIMS];
} point;

typedef double (*integrand)(point);

inline static double randf() {
    return ((double) rand() / (double) RAND_MAX);
}

inline static double randf_range(double a, double b) {
    assert(b > a);
    return randf() * (b - a) + a;
}


/* * * * * * * * * * * *
 * CUSTOM DEFINITIONS  *
 * * * * * * * * * * * */

inline static double f(point p) {
    // return 1.0 / ((1 + p.x[0] * p.x[0]) * (1 + p.x[1] * p.x[1]));
    return 1.0 / sqrt(2 * M_PI * 1) * exp(-0.5 * (p.x[0] * p.x[0]));
}

inline static const char *f_expr() {
    return "None";
}



/* * * * * * * * * * * * * * * * * * *
 * STRUCTURED DATA OUTPUT            *
 * * * * * * * * * * * * * * * * * * */

inline static void _write_point(FILE *out, point pt) {
    fprintf(out, "[");
    for (int k = 0; k < DIMS; k++) {
        fprintf(out, "%lf", pt.x[k]);
        if (k != DIMS - 1) fprintf(out, ",");
    }
    fprintf(out, "]");
}

inline static void write_trials(const char *out_name, double trials[]) {
    FILE *out = fopen(out_name, "w");
    fprintf(out, "%d\n", N_TRIALS);
    fprintf(out, "%s\n", f_expr());
    fprintf(out, "0.0 0.0\n");

    for (int i = 0; i < N_TRIALS; i++) {
        fprintf(out, "%lf\n", trials[i]);
    }

    fclose(out);
}

inline static void print_trials_summary(double trials[]) {
    double mean = 0.0;
    for (int i = 0; i < N_TRIALS; i++) {
        mean += trials[i] / (N_TRIALS);
    }

    double svar = 0.0;
    for (int i = 0; i < N_TRIALS; i++) {
        svar += (trials[i] - mean) * (trials[i] - mean) / (N_TRIALS - 1);
    }

    printf("[MEAN] = %lf\n", mean);
    printf("[SVAR] = %lf\n", svar);
}

#endif
