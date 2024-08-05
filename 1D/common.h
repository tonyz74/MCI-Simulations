#ifndef COMMON_H
#define COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define LOWER_BOUND (-2.0)
#define UPPER_BOUND (2.0)
#define N_TOTAL_SAMPS 600000

#define N_TRIALS 100000
#define N_PDF_STEPS 100

/* * * * * * * * * * * * * * * * * * *
 * UTILITY FUNCTIONS AND STRUCTURES  *
 * * * * * * * * * * * * * * * * * * */

typedef double (*integrand)(double);

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

inline static double f(double x) {
    return exp(-x * x);
    // return (3.0 / 16.0) * (-0.5 * x * x + 2);
}

inline static const char *f_expr() {
    return "None";
    // return "lambda x: (3 / 16) * (-0.5 * x * x + 2)";
}



inline static void write_trials(const char *out_name, double trials[]) {
    FILE *out = fopen(out_name, "w");
    fprintf(out, "%d\n", N_TRIALS);
    fprintf(out, "%s\n", f_expr());
    fprintf(out, "%lf %lf\n", LOWER_BOUND, UPPER_BOUND);

    for (int i = 0; i < N_TRIALS; i++) {
        fprintf(out, "%lf\n", trials[i]);
    }

    fclose(out);
}

inline static void write_pdf(
    const char *out_name,
    double values[]
) {
    FILE *out = fopen(out_name, "w");

    fprintf(out, "%d\n", N_PDF_STEPS);
    fprintf(out, "%s\n", f_expr());
    fprintf(out, "%lf %lf\n", LOWER_BOUND, UPPER_BOUND);

    for (int i = 0; i < N_PDF_STEPS; i++) {
        fprintf(out, "%lf\n", values[i]);
    }

    fclose(out);
}

#endif
