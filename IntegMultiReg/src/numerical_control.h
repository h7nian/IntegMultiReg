#ifndef IMR_NUMERICAL_CONTROL_H
#define IMR_NUMERICAL_CONTROL_H

#include <Rinternals.h>

enum imr_laplace_stage {
    IMR_LAPLACE_INITIAL, IMR_LAPLACE_SELECTION, IMR_LAPLACE_LATENT,
    IMR_LAPLACE_PREDICTION
};
enum imr_laplace_measure {
    IMR_LAPLACE_CALLS, IMR_LAPLACE_LIMIT, IMR_LAPLACE_NONFINITE,
    IMR_LAPLACE_FACTORIZATION_FAILURE
};
typedef struct {
    int n_subgroups;
    int subgroup;
    double *values;
} imr_laplace_diagnostics;

typedef struct {
    int historical_prior_index;
    int initial_max_iter;
    int selection_max_iter;
    int latent_max_iter;
    int prediction_max_iter;
    double tolerance;
    imr_laplace_diagnostics *diagnostics;
} imr_numerical_control;

/* Numerical choices are immutable; an optional fit-owned counter is mutable. */
imr_numerical_control imr_read_numerical_control(SEXP control);

/* Fit-owned counters only; prediction callers leave diagnostics NULL. */
void imr_record_laplace(const imr_numerical_control *control, int stage, int measure);

#endif
