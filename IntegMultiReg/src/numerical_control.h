#ifndef IMR_NUMERICAL_CONTROL_H
#define IMR_NUMERICAL_CONTROL_H

#include <Rinternals.h>

typedef struct {
    int historical_prior_index;
    int initial_max_iter;
    int selection_max_iter;
    int latent_max_iter;
    int prediction_max_iter;
    double tolerance;
} imr_numerical_control;

/* Parsed before native resources are allocated; immutable for one call. */
imr_numerical_control imr_read_numerical_control(SEXP control);

#endif
