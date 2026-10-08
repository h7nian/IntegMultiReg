#ifndef IMR_MARGINAL_SAMPLER_H
#define IMR_MARGINAL_SAMPLER_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* R's vector arithmetic rounds between operations. Do not fuse those steps. */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

/* Same capacity as .imr_check_mrf_capacity() and the joint sampler. */
#define IMR_MARGINAL_MAX_SUBGROUPS 16

enum {
    IMR_SWAP_PROPOSALS, IMR_SWAP_ACCEPTS,
    IMR_INTERACTION_PROPOSALS, IMR_INTERACTION_ACCEPTS,
    IMR_VARIANCE_PROPOSALS, IMR_VARIANCE_ACCEPTS,
    IMR_LATENT_PROPOSALS, IMR_LATENT_ACCEPTS
};

typedef struct {
    int n, p, forced, outcome;
    SEXP design;
    double *beta, *response, *observed, *tau, *mean, *other, *scratch;
    int *selected, *status;
    double shape, rate, variance, variance_step;
    double *coefficient_draws, *variance_draws, *latent_draws;
} imr_marginal_group;

typedef struct {
    int members, features, patterns_count;
    int *groups, *columns, *permutation, *pair_pool;
    SEXP patterns, theta, selected_matrix;
    double nu, *energy, *weights, *theta_draws;
} imr_marginal_platform;

typedef struct {
    int groups_count, platforms_count, draws, burnin, thin, sharing, verbose;
    int coefficients_marginalized;
    double theta_step, swap_rate, interaction_shape, interaction_rate;
    imr_marginal_group *groups;
    imr_marginal_platform *platforms;
    SEXP roots, result, multiply, crossprod, largest_eigenvalue, prepare_model;
    double *log_density, *acceptance;
} imr_marginal_state;

SEXP imr_marginal_element(SEXP list, const char *name);
SEXP imr_marginal_product(imr_marginal_state *state, SEXP x, SEXP y, int cross);
double imr_marginal_sum(const double *values, int n);
double imr_marginal_log_sum_exp(double *values, int n);
void imr_marginal_energy(imr_marginal_state *state, imr_marginal_platform *platform, double *out);
int imr_marginal_select(imr_marginal_state *state, imr_marginal_platform *platform, double *log_bf);
void imr_marginal_pair(imr_marginal_platform *platform, int *pair);
int imr_marginal_accept(double log_ratio);
void imr_marginal_interaction(imr_marginal_state *state, imr_marginal_platform *platform);
void imr_marginal_store(imr_marginal_state *state, int row, double score);
void imr_marginal_progress(imr_marginal_state *state, int iteration);
void imr_variance_run(imr_marginal_state *state);
void imr_coefficients_run(imr_marginal_state *state);

#endif
