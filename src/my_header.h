#ifndef INTEG_MULTI_REG_INTERNAL_H
#define INTEG_MULTI_REG_INTERNAL_H

#include <stdbool.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#define IMR_PI 3.1415926535897932384
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef enum {
    IMR_OUTCOME_SURVIVAL = 1,
    IMR_OUTCOME_BINARY = 2,
    IMR_OUTCOME_CONTINUOUS = 3
} imr_outcome_type;

/*
 * Internal C naming convention:
 *   - registered entry points use the `imr_` prefix;
 *   - C identifiers use snake_case;
 *   - variable-selection indicators are called gamma throughout because it
 *     matches the notation in the model.
 */

/* Sampler state updates. */
void initialize_sampler_state(
    int outcome_type, double **latent_y, double ***covariates,
    double ****features, _Bool ***gamma, int n_platforms, int *n_features,
    int n_subgroups, int **platform_models, int *n_platform_models,
    int **model_platforms, int *n_model_platforms, int *sample_size,
    double *log_likelihood, double *logdet, double *quadratic_form,
    double *slab_scale, double covariate_scale, double intercept_scale,
    double first_platform_scale, double alpha, double psi, int n_covariates);

void sample_gamma_indicators(
    int subgroup, int n_platforms, int *selected_platforms,
    int n_selected_platforms, int *n_features, int sample_size,
    double *latent_y, double **covariates, double ***features,
    _Bool ***gamma, double *log_likelihood, double *logdet,
    double *quadratic_form, double *nu, double ***theta,
    int *n_platform_models, int **platform_models, double **accept_gamma,
    gsl_rng *rng, const char *likelihood_type, double slab_scale,
    double covariate_scale, double intercept_scale, double first_platform_scale,
    int n_covariates, double alpha, double psi);

void sample_censored_latent_response(
    int subgroup, int n_platforms, int *selected_platforms,
    int n_selected_platforms, int *n_features, int sample_size,
    double *latent_y, double *observed_y, double **covariates,
    double ***features, _Bool ***gamma, double *quadratic_form,
    double *log_likelihood, int n_censored, int *censored_index,
    double logdet, gsl_rng *rng, int *n_platform_models,
    int **platform_models, double *accept_y, double slab_scale,
    double covariate_scale, double intercept_scale, double first_platform_scale,
    int n_covariates, double alpha, double psi);

void sample_binary_latent_response(
    int subgroup, int n_platforms, int *selected_platforms,
    int n_selected_platforms, int *n_features, int sample_size,
    double *latent_y, _Bool *observed_y, double **covariates,
    double ***features, _Bool ***gamma, double *log_likelihood,
    gsl_rng *rng, int *n_platform_models, int **platform_models,
    double *accept_y, double slab_scale, double covariate_scale,
    double intercept_scale, double first_platform_scale, int n_covariates,
    double alpha, double psi);

void sample_mrf_theta(
    int n_features, int n_models, double **theta, double **accept_theta,
    double *mrf_log_normalizer, _Bool **gamma, double nu, double alpha0,
    double **beta0, gsl_rng *rng);

/* Model probability, likelihood and design-matrix helpers. */
double **build_design_matrix(
    int n_covariates, int n_selected_platforms, int *n_selected_features,
    int **selected_feature_index, double **covariates, double ***features,
    int *selected_platforms, int sample_size);

double *build_posterior_precision(
    int n_coefficients, int n_covariates, int n_first_platform_features,
    int sample_size, double slab_scale, double covariate_scale,
    double intercept_scale, double first_platform_scale, double **design);

double log_likelihood_nonlocal(
    int n_coefficients, int n_covariates, int n_first_platform_features,
    int sample_size, double alpha, double psi, double *response,
    double **design, double *precision, const gsl_matrix *chol_precision,
    double *beta_mode, int moment_order, double slab_scale,
    double covariate_scale, double intercept_scale, double first_platform_scale,
    int max_iter, double tolerance, _Bool positive_beta);

void maximize_nonlocal_beta(
    double *xty, double nu, double sigma2, double *precision, int max_iter,
    double tolerance, double *beta_init, int n_coefficients, int moment_order,
    _Bool positive_beta);

double log_posterior(
    double *log_likelihood, _Bool ***gamma, double *nu, double ***theta,
    double *mrf_log_normalizer, double alpha0, double ***beta_theta,
    int n_subgroups, int n_platforms, int *n_features, int *n_platform_models);

void compute_mrf_log_normalizer(
    int n_models, double **theta, double nu, double *log_normalizer);

/* Prediction and cross-validation. */
double ***infer_posterior_models(
    double **latent_y, double ***covariates, double ****features,
    int n_draws, _Bool ****gamma_sample, double *nu, double ***theta,
    double *mrf_log_normalizer, double *slab_scale, double covariate_scale,
    double intercept_scale, double first_platform_scale, double alpha0,
    double alpha, double psi, int *n_features, int n_subgroups,
    int n_platforms, int *n_platform_models, int *n_model_platforms,
    int **model_platforms, int **platform_models, int *sample_size,
    int n_covariates, double ***beta_theta, const char *likelihood_type,
    double *posterior_weight, int *model_index, int *high_model_index,
    int *n_unique_models, int max_models);

double *predict_bma(
    int subgroup, int n_covariates, int n_selected_platforms, int sample_size,
    int *selected_platforms, int *n_platform_models, int **platform_models,
    int *n_features, double **covariates, double ***features,
    _Bool ****gamma_sample, double ***beta, double *posterior_weight,
    int max_models, int *model_index, int *high_model_index, int n_draws,
    int outcome_type);

/* Small numerical utilities shared by sampler and prediction code. */
void find_indices_not_equal(int n, _Bool *values, int excluded_value,
                            int *index, int *n_found);

void propose_gamma_update(int n_features, _Bool *current_gamma,
                          _Bool *proposal_gamma, float flip_probability,
                          gsl_rng *rng);

double cholesky_logdet(gsl_matrix *chol);

double cholesky_quadratic_form(
    int n, double matrix[n * n], double vector[n], double *logdet);

void sort_descending_index(int n, double *values, int *index);

#endif
