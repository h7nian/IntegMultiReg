/* Historical post-fit CV, adapted from package 0.1.0 (c523483).
 * legacy: ranked distinct full-fit selection models and historical scoring.
 * importance: all retained proposal draws, including their multiplicities.
 * Both deliberately condition on full-fit preprocessing and latent summaries.
 * Keep this compatibility implementation isolated from training-fold refits.
 */
#include <Rinternals.h>
#include <R.h>
#include <time.h>
#include <stdbool.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include "my_header.h"
#include "utils.h"

static double legacy_max(int n, double *values);

/* Keys refer to immutable draws owned by this call. Hash equality is never
 * sufficient: verify every selection indicator before reusing a result. */
typedef struct {
    uint64_t hash;
    int row;
} postfit_cache_entry;

static uint64_t postfit_state_hash(_Bool ***state, int n_platforms,
                                  const int *n_platform_models,
                                  const int *n_features)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (int p = 0; p < n_platforms; ++p)
        for (int g = 0; g < n_platform_models[p]; ++g)
            for (int j = 0; j < n_features[p]; ++j) {
                hash ^= (uint64_t)state[p][g][j];
                hash *= UINT64_C(1099511628211);
            }
    return hash;
}

static int postfit_same_state(_Bool ***left, _Bool ***right, int n_platforms,
                              const int *n_platform_models,
                              const int *n_features)
{
    for (int p = 0; p < n_platforms; ++p)
        for (int g = 0; g < n_platform_models[p]; ++g)
            if (memcmp(left[p][g], right[p][g],
                       (size_t)n_features[p] * sizeof(_Bool)) != 0) return 0;
    return 1;
}

/*
 * Identify repeated immutable selection states once per native call. The
 * temporary hash table and retained index together obey the cache budget.
 */
static int *postfit_model_representatives(_Bool ****gamma_sample, int n_draws,
                                         int n_platforms, const int *n_platform_models,
                                         const int *n_features, size_t cache_bytes,
                                         int cache_hash_bits)
{
    size_t index_bytes = (size_t)n_draws * sizeof(int);
    if (index_bytes > cache_bytes ||
        cache_bytes - index_bytes < 2 * sizeof(postfit_cache_entry)) return NULL;
    int *representatives = (int *)R_alloc(n_draws, sizeof(int));
    void *allocation_mark = vmaxget();
    size_t limit = (cache_bytes - index_bytes) / sizeof(postfit_cache_entry);
    size_t slots = 2, used = 0;
    while (slots < (size_t)n_draws * 2 && slots <= limit / 2) slots *= 2;
    postfit_cache_entry *cache = (postfit_cache_entry *)R_alloc(slots, sizeof(*cache));
    for (size_t slot = 0; slot < slots; ++slot) cache[slot].row = -1;
    for (int draw = 0; draw < n_draws; ++draw) {
        _Bool ***state = gamma_sample[n_draws - 1 - draw];
        uint64_t hash = postfit_state_hash(state, n_platforms, n_platform_models, n_features);
        if (cache_hash_bits == 0) hash = 0;
        else if (cache_hash_bits < 64) hash &= UINT64_MAX >> (64 - cache_hash_bits);
        size_t slot = (size_t)hash & (slots - 1);
        while (cache[slot].row >= 0) {
            if (cache[slot].hash == hash &&
                postfit_same_state(state, gamma_sample[n_draws - 1 - cache[slot].row],
                                   n_platforms, n_platform_models, n_features)) break;
            slot = (slot + 1) & (slots - 1);
        }
        representatives[draw] = cache[slot].row >= 0 ? cache[slot].row : draw;
        /* Bound probing at half capacity; uncached states are still computed. */
        if (cache[slot].row < 0 && used < slots / 2) {
            cache[slot].hash = hash;
            cache[slot].row = draw;
            ++used;
        }
    }
    vmaxset(allocation_mark);
    return representatives;
}

/* Fit fold coefficients and retain the original draw-wise averaging order. */
static double *postfit_predict_fold(int outcome_type, int subgroup, int n_covariates, int n_selected_platforms, int *selected_platforms,
                int *n_platform_models, int **platform_models, int *n_features,
                int model_sample_size, int test_sample_size, int *test_index, int *train_index,
                double *latent_response, double **covariates, double ***features, _Bool ****gamma_sample, double residual_shape, double residual_rate,
                int max_models, int *model_index, int *high_model_index, int n_draws,
                int importance, const int *model_representatives)
{

  int l, i, j, j1, in, i1, i2;
  double *weight = malloc(max_models * sizeof(double));
  // Latent predictions, transformed to probabilities for binary outcomes.
  double *prediction = dvector(0, test_sample_size - 1);
  for (i = 0; i < test_sample_size; i++)
    prediction[i] = 0;
  double **model_predictions = dmatrix(0, max_models - 1, 0, test_sample_size - 1);
  for (l = 0; l < max_models; l++)
  {
    int ranked_model = high_model_index[l];
    int draw_index = model_index[ranked_model];
    if (model_representatives != NULL) {
      int row = model_representatives[l];
      if (row < l && R_FINITE(weight[row])) {
        weight[l] = weight[row];
        memcpy(model_predictions[l], model_predictions[row],
               (size_t)test_sample_size * sizeof(double));
        continue;
      }
    }

    int **selected_feature_index = malloc(n_selected_platforms * sizeof(int *));
    int n_selected_features[n_selected_platforms];
    for (int i = 0; i < n_selected_platforms; i++)
    {
      int platform_index = selected_platforms[i];
      int platform_model_index = -1;
      for (int ss = 0; ss < n_platform_models[platform_index]; ss++)
      {
        if (platform_models[platform_index][ss] == subgroup)
        {
          platform_model_index = ss;
          break;
        }
      }
      if (platform_model_index == -1)
      {
        Rf_error("Subgroup %d not found for platform %d\n", 1 + subgroup, 1 + platform_index);
      }
      selected_feature_index[i] = malloc(n_features[platform_index] * sizeof(int));
      if (!selected_feature_index[i])
      {
        Rf_error("malloc failed for selected_feature_index[%d]\n", i);
      }
      n_selected_features[i] = 0;
      find_indices_not_equal(n_features[platform_index], gamma_sample[n_draws - 1 - draw_index][platform_index][platform_model_index], 0, selected_feature_index[i], &n_selected_features[i]);
    }

    double **design = build_design_matrix(n_covariates, n_selected_platforms, n_selected_features, selected_feature_index, covariates, features, selected_platforms, model_sample_size);

    int train_sample_size = model_sample_size - test_sample_size;
    int n_selected_total = 0;
    for (int p = 0; p < n_selected_platforms; p++)
    {
      n_selected_total += n_selected_features[p];
    }
    int n_coefficients = 1 + n_covariates + n_selected_total;

    double *precision = malloc(n_coefficients * n_coefficients * sizeof(double));
    for (j = 0; j < n_coefficients; j++)
    {
      for (j1 = 0; j1 <= j; j1++)
      {
        double a = 0;
        if ((j == 0) && (j1 == 0))
        {
          a = train_sample_size;
        }
        else if (j1 == 0)
        {
          for (i2 = 0; i2 < train_sample_size; i2++)
          {
            i1 = train_index[i2];
            a += design[i1][j - 1];
          }
        }
        else if ((j != 0) && (j1 != 0))
        {
          for (i2 = 0; i2 < train_sample_size; i2++)
          {
            i1 = train_index[i2];
            a += design[i1][j - 1] * design[i1][j1 - 1];
          }
        }
        if (j == j1)
          a += .001; // To always make the matrix positive definite
        precision[j * n_coefficients + j1] = precision[j1 * n_coefficients + j] = a;
      }
    }

    double *xty = calloc(n_coefficients, sizeof(double));
    for (j = 0; j < n_coefficients; j++)
    {
      double a1 = 0;
      for (i2 = 0; i2 < train_sample_size; i2++)
      {
        i1 = train_index[i2];
        if (j == 0)
          a1 += latent_response[i1];
        else
          a1 += design[i1][j - 1] * latent_response[i1];
      }
      xty[j] = a1;
    }
    gsl_vector_view b = gsl_vector_view_array(xty, n_coefficients);

    gsl_vector *x = gsl_vector_alloc(n_coefficients);
    gsl_matrix_view Aip = gsl_matrix_view_array(precision, n_coefficients, n_coefficients);
    int status = gsl_linalg_cholesky_decomp(&Aip.matrix);
    if (status)
    {
      Rprintf("Cholesky failed (subgroup %d): %s\n", subgroup, gsl_strerror(status));
      weight[l] = -INFINITY;
      gsl_vector_free(x);
      free(xty);
      free(precision);
      for (i = 0; i < model_sample_size; i++)
        free(design[i]);
      free(design);
      for (int platform = 0; platform < n_selected_platforms; platform++)
        free(selected_feature_index[platform]);
      free(selected_feature_index);
      continue;
    }

    gsl_linalg_cholesky_solve(&Aip.matrix, &b.vector, x);

    double *beta = malloc(n_coefficients * sizeof(double));
    for (j = 0; j < n_coefficients; j++)
    {
      beta[j] = gsl_vector_get(x, j);
    }
    gsl_vector_free(x);
    double test_sse = 0;
    for (in = 0; in < test_sample_size; in++)
    {
      model_predictions[l][in] = 0;
      i = test_index[in];
      for (j = 0; j < n_coefficients; j++)
      {
        if (j == 0)
          model_predictions[l][in] += beta[0];
        else
          model_predictions[l][in] += design[i][j - 1] * beta[j];
      }
      test_sse += pow(latent_response[i] - model_predictions[l][in], 2);
    }
    double train_sse = 0;
    double *train_predictions = malloc(train_sample_size * sizeof(double));
    for (in = 0; in < train_sample_size; in++)
    {
      train_predictions[in] = 0;
      i = train_index[in];
      for (j = 0; j < n_coefficients; j++)
      {
        if (j == 0)
          train_predictions[in] += beta[0];
        else
          train_predictions[in] += design[i][j - 1] * beta[j];
      }
      train_sse += pow(latent_response[i] - train_predictions[in], 2);
    }
    free(beta);
    free(train_predictions);
    free(xty);
    /* The historical code truncates this value to an integer. The paper
     * density uses the fractional degrees of freedom implied by the prior. */
    double degrees_freedom = 2 * residual_shape + train_sample_size;
    if (!importance) degrees_freedom = (int)degrees_freedom;
    double residual_scale = (residual_rate + train_sse) / degrees_freedom;
    weight[l] = (test_sample_size / 2.0) * log(residual_scale) + 0.5 * (2 * residual_shape + model_sample_size) * log(1 + test_sse / (degrees_freedom * residual_scale));
    for (i = 0; i < model_sample_size; i++)
      free(design[i]);
    free(design);
    free(precision);
    for (int i = 0; i < n_selected_platforms; i++)
    {
      free(selected_feature_index[i]);
    }
    free(selected_feature_index);
  }
  double max_log_weight = legacy_max(max_models, weight);
  double weight_sum = 0;
  for (l = 0; l < max_models; l++)
  {
    weight[l] = exp(weight[l] - max_log_weight);
    weight_sum += weight[l];
  }
  double *probability = NULL;
  if (outcome_type == IMR_OUTCOME_BINARY) // binary outcome
  {
    // Normalize the model-averaging weights once (not inside the per-test-point
    // loop, which would repeatedly divide the shared weights by weight_sum).
    for (l = 0; l < max_models; l++)
    {
      weight[l] = weight[l] / weight_sum;
    }
    probability = malloc(test_sample_size * sizeof(double));
    for (i = 0; i < test_sample_size; i++)
    {
      probability[i] = 0.0; // must be initialized before accumulating below
      double *log_prob = malloc(max_models * sizeof(double));

      for (l = 0; l < max_models; l++)
      {
        log_prob[l] = -log(2) + gsl_sf_log_erfc(-model_predictions[l][i] / sqrt(2));
      }
      double max_log_probability = legacy_max(max_models, log_prob);
      for (l = 0; l < max_models; l++)
      {
        probability[i] += weight[l] * exp(log_prob[l] - max_log_probability);
      }
      probability[i] = probability[i] * exp(max_log_probability);
      free(log_prob);
    }
  }
  else if ((outcome_type == IMR_OUTCOME_SURVIVAL)||(outcome_type == IMR_OUTCOME_CONTINUOUS)) // survival or continuous outcome
  {
    for (in = 0; in < test_sample_size; in++)
    {
      double a = 0;
      for (l = 0; l < max_models; l++)
      {
        a += model_predictions[l][in] * weight[l];
      }
      prediction[in] = a / weight_sum;
    }
  }

  free(weight);

  free_dmatrix(model_predictions, 0, max_models - 1, 0, test_sample_size - 1);

  if ((outcome_type == IMR_OUTCOME_SURVIVAL) || (outcome_type == IMR_OUTCOME_CONTINUOUS)) // survival outcome or continuous
  {
    return prediction;
  }
  else
  {
    free(prediction);
    return probability;
  }
}

/* Historical 0.1.0 concordance, including its original tie convention. */
static double legacy_concordance(int n, double *prediction, double *observed_time, _Bool *event)
{
  int i, j;
  double concordance_denominator = 0;
  double concordance_numerator = 0;
  double time1, time2, prediction1, prediction2;
  for (i = 0; i < n; i++)
  {
    time1 = observed_time[i];
    prediction1 = prediction[i];
    for (j = 0; j < n; j++)
    {
      if (i != j)
      {
        time2 = observed_time[j];
        prediction2 = prediction[j];
        concordance_numerator +=
            (prediction2 > prediction1) * (time2 > time1) * (event[i] == 1) +
            (prediction2 < prediction1) * (time2 < time1) * (event[j] == 1) +
            0.5 * ((prediction2 == prediction1) || (time2 == time1)) * (event[i] == 1) * (event[j] == 0) +
            0.5 * ((prediction2 == prediction1) || (time2 == time1)) * (event[j] == 1) * (event[i] == 0);
        concordance_denominator +=
            (time2 > time1) * (event[i] == 1) +
            (time2 < time1) * (event[j] == 1) +
            (time2 == time1) * (event[i] == 1) * (event[j] == 0) +
            (time2 == time1) * (event[i] == 0) * (event[j] == 1);
      }
    }
  }
  return concordance_denominator > 0 ?
      concordance_numerator / concordance_denominator : NA_REAL;
}


static void postfit_partition(int fold, int n_folds, int model_sample_size, int *test_sample_size, int *censored_index, int n_censored, int *uncensored_index, int *test_index, int *train_index)
{
  int n_uncensored = model_sample_size - n_censored;
  int i = 0;
  int i1 = 0;
  int i2 = 0;

  for (i = 0; i < n_censored; i++)
  {
    if ((fold * n_censored / n_folds > i) || ((fold + 1) * n_censored / n_folds <= i))
    {
      train_index[i1] = censored_index[i];
      i1++;
    }
    else
    {
      test_index[i2] = censored_index[i];
      i2++;
    }
  }
  for (i = 0; i < n_uncensored; i++)
  {
    if ((fold * n_uncensored / n_folds > i) || ((fold + 1) * n_uncensored / n_folds <= i))
    {
      train_index[i1] = uncensored_index[i];
      i1++;
    }
    else
    {
      test_index[i2] = uncensored_index[i];
      i2++;
    }
  }
  *test_sample_size = i2;
}

static double legacy_auc(int n, double *esti, _Bool * class)
{
    if (n == 0) return NA_REAL;
    int positives = 0;
    for (int i = 0; i < n; ++i) positives += class[i];
    if (positives == 0 || positives == n) return NA_REAL;
    double fpr[n + 2], tpr[n + 2];
    double auc1 = 0;
    int i, j;
    double esti1[n];
    for (i = 0; i < n; i++)
    {
        esti1[i] = esti[i];
    }
    int idx[n];
    sort_descending_index(n, esti1, idx);

    fpr[n + 1] = 1;
    tpr[n + 1] = 1;
    fpr[0] = 0;
    tpr[0] = 0;
    for (i = n; i >= 1; --i)
    {
        double af = 0;
        double at = 0;
        for (j = 0; j < n; j++)
        {
            if (esti[j] > esti1[i - 1])
            {
                if (class[j] == 0)
                {
                    af += 1;
                }
                else
                {
                    at += 1;
                }
            }
        }
        tpr[i] = at / positives;
        fpr[i] = af / (n - positives);
        auc1 += (fpr[i + 1] - fpr[i]) * (tpr[i + 1] + tpr[i]);
    }
    auc1 += (fpr[1] - fpr[0]) * (tpr[1] + tpr[0]);
    auc1 = 0.5 * (auc1);
    return auc1;
}


static double legacy_max(int n, double *values) {
    double result = -INFINITY;
    for (int i = 0; i < n; ++i) if (values[i] > result) result = values[i];
    return result;
}
static double legacy_mse(int n, double *prediction, double *response) {
    if (n == 0) return NA_REAL;
    double total = 0;
    for (int i = 0; i < n; ++i)
        total += (prediction[i] - response[i]) * (prediction[i] - response[i]);
    return total / n;
}
/*
 * Cross-validation prediction entry point called from R.
 *
 * Consume full-fit draws and preprocessing. Legacy ranks distinct models;
 * importance retains every proposal draw. Both use the historical GSL
 * partitions. R reconstructs subject records and scores importance predictions
 * with the current metric definitions.
 */
SEXP imr_cv_postfit(SEXP forced_scale_R, SEXP molecular_scale_R,
                            SEXP residual_shape_R, SEXP residual_rate_R,
                            SEXP interaction_shape_R, SEXP interaction_rate_R,
                            SEXP seed_R, SEXP nu_R,
                            SEXP latent_y_R, SEXP gamma_sample_R, SEXP theta_R,
                            SEXP n_platforms_R,
                            SEXP platform_models_R, SEXP model_platforms_R, SEXP n_subgroups_R,
                            SEXP sample_sizes_R, SEXP n_features_R, SEXP n_covariates_R,
                            SEXP features_R, SEXP response_R, SEXP covariates_R,
                            SEXP outcome_type_R,
                            SEXP draws_R, SEXP folds_R_input, SEXP rounds_R,
                            SEXP max_models_R, SEXP importance_R, SEXP cache_control_R)
{
    if (!isReal(cache_control_R) || XLENGTH(cache_control_R) != 2 ||
        !R_FINITE(REAL(cache_control_R)[0]) || !R_FINITE(REAL(cache_control_R)[1]) ||
        REAL(cache_control_R)[0] < 0 || REAL(cache_control_R)[0] > 128.0 * 1024 * 1024 ||
        REAL(cache_control_R)[1] < 0 || REAL(cache_control_R)[1] > 64 ||
        floor(REAL(cache_control_R)[0]) != REAL(cache_control_R)[0] ||
        floor(REAL(cache_control_R)[1]) != REAL(cache_control_R)[1])
        Rf_error("Invalid internal post-fit cache controls");
    size_t cache_bytes = (size_t)REAL(cache_control_R)[0];
    int cache_hash_bits = (int)REAL(cache_control_R)[1];
    clock_t started = clock();
    int importance = asLogical(importance_R);
    int protect_count = 0;

    double forced_scale = REAL(forced_scale_R)[0];
    double molecular_scale = REAL(molecular_scale_R)[0];
    double residual_shape = REAL(residual_shape_R)[0];
    double residual_rate = REAL(residual_rate_R)[0];
    double interaction_shape = REAL(interaction_shape_R)[0];
    double interaction_rate = REAL(interaction_rate_R)[0];

    PROTECT(platform_models_R);
    protect_count++;
    PROTECT(model_platforms_R);
    protect_count++;
    int n_cv_rounds = asInteger(rounds_R);
    int n_draws = asInteger(draws_R);
    int n_subgroups = asInteger(n_subgroups_R);
    int max_models_requested = asInteger(max_models_R); // Maximum number of models to be used for prediction and Bayesian model averaging

    int n_unique_models;
    int outcome_type = asInteger(outcome_type_R);
    double *nu = REAL(nu_R);
    int n_platforms = asInteger(n_platforms_R);
    int n_folds = asInteger(folds_R_input);
    int fold;
    double *post = dvector(0, n_draws - 1);
    int *model_index = malloc(n_draws * sizeof(int));
    int *high_model_index = malloc(n_draws * sizeof(int));
    double ***response = r_list_matrix_to_c(n_subgroups, response_R);
    double **latent_response = r_list_vector_double_to_c(n_subgroups, latent_y_R);
    double ***covariates = r_list_matrix_to_c(n_subgroups, covariates_R);
    double ****features = r_list_list_matrix_to_c(n_subgroups, features_R);
    _Bool ****gamma_sample = r_list_list_matrix_to_c_bool(n_draws, gamma_sample_R);

    double ***theta = r_list_matrix_to_c(n_platforms, theta_R);
    int **platform_models_c = malloc(n_platforms * sizeof(int *));
    int *n_platform_models_c = malloc(n_platforms * sizeof(int));

    for (int i = 0; i < n_platforms; i++)
    {
        SEXP mapping_R = VECTOR_ELT(platform_models_R, i);
        int mapping_size = LENGTH(mapping_R);
        n_platform_models_c[i] = mapping_size;
        platform_models_c[i] = INTEGER(mapping_R);
    }

    int **model_platforms_c = malloc(n_subgroups * sizeof(int *));
    int *n_model_platforms_c = malloc(n_subgroups * sizeof(int));

    for (int i = 0; i < n_subgroups; i++)
    {
        SEXP mapping_R = VECTOR_ELT(model_platforms_R, i);
        int mapping_size = LENGTH(mapping_R);
        n_model_platforms_c[i] = mapping_size;
        model_platforms_c[i] = INTEGER(mapping_R);
    }
    double *mrf = calloc(n_platforms, sizeof(double));
    for (int l = 0; l < n_platforms; l++)
    {
        compute_mrf_log_normalizer(n_platform_models_c[l], theta[l], nu[l], &mrf[l]);
    }
    double *slab_scales = malloc(n_subgroups * sizeof(double));
    for (int i = 0; i < n_subgroups; i++)
    {
        slab_scales[i] = molecular_scale;
    }
    double covariate_scale = forced_scale;
    double first_platform_scale = molecular_scale;
    int *n_features = INTEGER(n_features_R);
    int *sample_size_ptr = INTEGER(sample_sizes_R);
    int n_covariates = asInteger(n_covariates_R);
    double ***interaction_rates = calloc(n_platforms, sizeof(double **));
    for (int l = 0; l < n_platforms; l++)
    {
        interaction_rates[l] = calloc(n_platform_models_c[l], sizeof(double *));
        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            interaction_rates[l][m] = calloc(n_platform_models_c[l], sizeof(double));
            for (int m1 = 0; m1 < n_platform_models_c[l]; m1++)
            {
                interaction_rates[l][m][m1] = interaction_rate;
            }
        }
    }

    const char likelihood_type[] = "NonLocal";
    double ***beta = NULL;
    if (!importance) beta = infer_posterior_models(latent_response, covariates, features, n_draws, gamma_sample,
                                     nu, theta, mrf, slab_scales, covariate_scale, forced_scale,
                                     first_platform_scale, interaction_shape,
                                     residual_shape, residual_rate, n_features, n_subgroups, n_platforms,
                                     n_platform_models_c, n_model_platforms_c, model_platforms_c,
                                     platform_models_c, sample_size_ptr, n_covariates,
                                     interaction_rates, likelihood_type, post, model_index,
                                     high_model_index, &n_unique_models, max_models_requested);

    int max_models;
    if (importance) {
        /* Retain empirical multiplicities: each retained MCMC state is one
         * proposal draw in the importance average (paper equations 6--7). */
        n_unique_models = 0;
        max_models = n_draws;
        for (int s = 0; s < n_draws; ++s)
            model_index[s] = high_model_index[s] = s;
    } else {
        max_models = MIN(n_unique_models, max_models_requested);
    }
    /* Share only immutable state identities, never fold-specific estimates. */
    int *model_representatives = importance ? postfit_model_representatives(
        gamma_sample, n_draws, n_platforms, n_platform_models_c, n_features,
        cache_bytes, cache_hash_bits) : NULL;
    int n_subjects = 0;
    for (int l = 0; l < n_subgroups; l++)
        n_subjects += sample_size_ptr[l];
    SEXP predictions_R = PROTECT(allocMatrix(REALSXP, n_subjects, n_cv_rounds));
    SEXP folds_R = PROTECT(allocMatrix(INTSXP, n_subjects, n_cv_rounds));
    protect_count += 2;
    double fold_score_sum[n_subgroups + 1], pooled_score[n_subgroups + 1];
    double **round_predictions = malloc((n_subgroups + 1) * sizeof(double *));
    double **round_response = NULL;
    _Bool **round_binary = NULL;
    if (outcome_type == IMR_OUTCOME_BINARY)
        round_binary = malloc((n_subgroups + 1) * sizeof(_Bool *));
    else
        round_response = malloc((n_subgroups + 1) * sizeof(double *));

    _Bool **round_event = malloc((n_subgroups + 1) * sizeof(_Bool *));
    int round_offsets[n_subgroups + 1];
    for (int m = 0; m < n_subgroups + 1; m++)
    {
        if (m < n_subgroups)
        {
            round_predictions[m] = malloc(sample_size_ptr[m] * sizeof(double));
            if (outcome_type == IMR_OUTCOME_BINARY)
                round_binary[m] = malloc(sample_size_ptr[m] * sizeof(_Bool));
            else
                round_response[m] = malloc(sample_size_ptr[m] * sizeof(double));
            round_event[m] = malloc(sample_size_ptr[m] * sizeof(_Bool));
        }
        else
        {
            round_predictions[m] = malloc(n_subjects * sizeof(double));
            if (outcome_type == IMR_OUTCOME_BINARY)
                round_binary[m] = malloc(n_subjects * sizeof(_Bool));
            else
                round_response[m] = malloc(n_subjects * sizeof(double));
            round_event[m] = malloc(n_subjects * sizeof(_Bool));
        }
        fold_score_sum[m] = 0;
        round_offsets[m] = 0;
    }
    int **censored_index = calloc(n_subgroups, sizeof(int *));
    int *n_censored = calloc(n_subgroups, sizeof(int));
    _Bool **event = calloc(n_subgroups, sizeof(_Bool *));
    for (int i = 0; i < n_subgroups; i++)
    {
        censored_index[i] = calloc(sample_size_ptr[i], sizeof(int));
        event[i] = calloc(sample_size_ptr[i], sizeof(_Bool));
        n_censored[i] = 0;

        for (int j = 0; j < sample_size_ptr[i]; j++)
        {
            if (outcome_type == IMR_OUTCOME_SURVIVAL) // for survival outcome
            {
                event[i][j] = (_Bool)response[i][j][1];
            }
            else
            {
                event[i][j] = 1; // for continuous and binary outcome
            }
        }
        find_indices_not_equal(sample_size_ptr[i], event[i], 1, censored_index[i], &n_censored[i]);
    }

    long seed = (long)REAL(seed_R)[0];
    gsl_rng *r = gsl_rng_alloc(gsl_rng_rand48);
    gsl_rng_set(r, seed);
    int n_uncensored[n_subgroups];
    int **uncensored_index = malloc(n_subgroups * sizeof(int *));
    for (int m = 0; m < n_subgroups; m++)
    {
        uncensored_index[m] = malloc((sample_size_ptr[m] - n_censored[m]) * sizeof(int));
    }

    int cv_round;
    double **c_index_list = calloc(n_cv_rounds, sizeof(double *));
    double **total_c_index_list = calloc(n_cv_rounds, sizeof(double *));
    for (cv_round = 0; cv_round < n_cv_rounds; cv_round++)
    {
        c_index_list[cv_round] = calloc((n_subgroups + 1), sizeof(double));
        total_c_index_list[cv_round] = calloc((n_subgroups + 1), sizeof(double));
        for (int m = 0; m < n_subgroups; m++)
        {
            if (n_censored[m] > 0)
                gsl_ran_shuffle(r, censored_index[m], n_censored[m], sizeof(int));
            find_indices_not_equal(sample_size_ptr[m], event[m], 0, uncensored_index[m], &n_uncensored[m]);
            if (n_uncensored[m] > 0)
                gsl_ran_shuffle(r, uncensored_index[m], n_uncensored[m], sizeof(int));
            round_offsets[m] = 0;
            fold_score_sum[m] = 0;
        }
        fold_score_sum[n_subgroups] = 0;
        round_offsets[n_subgroups] = 0;

        for (fold = 0; fold < n_folds; fold++)
        {
            int j = 0;
            double *fold_response = NULL;
            _Bool *fold_binary = NULL;
            if ((outcome_type == IMR_OUTCOME_SURVIVAL) || (outcome_type == IMR_OUTCOME_CONTINUOUS)) // survival outcome or continuous
                fold_response = dvector(0, n_subjects - 1);
            else // binary outcome
                fold_binary = calloc(n_subjects, sizeof(_Bool));

            _Bool *fold_event = malloc(n_subjects * sizeof(_Bool));
            double *fold_prediction = dvector(0, n_subjects - 1);

            for (int m = 0; m < n_subgroups; m++)
            {
                int test_sample_size;
                _Bool *test_binary = NULL;
                double *test_response = NULL;
                _Bool test_event[sample_size_ptr[m]];
                int test_index[sample_size_ptr[m]], train_index[sample_size_ptr[m]];
                postfit_partition(fold, n_folds, sample_size_ptr[m], &test_sample_size, censored_index[m], n_censored[m],
                          uncensored_index[m], test_index, train_index);

                int subgroup_offset = 0;
                for (int g = 0; g < m; ++g) subgroup_offset += sample_size_ptr[g];
                double *prediction = postfit_predict_fold(outcome_type, m, n_covariates, n_model_platforms_c[m], model_platforms_c[m],
                                        n_platform_models_c, platform_models_c, n_features,
                                        sample_size_ptr[m], test_sample_size, test_index, train_index,
                                        latent_response[m], covariates[m], features[m], gamma_sample, residual_shape, residual_rate,
                                        max_models, model_index, high_model_index, n_draws, importance,
                                        model_representatives);

                /* A fold cannot exceed its validated subgroup size. Allocate
                 * that capacity so allocation does not depend on GCC's range
                 * inference for the partition helper's signed output. */
                if (outcome_type == IMR_OUTCOME_BINARY) // binary outcome
                    test_binary = calloc(sample_size_ptr[m], sizeof(_Bool));
                else
                    test_response = calloc(sample_size_ptr[m], sizeof(double));

                for (int i = 0; i < test_sample_size; i++)
                {
                    test_event[i] = event[m][test_index[i]];
                    if (outcome_type == IMR_OUTCOME_BINARY) // binary outcome
                    {
                        test_binary[i] = (_Bool)response[m][test_index[i]][0];
                    }
                    else
                    {
                        test_response[i] = response[m][test_index[i]][0];
                    }
                }
                double ci;
                if (outcome_type == IMR_OUTCOME_SURVIVAL) // survival outcome
                    ci = legacy_concordance(test_sample_size, prediction, test_response, test_event);
                else if (outcome_type == IMR_OUTCOME_BINARY)                      // binary outcome
                    ci = legacy_auc(test_sample_size, prediction, test_binary); // this is the AUC
                else
                    ci = legacy_mse(test_sample_size, prediction, test_response); // continuous outcome

                fold_score_sum[m] += ci;
                for (int i = 0; i < test_sample_size; i++)
                {

                    if (outcome_type == IMR_OUTCOME_BINARY)
                    {
                        round_binary[m][i + round_offsets[m]] = test_binary[i];
                        fold_binary[i + j] = test_binary[i];
                    }
                    else
                    {
                        fold_response[i + j] = test_response[i];
                        round_response[m][i + round_offsets[m]] = test_response[i];
                    }
                    round_event[m][i + round_offsets[m]] = test_event[i];
                    fold_event[i + j] = test_event[i];
                    fold_prediction[i + j] = prediction[i];
                    round_predictions[m][i + round_offsets[m]] = prediction[i];
                }
                for (int i = 0; i < test_sample_size; ++i) {
                    int row = subgroup_offset + test_index[i] + n_subjects * cv_round;
                    REAL(predictions_R)[row] = prediction[i];
                    INTEGER(folds_R)[row] = fold + 1;
                }
                j += test_sample_size;
                round_offsets[m] += test_sample_size;
                if (outcome_type == IMR_OUTCOME_BINARY)
                    free(test_binary);
                else
                    free(test_response);
                free(prediction);
            }
            for (int i = 0; i < j; i++)
            {
                round_predictions[n_subgroups][i + round_offsets[n_subgroups]] = fold_prediction[i];
                if (outcome_type == IMR_OUTCOME_BINARY)
                {
                    round_binary[n_subgroups][i + round_offsets[n_subgroups]] = fold_binary[i];
                }
                else
                    round_response[n_subgroups][i + round_offsets[n_subgroups]] = fold_response[i];

                round_event[n_subgroups][i + round_offsets[n_subgroups]] = fold_event[i];
            }
            round_offsets[n_subgroups] += j;
            double ci4;

            if (outcome_type == IMR_OUTCOME_SURVIVAL) // survival outcome
                ci4 = legacy_concordance(j, fold_prediction, fold_response, fold_event);
            else if (outcome_type == IMR_OUTCOME_BINARY)                // binary outcome
                ci4 = legacy_auc(j, fold_prediction, fold_binary); // historical binary score
            else
                ci4 = legacy_mse(j, fold_prediction, fold_response); // continuous outcome

            fold_score_sum[n_subgroups] += ci4;
            if (outcome_type == IMR_OUTCOME_BINARY)
                free(fold_binary);
            else
                free(fold_response);

            free(fold_event);
            free(fold_prediction);
        }
        Rprintf("\nCV round %d/%d (%d folds)\n", cv_round + 1, n_cv_rounds, n_folds);
        const char *metric_name[] = {"C-index", "AUC", "MSE"};
        for (int m = 0; m < n_subgroups + 1; m++)
        {
            fold_score_sum[m] = fold_score_sum[m] / n_folds;
            if (m < n_subgroups)
            {
                Rprintf("Average of %s from %d-CVfold of group %d is  %f\n", metric_name[outcome_type - 1], n_folds, m + 1, fold_score_sum[m]);
            }
            else
            {
                Rprintf("Average of %s from %d-CVfold of all samples  is  %f\n", metric_name[outcome_type - 1], n_folds, fold_score_sum[m]);
            }
            int subset_size;
            if (m == n_subgroups)
                subset_size = n_subjects;
            else
                subset_size = sample_size_ptr[m];
            if (outcome_type == IMR_OUTCOME_SURVIVAL) // survival outcome
                pooled_score[m] = legacy_concordance(subset_size, round_predictions[m], round_response[m], round_event[m]);
            else if (outcome_type == IMR_OUTCOME_BINARY)                      // binary outcome
                pooled_score[m] = legacy_auc(subset_size, round_predictions[m], round_binary[m]);
            else if (outcome_type == IMR_OUTCOME_CONTINUOUS)                      // continuous outcome
                pooled_score[m] = legacy_mse(subset_size, round_predictions[m], round_response[m]);
            if (m < n_subgroups)
                Rprintf("\nTotal %s from group %d is %f\n", metric_name[outcome_type - 1], m + 1, pooled_score[m]);
            else
                Rprintf("\nTotal %s from all samples is %f\n", metric_name[outcome_type - 1], pooled_score[m]);
            c_index_list[cv_round][m] = fold_score_sum[m];
            total_c_index_list[cv_round][m] = pooled_score[m];
        }

    }

    // Convert to R outputs
    SEXP total_c_index_R = PROTECT(c_array_to_r_matrix(total_c_index_list, n_cv_rounds, n_subgroups + 1));
    protect_count++;
    SEXP c_index_R = PROTECT(c_array_to_r_matrix(c_index_list, n_cv_rounds, n_subgroups + 1));
    protect_count++;
    int list_size = 4;
    SEXP list;
    SEXP list_names;
    PROTECT(list = allocVector(VECSXP, list_size));
    protect_count++;
    SET_VECTOR_ELT(list, 0, total_c_index_R);
    SET_VECTOR_ELT(list, 1, c_index_R);
    SET_VECTOR_ELT(list, 2, predictions_R);
    SET_VECTOR_ELT(list, 3, folds_R);
    PROTECT(list_names = allocVector(STRSXP, list_size));
    protect_count++;
    SET_STRING_ELT(list_names, 0, mkChar("total_cindex"));
    SET_STRING_ELT(list_names, 1, mkChar("subset_cindex"));
    SET_STRING_ELT(list_names, 2, mkChar("predictions"));
    SET_STRING_ELT(list_names, 3, mkChar("folds"));
    setAttrib(list, R_NamesSymbol, list_names);

    /* We free allocated memories */
    free_r_list_list_matrix_to_c(features, n_subgroups, features_R);
    features = NULL;
    for (int m = 0; m < n_subgroups; m++)
    {
        free(latent_response[m]);
        free_dmatrix(covariates[m], 0, sample_size_ptr[m] - 1, 0, n_covariates - 1);
        free_dmatrix(response[m], 0, sample_size_ptr[m] - 1, 0, 2 - 1);
    }
    free(covariates);
    free(response);
    free(latent_response);
    for (int m = 0; m < n_subgroups; m++)
    {
        free(censored_index[m]);
        free(uncensored_index[m]);
        free(event[m]);
    }
    free(n_censored);
    free(censored_index);
    free(uncensored_index);
    free(event);
    gsl_rng_free(r);
    for (int m = 0; m < n_subgroups + 1; m++)
    {
        free(round_predictions[m]);
        if (outcome_type == IMR_OUTCOME_BINARY)
            free(round_binary[m]);
        else
            free(round_response[m]);
        free(round_event[m]);
    }
    free(round_predictions);
    if (outcome_type == IMR_OUTCOME_BINARY)
        free(round_binary);
    else
        free(round_response);
    free(round_event);
    for (int s = 0; s < n_unique_models; s++)
    {
        for (int m = 0; m < n_subgroups; m++)
        {
            free(beta[s][m]);
        }
        free(beta[s]);
    }
    free(beta);

    for (int l = 0; l < n_platforms; l++)
    {
        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            free(theta[l][m]);
            free(interaction_rates[l][m]);
        }
        free(theta[l]);
        free(interaction_rates[l]);
    }
    free(theta);
    free(interaction_rates);
    for (int s = 0; s < n_draws; s++)
    {
        for (int l = 0; l < n_platforms; l++)
        {
            for (int m = 0; m < n_platform_models_c[l]; m++)
            {
                free(gamma_sample[s][l][m]);
            }
            free(gamma_sample[s][l]);
        }
        free(gamma_sample[s]);
    }
    free(gamma_sample);
    free(model_index);
    free(high_model_index);
    free(slab_scales);
    free(post);
    free(platform_models_c);
    free(n_platform_models_c);
    free(model_platforms_c);
    free(n_model_platforms_c);
    for (cv_round = 0; cv_round < n_cv_rounds; cv_round++)
    {
        free(c_index_list[cv_round]);
        free(total_c_index_list[cv_round]);
    }
    free(c_index_list);
    free(total_c_index_list);
    free(mrf);

    started = clock() - started;
    double elapsed = ((double)started) / CLOCKS_PER_SEC; // in seconds
    Rprintf("\n\nTime taken for assessing prediction in seconds is %f\n", elapsed);
    Rprintf("\nTime taken for assessing prediction in minutes is %f\n", elapsed / 60);
    Rprintf("\nTime taken for assessing prediction in hours is %f\n", elapsed / 3600);
    /* Printing can allocate through R's output connection. Keep the result
     * protected until the last allocation-capable operation is finished. */
    UNPROTECT(protect_count);
    return list;
}
