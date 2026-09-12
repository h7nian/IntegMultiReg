#include <stdbool.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include "my_header.h"
#include "utils.h"
#include <Rmath.h>

/*
 * Predict all samples for one subgroup by Bayesian model averaging over the
 * posterior model list inferred from MCMC gamma samples.
 */
double *predict_bma(int model, int K, int n_selected_platforms, int sample_size, int *selected_platforms,
                 int *n_platform_models, int **platform_models, int *G,
                 double **C, double ***X, _Bool ****gamma_sample,
                 double ***beta, double *post, int max_models,
                 int *model_index, int *high_model_index, int n_draws, int outcome_type)
{
  int l, i, j;
  double *yhat = dvector(0, sample_size - 1);
  for (i = 0; i < sample_size; i++)
    yhat[i] = 0;
  for (l = 0; l < max_models; l++)
  {
    int l0 = high_model_index[l];
    int l1 = model_index[l0];

    int *selected_feature_index[n_selected_platforms];
    int n_selected_features[n_selected_platforms];
    for (int i = 0; i < n_selected_platforms; i++)
    {
      int platform_index = selected_platforms[i];
      int platform_model_index = -1;
      for (int ss = 0; ss < n_platform_models[platform_index]; ss++)
      {
        if (platform_models[platform_index][ss] == model)
        {
          platform_model_index = ss;
          break;
        }
      }
      if (platform_model_index == -1)
      {
        Rf_error("Subgroup %d not found for platform %d\n", 1 + model, 1 + platform_index);
      }
      selected_feature_index[i] = malloc(G[platform_index] * sizeof(int));
      if (!selected_feature_index[i])
      {
        Rf_error("malloc failed for selected_feature_index[%d]\n", i);
      }
      n_selected_features[i] = 0;
      find_indices_not_equal(G[platform_index], gamma_sample[n_draws - 1 - l1][platform_index][platform_model_index], 0, selected_feature_index[i], &n_selected_features[i]);
    }

    double **PG = build_design_matrix(K, n_selected_platforms, n_selected_features, selected_feature_index, C, X, selected_platforms, sample_size);
    int total_nx = 0;
    for (int p = 0; p < n_selected_platforms; p++)
    {
      total_nx += n_selected_features[p];
    }
    int k = 1 + K + total_nx;
    for (i = 0; i < sample_size; i++)
    {
      double yh = 0;
      for (j = 0; j < k; j++)
      {
        if (j == 0)
          yh += beta[l0][model][0];
        else
          yh += PG[i][j - 1] * beta[l0][model][j];
      }
      yhat[i] += (outcome_type == IMR_OUTCOME_BINARY ?
                  pnorm5(yh, 0.0, 1.0, 1, 0) : yh) * post[l];
    }
    for (i = 0; i < sample_size; i++)
      free(PG[i]);
    free(PG);
    for (int i = 0; i < n_selected_platforms; i++)
      free(selected_feature_index[i]);
  }
  return yhat;
}

/*
 * Collapse duplicated gamma samples, keep the highest-posterior models, and
 * compute regression coefficients used by prediction and cross-validation.
 */
double ***infer_posterior_models(double **y, double ***C, double ****X, int n_draws, _Bool ****gamma_sample,
                          double *nu, double ***theta, double *mrf, double *h, double h1, double h0,
                          double hg, double alpha0,
                          double alpha, double psi, int *G, int n_subgroups, int n_platforms,
                          int *n_platform_models_c, int *n_model_platforms_c, int **model_platforms_c,
                          int **platform_models_c, int *sample_size_ptr, int K,
                          double ***betaTh, const char *likelihood_type, double *post, int *model_index,
                          int *high_model_index, int *n_unique_models_out, int max_models)

{

  int i, j;
  int n_unique_models = 0;
  model_index[n_unique_models] = 0;
  for (i = 0; i < n_draws; i++)
  {
    for (j = 0; j < n_unique_models; j++)
    {
      int model_differs = 0;
      for (int l = 0; l < n_platforms; l++)
      {
        for (int k = 0; k < n_platform_models_c[l]; k++)
        {
          if (!bool_vectors_equal(G[l], gamma_sample[n_draws - 1 - i][l][k], gamma_sample[n_draws - 1 - model_index[j]][l][k]))
          {
            model_differs = 1;
            break;
          }
        }
        if (model_differs)
        {
          break; // breaks outer loop
        }
      }
      if (!model_differs)
        break; // The complete selection state matches this representative.

    } // end of the j loop

    if (j == n_unique_models)
    {
      model_index[n_unique_models] = i;
      n_unique_models++;
    }
    if (n_unique_models == 100 * max_models)
    { // we limit the number of models to maximum of 100 times maxmodels
      break;
    }
  } // end of the i loop
  *n_unique_models_out = n_unique_models;
  Rprintf("\n\n\nNumber of different gamma values for variable selection (nbr of models) = %d", n_unique_models);
  Rprintf("\n");
  double ***beta = malloc(n_unique_models * sizeof(double **));
  int l;
  for (l = 0; l < n_unique_models; l++)
  {
    beta[l] = malloc(n_subgroups * sizeof(double *));
    int l1 = model_index[l];
    double loglik[n_subgroups];
    for (int m = 0; m < n_subgroups; m++)
    {
      int N = sample_size_ptr[m];
      int **selected_feature_index = calloc(n_model_platforms_c[m], sizeof(int *));
      int *n_selected_features = calloc(n_model_platforms_c[m], sizeof(int));
      for (int ll = 0; ll < n_model_platforms_c[m]; ll++)
      {
        int platform_index = model_platforms_c[m][ll];
        int platform_model_index = -1;
        for (int ss = 0; ss < n_platform_models_c[platform_index]; ss++)
        {
	          if (platform_models_c[platform_index][ss] == m)
	          {
	            platform_model_index = ss;
	            break;
	          }
	        }
	        if (platform_model_index == -1)
	        {
	          Rf_error("Subgroup %d not found for platform %d\n", m, platform_index);
	        }
        selected_feature_index[ll] = calloc(G[platform_index], sizeof(int));
        find_indices_not_equal(G[platform_index], gamma_sample[n_draws - 1 - l1][platform_index][platform_model_index], 0, selected_feature_index[ll], &n_selected_features[ll]);
      }
      double **PG = build_design_matrix(K, n_model_platforms_c[m], n_selected_features, selected_feature_index,
                            C[m], X[m], model_platforms_c[m], N);
      int total_selected_features = 0;
      for (int ll = 0; ll < n_model_platforms_c[m]; ll++)
      {
	        total_selected_features += n_selected_features[ll];
	      }

      int s;

      if (strcmp(likelihood_type, "Local") == 0)
      {
        double a;
        double *Sigma = malloc((size_t) N * N * sizeof(double));
        if (!Sigma) Rf_error("malloc failed for Sigma");
        for (int i = 0; i < N; i++)
        {
          for (int j = 0; j <= i; j++)
          {
            a = 0;
            for (s = 0; s < total_selected_features; s++)
            {
              a += PG[i][s] * PG[j][s];
            }
            Sigma[i * N + j] = Sigma[j * N + i] = h0 + h[m] * a;
          }
          Sigma[i * N + i] += 1;
        }
        double logdet = 0;
        double scal = cholesky_quadratic_form(N, Sigma, y[m], &logdet);
        loglik[m] = -(N / 2.0) * log(IMR_PI * 2 * alpha) + gsl_sf_lngamma(N / 2.0 + alpha) - gsl_sf_lngamma(alpha) - 0.5 * logdet - ((N / 2.0) + alpha) * log(1 + scal / (2 * psi));
        free(Sigma);
      }
      else
      {
        int maxiter = 40;
        double stop = 1e-3;
        int rr = 1;
        int k = 1 + K + total_selected_features;
        double *precision = build_posterior_precision(k, K, n_selected_features[0], N, h[m], h1, h0, hg, PG);
        double *precision_copy = malloc((size_t) k * k * sizeof(double));
        if (!precision_copy) Rf_error("malloc failed for precision_copy");
        for (int i = 0; i < k; i++)
        {
          for (int j = 0; j <= i; j++)
          {
	            precision_copy[i * k + j] = precision_copy[j * k + i] = precision[i * k + j];
	          }
	        }
	        gsl_matrix_view m11 = gsl_matrix_view_array(precision, k, k);
        gsl_linalg_cholesky_decomp(&m11.matrix);
        beta[l][m] = malloc(k * sizeof(double));
	        loglik[m] = log_likelihood_nonlocal(k, K, n_selected_features[0], N, alpha, psi, y[m], PG, precision_copy, &m11.matrix,
	                                   beta[l][m], rr, h[m], h1, h0, hg, maxiter, stop, 0);
	        free(precision);
	        free(precision_copy);
	      }
      for (int i = 0; i < N; i++)
        free(PG[i]);
      free(PG);
      free(n_selected_features);
      for (int ll = 0; ll < n_model_platforms_c[m]; ll++)
        free(selected_feature_index[ll]);
      free(selected_feature_index);
    }
    post[l] = log_posterior(loglik, gamma_sample[n_draws - 1 - l1], nu, theta, mrf, alpha0, betaTh, n_subgroups,
                      n_platforms, G, n_platform_models_c);
  } // end number of models l=0

  sort_descending_index(n_unique_models, post, high_model_index);
  double max_log_post = post[0];

  int n_top_models = MIN(max_models, n_unique_models);
  for (l = 0; l < n_top_models; l++)
  {
    post[l] = exp(post[l] - max_log_post);
  }
  double sum_post = n_top_models * mean(n_top_models, post);

	for (l = 0; l < n_top_models; l++)
	{
	    post[l] = post[l] / sum_post;
	}
  return beta;
}
