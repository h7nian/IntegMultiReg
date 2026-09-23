#include <stdbool.h>
#include <stdint.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include "my_header.h"
#include "utils.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

static double ****X1 = NULL;
static double ***newYY = NULL;
static double ***newCC = NULL;

/* Export every draw in its original position, sharing only exactly equal
 * platform matrices. Hash collisions are resolved by full comparison; neither
 * sampler state nor draw multiplicity changes. R owns the temporary hash table
 * so it is reclaimed even if a subsequent R allocation raises an error. */
static SEXP export_selection_history(_Bool ****samples, int draws, int platforms,
                                     const int *rows, const int *columns)
{
    SEXP result = PROTECT(allocVector(VECSXP, draws));
    for (int s = 0; s < draws; ++s) {
        SEXP state = PROTECT(allocVector(VECSXP, platforms));
        SET_VECTOR_ELT(result, s, state);
        UNPROTECT(1);
    }
    size_t capacity = 1;
    while (capacity < (size_t)draws * 2 && capacity <= SIZE_MAX / 2)
        capacity *= 2;
    int *table = (int *)R_alloc(capacity, sizeof(int));
    for (int p = 0; p < platforms; ++p) {
        memset(table, 0, capacity * sizeof(int));
        for (int s = 0; s < draws; ++s) {
            size_t slot = 0;
            int previous = -1;
            if (table) {
                uint64_t hash = UINT64_C(14695981039346656037);
                for (int r = 0; r < rows[p]; ++r)
                    for (int c = 0; c < columns[p]; ++c) {
                        hash ^= (uint64_t)samples[s][p][r][c];
                        hash *= UINT64_C(1099511628211);
                    }
                slot = (size_t)hash & (capacity - 1);
                while (table[slot]) {
                    int candidate = table[slot] - 1, equal = 1;
                    for (int r = 0; r < rows[p] && equal; ++r)
                        for (int c = 0; c < columns[p]; ++c)
                            if (samples[s][p][r][c] != samples[candidate][p][r][c]) {
                                equal = 0;
                                break;
                            }
                    if (equal) { previous = candidate; break; }
                    slot = (slot + 1) & (capacity - 1);
                }
            }
            if (previous >= 0) {
                SEXP matrix = VECTOR_ELT(VECTOR_ELT(result, previous), p);
                MARK_NOT_MUTABLE(matrix); /* preserve R copy-on-modify semantics */
                SET_VECTOR_ELT(VECTOR_ELT(result, s), p, matrix);
            } else {
                SEXP matrix = PROTECT(c_array_to_r_matrix_int(samples[s][p], rows[p], columns[p]));
                SET_VECTOR_ELT(VECTOR_ELT(result, s), p, matrix);
                UNPROTECT(1);
                if (table) table[slot] = s + 1;
            }
        }
    }
    UNPROTECT(1);
    return result;
}

/*
 * Main training entry point called from R.
 *
 * The routine converts nested R lists into row-addressable C arrays, initializes
 * the latent responses and selection state, runs the MCMC updates, and returns
 * posterior samples/summaries in the historical R list layout.
 */
SEXP imr_fit(SEXP h0_R, SEXP hh_R, SEXP alpha_R, SEXP psi_R, SEXP alpha0_R, SEXP beta0_R,
                  SEXP seed_R, SEXP nu_R,
                  SEXP method1_R, SEXP n_platforms_R,
                  SEXP platform_models_R, SEXP model_platforms_R, SEXP n_subgroups_R,
                  SEXP sample_size, SEXP n_features_R, SEXP n_covariates_R,
                  SEXP X1_filtered, SEXP newYY_list, SEXP type_outcome,
                  SEXP newCC_list,
                  SEXP draws_R, SEXP burnin_R, SEXP sampler_method_R, SEXP numerical_R, SEXP initial_R)
{
    const imr_numerical_control numerical = imr_read_numerical_control(numerical_R);
    if (!isInteger(sampler_method_R) || XLENGTH(sampler_method_R) != 1 ||
        INTEGER(sampler_method_R)[0] < IMR_SAMPLER_LEGACY ||
        INTEGER(sampler_method_R)[0] > IMR_SAMPLER_PAPER)
        Rf_error("Invalid sampler method");
    int sampler_method = INTEGER(sampler_method_R)[0];
    /* Reject malformed explicit starts before allocating native workspaces. */
    if (initial_R != R_NilValue) {
        int np = asInteger(n_platforms_R);
        if (TYPEOF(initial_R) != VECSXP || XLENGTH(initial_R) != 2)
            Rf_error("Invalid initial state");
        for (int component = 0; component < 2; ++component) {
            SEXP values = VECTOR_ELT(initial_R, component);
            if (values == R_NilValue) continue;
            if (component == 1 && strcmp(CHAR(STRING_ELT(method1_R, 0)), "BMS") == 0)
                Rf_error("Initial interactions are not applicable to BMS");
            if (TYPEOF(values) != VECSXP || XLENGTH(values) != np)
                Rf_error("Invalid initial platform list");
            for (int p = 0; p < np; ++p) {
                int nr = LENGTH(VECTOR_ELT(platform_models_R, p));
                int nc = component ? nr : INTEGER(n_features_R)[p];
                SEXP x = VECTOR_ELT(values, p), dim = getAttrib(x, R_DimSymbol);
                if (TYPEOF(x) != (component ? REALSXP : INTSXP) ||
                    TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
                    INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc)
                    Rf_error("Invalid initial matrix dimensions or type");
                for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i) {
                    R_xlen_t index = i + (R_xlen_t)nr * j;
                    if (!component) {
                        if (INTEGER(x)[index] != 0 && INTEGER(x)[index] != 1)
                            Rf_error("Invalid initial selection value");
                    } else {
                        double v = REAL(x)[index];
                        if (!R_FINITE(v) || (i == j ? v != 0 : v <= 0) ||
                            v != REAL(x)[j + (R_xlen_t)nr * i])
                            Rf_error("Invalid initial interaction value");
                    }
                }
            }
        }
    }
    clock_t t = clock();
    int protect_count = 0;
    /* platform_models_R maps each platform to the subgroups using it;
     * model_platforms_R is the inverse mapping from subgroup to platforms. */

    double h0 = REAL(h0_R)[0]; // scaling factor
    double h11 = REAL(hh_R)[0];
    double alpha = REAL(alpha_R)[0];   // weight of prior beliefs
    double psi = REAL(psi_R)[0];       // control var of prior distributions
    double alpha0 = REAL(alpha0_R)[0]; // prior
    double beta0 = REAL(beta0_R)[0];   // prior scaling factor

    long seed = (long)REAL(seed_R)[0];

    PROTECT(method1_R);
    protect_count++;
    const char *model_method = CHAR(STRING_ELT(method1_R, 0));

    PROTECT(n_platforms_R);
    protect_count++;
    PROTECT(platform_models_R);
    protect_count++;
    PROTECT(model_platforms_R);
    protect_count++;
    PROTECT(n_subgroups_R);
    protect_count++;
    PROTECT(sample_size);
    protect_count++;
    PROTECT(n_features_R);
    protect_count++;
    PROTECT(n_covariates_R);
    protect_count++;
    PROTECT(X1_filtered);
    protect_count++;
    PROTECT(newYY_list);
    protect_count++;
    PROTECT(newCC_list);
    protect_count++;


    int K = asInteger(n_covariates_R);
    int outcome_type = asInteger(type_outcome);
    int n_draws = asInteger(draws_R);
    int n_burnin = asInteger(burnin_R);

    int n_platforms = asInteger(n_platforms_R);
    Rprintf("We have %d platforms  in total \n", n_platforms);

    // We read model indices for each plaform
    int **platform_models_c = malloc(n_platforms * sizeof(int *));
    int *n_platform_models_c = malloc(n_platforms * sizeof(int));

    for (int i = 0; i < n_platforms; i++)
    {
        SEXP mPM = VECTOR_ELT(platform_models_R, i);
        int sizeM = LENGTH(mPM);
        n_platform_models_c[i] = sizeM;
        platform_models_c[i] = INTEGER(mPM);
    }

    Rprintf("\n");
    for (int i = 0; i < n_platforms; i++)
    {
        Rprintf("Platform %d is involved in  %d subgroups\n", i + 1, n_platform_models_c[i]);
        Rprintf("Platform %d is involved in subgroups: ", i + 1);
        for (int j = 0; j < n_platform_models_c[i]; j++)
        {
            Rprintf("%d ", 1 + platform_models_c[i][j]);
        }
        Rprintf("\n\n");
    }

    int n_subgroups = asInteger(n_subgroups_R);
    // We read platform indices for each model/subgroup
    int **model_platforms_c = malloc(n_subgroups * sizeof(int *));
    int *n_model_platforms_c = malloc(n_subgroups * sizeof(int));

    for (int i = 0; i < n_subgroups; i++)
    {
        SEXP mPM = VECTOR_ELT(model_platforms_R, i);
        int sizeM = LENGTH(mPM);
        n_model_platforms_c[i] = sizeM;
        model_platforms_c[i] = INTEGER(mPM);
    }

    Rprintf("\n");

    for (int i = 0; i < n_subgroups; i++)
    {
        Rprintf("\nNumber of platforms for subgroup %d is %d: \n\n", i + 1, n_model_platforms_c[i]);
        Rprintf("Platforms for subgroup %d are: ", i + 1);
        for (int j = 0; j < n_model_platforms_c[i]; j++)
        {
            Rprintf("Platform/View %d, ", 1 + model_platforms_c[i][j]);
        }
        Rprintf("\n");
    }

    int *sample_size_ptr = INTEGER(sample_size);

    Rprintf("\nSample sizes for each selected subgroup:\n");
    for (int i = 0; i < n_subgroups; i++)
    {
        Rprintf(" %d", sample_size_ptr[i]);
    }
    Rprintf("\n");

    int *G = INTEGER(n_features_R);
    Rprintf("Number of features for each platform:\n");
    for (int i = 0; i < n_platforms; i++)
    {
        Rprintf(" %d", G[i]);
    }
    Rprintf("\n");

    double ****X0 = r_list_list_matrix_to_c(n_subgroups, X1_filtered);

    X1 = X0;
    for (int i = 0; i < n_subgroups; i++)
    {
        SEXP subgroup = VECTOR_ELT(X1_filtered, i);
        int n_platforms = LENGTH(subgroup);
        Rprintf("Subgroup %d: %d platforms\n", i + 1, n_platforms);
        for (int j = 0; j < n_platforms; j++)
        {
            SEXP df = VECTOR_ELT(subgroup, j);
            SEXP dims = getAttrib(df, R_DimSymbol);
            int n_rows = INTEGER(dims)[0];
            int n_cols = INTEGER(dims)[1];
            Rprintf("  Platform %d: %d rows, %d cols\n", j + 1, n_rows, n_cols);
        }
    }

    /* Read outcome and covariate arrays from R list storage. */
    double ***newYY_arr = r_list_matrix_to_c(n_subgroups, newYY_list);

    newYY = newYY_arr;

    double ***newCC_ptrs = r_list_matrix_to_c(n_subgroups, newCC_list);
    newCC = newCC_ptrs;
    const char likelihood_type[] = "NonLocal";

    double *h = malloc(n_subgroups * sizeof(double));
    for (int i = 0; i < n_subgroups; i++)
    {
        h[i] = h11;
    }
    double h1 = h0, hg = h11;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_rand48);
    gsl_rng_set(r, seed);

    // Initialize censor index, acceptance ratios, and other variables
    int **censored_index = malloc(n_subgroups * sizeof(int *));
    double **ylatent = malloc(n_subgroups * sizeof(double *));

    for (int m = 0; m < n_subgroups; m++)
    {
        censored_index[m] = malloc(sample_size_ptr[m] * sizeof(int));
        ylatent[m] = dvector(0, sample_size_ptr[m] - 1);
    }
    int n_censored[n_subgroups];
    double **ymean = malloc(n_subgroups * sizeof(double *));
    double **yobs = NULL;
    _Bool **yobsb = NULL;
    if (outcome_type == IMR_OUTCOME_BINARY)
    {
        yobsb = malloc(n_subgroups * sizeof(_Bool *));
    }
    else // survival or continuous
    {
        yobs = malloc(n_subgroups * sizeof(double *));
    }

    for (int i = 0; i < n_subgroups; i++)
    {
        n_censored[i] = 0;
        ymean[i] = dvector(0, sample_size_ptr[i] - 1);
        if ((outcome_type == IMR_OUTCOME_SURVIVAL) ||
            (outcome_type == IMR_OUTCOME_CONTINUOUS))
            yobs[i] = dvector(0, sample_size_ptr[i] - 1);
        else // binary
            yobsb[i] = (_Bool *)malloc(sample_size_ptr[i] * sizeof(_Bool));

        _Bool Delta[sample_size_ptr[i]];
        for (int j = 0; j < sample_size_ptr[i]; j++)
        {
            if (outcome_type == IMR_OUTCOME_SURVIVAL)
                Delta[j] = (_Bool)newYY_arr[i][j][1];
            else
                Delta[j] = 1;
            if ((outcome_type == IMR_OUTCOME_SURVIVAL) ||
                (outcome_type == IMR_OUTCOME_CONTINUOUS))
                ymean[i][j] = ylatent[i][j] =yobs[i][j] = newYY_arr[i][j][0];
            else //binary
                ylatent[i][j] = yobsb[i][j] = (_Bool)newYY_arr[i][j][0];
            // To check for binary
           // ymean[i][j] = ylatent[i][j] = yobs[i][j];
            // ylatent[i][j] = gsl_ran_exponential(r, ylatent[i]);

            if ((Delta[j] == 0) && (outcome_type == IMR_OUTCOME_SURVIVAL))
            {
                ylatent[i][j] += 0.01;
                ymean[i][j] = 0;
            }
             if (outcome_type == IMR_OUTCOME_BINARY)
                ymean[i][j] = 0;
            
        }
        find_indices_not_equal(sample_size_ptr[i], Delta, 1, censored_index[i], &n_censored[i]);
        if (outcome_type == IMR_OUTCOME_SURVIVAL)
            Rprintf("\nNumber of censored values for subgroup %d is %d\n", i + 1, n_censored[i]);
    }
    // Memory for MCMC acceptance (only needed for censored or binary latent
    // updates; for continuous outcomes accept_y is left NULL).
    double **accept_y = NULL;
    if ((outcome_type == IMR_OUTCOME_SURVIVAL) ||
        (outcome_type == IMR_OUTCOME_BINARY))
    {
        accept_y = malloc(n_subgroups * sizeof(double *));
        for (int m = 0; m < n_subgroups; m++)
        {
            accept_y[m] = calloc(sample_size_ptr[m], sizeof(double));
        }
    }

    double ***theta = calloc(n_platforms, sizeof(double **));
    double ***betaTh = calloc(n_platforms, sizeof(double **));

    double ***accept_theta = calloc(n_platforms, sizeof(double **));
    double **accept_gamma = calloc(n_platforms, sizeof(double *));
    _Bool ***gamma = calloc(n_platforms, sizeof(_Bool **));
    double ***gamma_mean = calloc(n_platforms, sizeof(double **));

    for (int l = 0; l < n_platforms; l++)
    {
        theta[l] = calloc(n_platform_models_c[l], sizeof(double *));
        betaTh[l] = calloc(n_platform_models_c[l], sizeof(double *));
        accept_theta[l] = calloc(n_platform_models_c[l], sizeof(double *));
        accept_gamma[l] = calloc(n_platform_models_c[l], sizeof(double));
        gamma[l] = calloc(n_platform_models_c[l], sizeof(_Bool *));
        gamma_mean[l] = calloc(n_platform_models_c[l], sizeof(double *));
        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            gamma[l][m] = calloc(G[l], sizeof(_Bool));
            gamma_mean[l][m] = calloc(G[l], sizeof(double));
            accept_theta[l][m] = calloc(n_platform_models_c[l], sizeof(double));
            theta[l][m] = calloc(n_platform_models_c[l], sizeof(double));
            betaTh[l][m] = calloc(n_platform_models_c[l], sizeof(double));
            for (int m1 = 0; m1 < n_platform_models_c[l]; m1++)
            {
                betaTh[l][m][m1] = beta0;
            }
        }
    }
    // double *nu=calloc(n_platforms,sizeof(double*));
    double *nu = REAL(nu_R);
    for (int l = 0; l < n_platforms; l++)
    {
        /* Evaluate the logistic inclusion probability without forming
         * exp(nu) / (1 + exp(nu)), which becomes Inf / Inf for a large
         * positive prior and is undefined when converted to an integer. */
        double exp_term = exp(nu[l] < 0.0 ? nu[l] : -nu[l]);
        double inclusion_probability = nu[l] < 0.0
            ? exp_term / (1.0 + exp_term)
            : 1.0 / (1.0 + exp_term);
        int initial_inclusions = (int)(inclusion_probability * G[l]);
        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            SEXP selection = initial_R == R_NilValue ? R_NilValue : VECTOR_ELT(initial_R, 0);
            if (selection != R_NilValue) {
                SEXP values = VECTOR_ELT(selection, l);
                for (int j = 0; j < G[l]; ++j)
                    gamma[l][m][j] = INTEGER(values)[m + (R_xlen_t)n_platform_models_c[l] * j];
            }
            for (int i = 0; selection == R_NilValue && i < initial_inclusions; i++)
            {
                int ii = gsl_rng_uniform_int(r, G[l]);
                gamma[l][m][ii] = 1;
            }
            for (int m1 = 0; m1 < n_platform_models_c[l]; m1++)
            {
                if (strcmp(model_method, "BMS") != 0)
                {
                    SEXP interaction = initial_R == R_NilValue ? R_NilValue : VECTOR_ELT(initial_R, 1);
                    theta[l][m][m1] = interaction == R_NilValue ? 0.1 * (m != m1) :
                        REAL(VECTOR_ELT(interaction, l))[m + (R_xlen_t)n_platform_models_c[l] * m1];
                }
                else
                {
                    theta[l][m][m1] = 0;
                }
            }
        }
    }

    // int model;
    //srand(seed);
    GetRNGstate();
    double log_likelihood[n_subgroups], logdet[n_subgroups], scal[n_subgroups];
    Rprintf("\n");

    initialize_sampler_state(outcome_type, ylatent, newCC, X1,
                gamma, n_platforms, G, n_subgroups,
                platform_models_c, n_platform_models_c,
                model_platforms_c,
                n_model_platforms_c, sample_size_ptr,
                log_likelihood, logdet, scal, h, h1, h0, hg, alpha, psi, K, &numerical);

    double *mrf = calloc(n_platforms, sizeof(double));
    for (int l = 0; l < n_platforms; l++)
    {
        compute_mrf_log_normalizer(n_platform_models_c[l], theta[l], nu[l], &mrf[l]);
        Rprintf("%.3lf ", mrf[l]);
    }

    // int s, su, su1;

    double *log_posterior_sample = dvector(0, n_burnin + n_draws - 1);
    _Bool ****gamma_sample = malloc(n_draws * sizeof(_Bool ***));

    for (int s = 0; s < n_draws; s++)
    {
        gamma_sample[s] = malloc(n_platforms * sizeof(_Bool **)); // n_platforms
        for (int l = 0; l < n_platforms; l++)
        {
            gamma_sample[s][l] = bmatrix(0, n_platform_models_c[l] - 1, 0, G[l] - 1);
        }
    }

    if (!gamma_sample)
    {
        nrerror("allocation failure; take appropriate action");
    }

    double ***theta_sample = malloc(n_platforms * sizeof(double **));
    for (int l = 0; l < n_platforms; l++)
    {
        int n_theta_pairs = n_platform_models_c[l] *
                            (n_platform_models_c[l] - 1) / 2;
        theta_sample[l] = n_theta_pairs > 0
            ? dmatrix(0, n_draws - 1, 0, n_theta_pairs - 1)
            : NULL;
    }
    _Bool thetaFreed = false;
    if (strcmp(model_method, "BMS") == 0)
    {
        for (int l = 0; l < n_platforms; l++)
        {
            int n_theta_pairs = n_platform_models_c[l] *
                                (n_platform_models_c[l] - 1) / 2;
            if (n_theta_pairs > 0)
                free_dmatrix(theta_sample[l], 0, n_draws - 1,
                             0, n_theta_pairs - 1);
        }
        free(theta_sample);
        thetaFreed = true;
    }

    Rprintf("\n");

    /* Report progress ~10 times; guard against a zero interval (and the
       resulting division by zero) when the chain is shorter than 10. */
    int report_every = (n_burnin + n_draws) / 10;
    if (report_every < 1)
        report_every = 1;

    for (int s = 0; s < n_burnin + n_draws; s++)
    {
        for (int m = 0; m < n_subgroups; m++)
        {
             sample_gamma_indicators(m, n_platforms, model_platforms_c[m], n_model_platforms_c[m], G, sample_size_ptr[m],
                        ylatent[m], newCC[m], X1[m], gamma, &log_likelihood[m], &logdet[m], &scal[m], nu, theta,
                        n_platform_models_c, platform_models_c, accept_gamma, r, likelihood_type, h[m], h1, h0, hg, K, alpha, psi, sampler_method, &numerical);
            if ((outcome_type == IMR_OUTCOME_SURVIVAL) && (n_censored[m] > 0))
            {
                sample_censored_latent_response(m, n_platforms, model_platforms_c[m], n_model_platforms_c[m], G, sample_size_ptr[m],
                             ylatent[m], yobs[m], newCC[m], X1[m], gamma, &scal[m], &log_likelihood[m],
                             n_censored[m], censored_index[m], logdet[m], r, n_platform_models_c, platform_models_c,
                             accept_y[m], h[m], h1, h0, hg, K, alpha, psi, &numerical);
                if (s >= n_burnin)
                {
                    for (int i = 0; i < n_censored[m]; i++)
                    {
                        int jj = censored_index[m][i];
                        ymean[m][jj] += ylatent[m][jj] / n_draws;
                    }
                }
            }

            if (outcome_type == IMR_OUTCOME_BINARY)
            {
                 sample_binary_latent_response(m, n_platforms, model_platforms_c[m], n_model_platforms_c[m], G, sample_size_ptr[m],
                                  ylatent[m], yobsb[m], newCC[m], X1[m], gamma, &log_likelihood[m],
                                r, n_platform_models_c, platform_models_c,
                               accept_y[m], h[m], h1, h0, hg, K, alpha, psi, &numerical);
                if (s >= n_burnin)
                {
                    for (int i = 0; i < sample_size_ptr[m]; i++)
                    { 
                        // printf(" yyy= %lf",ylatent[m][i]);
                        ymean[m][i] += ylatent[m][i] / n_draws;
                    }   
                }
            }
        } // end of loop with m

        if (strcmp(model_method, "BMS") != 0)
        {
            for (int l = 0; l < n_platforms; l++)
            {
                 sample_mrf_theta(G[l], n_platform_models_c[l], theta[l], accept_theta[l], &mrf[l],
                          gamma[l], nu[l], alpha0, betaTh[l], r);
            }
            if (s >= n_burnin)
            {
                for (int l = 0; l < n_platforms; l++)
                {
                    int m1 = 0;
                    for (int i = 1; i < n_platform_models_c[l]; i++)
                    {
                        for (int j = 0; j < i; j++)
                        {
                            theta_sample[l][s - n_burnin][m1] = theta[l][i][j];
                            m1++;
                        }
                    }
                }
            }
        }
        if (s >= n_burnin)
        {
            for (int l = 0; l < n_platforms; l++)
            {
                for (int m = 0; m < n_platform_models_c[l]; m++)
                {
                    for (int j = 0; j < G[l]; j++)
                    {
                        gamma_mean[l][m][j] += gamma[l][m][j] / (double)n_draws;
                        gamma_sample[s - n_burnin][l][m][j] = gamma[l][m][j];
                    }
                }
            }
        }

        log_posterior_sample[s] = log_posterior(log_likelihood, gamma, nu, theta, mrf, alpha0, betaTh, n_subgroups,
                              n_platforms, G, n_platform_models_c, sampler_method);

        // Print status every 10% of the MCMC samples
        if (s % report_every == 1)
        {
            Rprintf("\nNbr of MCMC samples = %d\n", s);
            Rprintf("LogPosterior=%f\n", log_posterior_sample[s]);
            for (int l = 0; l < n_platforms; l++)
            {
                Rprintf("Nbr of selected features in platform %d = ", l + 1);
                for (int m = 0; m < n_platform_models_c[l]; m++)
                {
                    double sum = 0;
                    for (int g = 0; g < G[l]; g++)
                        sum += gamma[l][m][g];
                    Rprintf("%.4f ", sum);
                }
                Rprintf("\n");
            }
            for (int l = 0; l < n_platforms; l++)
            {
                Rprintf("theta in platform %d = ", l + 1);
                for (int m = 0; m < n_platform_models_c[l]; m++)
                {
                    for (int m1 = 0; m1 < m; m1++)
                    {
                        Rprintf("%f ", theta[l][m][m1]);
                    }
                }
                Rprintf("\n");
            }
        }
    } // end of loop with s index for the number of mcmc samples
    Rprintf("\nAcceptance ratio\n");
    for (int l = 0; l < n_platforms; l++)
    {
        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            Rprintf("%.4f ", accept_gamma[l][m] / (n_draws + n_burnin));
        }
        Rprintf("\n");
    }

    if (accept_y != NULL)
    {
        Rprintf("\nAcceptance ratio for latent Y\n");
        for (int m = 0; m < n_subgroups; m++)
        {
            for (int i = 0; i < sample_size_ptr[m]; i++)
                Rprintf("%.4f ", accept_y[m][i] / (n_draws + n_burnin));
            Rprintf("\n\n");
        }
    }
    int i0, j0;

    // Preserve the historical nested integer-matrix layout and every draw.
    SEXP gamma_sample_R = PROTECT(export_selection_history(gamma_sample,
        n_draws, n_platforms, n_platform_models_c, G));
    protect_count++;

    // export gamma_mean
    SEXP GamMean_R;
    PROTECT(GamMean_R = allocVector(VECSXP, n_platforms));
    protect_count++;
    /// return gamma_mean
    for (int l = 0; l < n_platforms; l++)
    {
        SEXP gamXmeanMatrix = PROTECT(c_array_to_r_matrix(gamma_mean[l], n_platform_models_c[l], G[l]));
        SET_VECTOR_ELT(GamMean_R, l, gamXmeanMatrix);
        UNPROTECT(1);
    }

    if (strcmp(model_method, "BMS") != 0)
    {
        for (int l = 0; l < n_platforms; l++)
        {
            int n_theta_pairs = n_platform_models_c[l] *
                                (n_platform_models_c[l] - 1) / 2;
            if (n_theta_pairs == 0)
                continue;

            double *ThetaXSM = malloc((size_t)n_theta_pairs * sizeof(double));
            if (!ThetaXSM)
                nrerror("allocation failure for theta posterior means");

            mean_array_columns(n_draws, n_theta_pairs,
                               theta_sample[l], ThetaXSM);
            int m1 = 0;
            for (int i = 1; i < n_platform_models_c[l]; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    theta[l][i][j] = theta[l][j][i] = ThetaXSM[m1];
                    m1++;
                }
            }
            free(ThetaXSM);
        }
    }

    // return ThetaMean
    SEXP thetaMean_R;
    PROTECT(thetaMean_R = allocVector(VECSXP, n_platforms));
    protect_count++;
    for (int l = 0; l < n_platforms; l++)
    {
        SEXP thetaMatrix = PROTECT(c_array_to_r_matrix(theta[l], n_platform_models_c[l], n_platform_models_c[l]));
        SET_VECTOR_ELT(thetaMean_R, l, thetaMatrix);
        UNPROTECT(1);
    }

    // return theta_sample
    // --- Convert theta_sample to an R matrix ---
    // For the BMS method theta_sample has already been freed (thetaFreed); in
    // that case we return an empty list element per platform instead of
    // dereferencing freed memory.
    SEXP thetaSampleMatrix_R;
    PROTECT(thetaSampleMatrix_R = allocVector(VECSXP, n_platforms));
    protect_count++;
    for (int l = 0; (l < n_platforms) && (!thetaFreed); l++)
    {
        SEXP thetaSampleMatrix;
        int nrowThetaSample = n_draws;                                             // Rows for ThetaXSample
        int ncolThetaSample = n_platform_models_c[l] * (n_platform_models_c[l] - 1) / 2; // Columns for ThetaXSample
        PROTECT(thetaSampleMatrix = allocVector(REALSXP, nrowThetaSample * ncolThetaSample));
        double *thetaSamplePtr = REAL(thetaSampleMatrix);
        for (i0 = 0; i0 < nrowThetaSample; i0++)
        {
            for (j0 = 0; j0 < ncolThetaSample; j0++)
            {
                thetaSamplePtr[i0 + nrowThetaSample * j0] = theta_sample[l][i0][j0];
            }
        }
        SEXP dimThetaSample;
        PROTECT(dimThetaSample = allocVector(INTSXP, 2));
        INTEGER(dimThetaSample)
        [0] = nrowThetaSample; // Rows
        INTEGER(dimThetaSample)
        [1] = ncolThetaSample; // Columns
        setAttrib(thetaSampleMatrix, R_DimSymbol, dimThetaSample);

        SET_VECTOR_ELT(thetaSampleMatrix_R, l, thetaSampleMatrix);
        UNPROTECT(2);
    }
    // We convert ymean to an R subject
    SEXP YMean_R = PROTECT(array_to_r_list(ymean, n_subgroups, sample_size_ptr));
    protect_count++;

    SEXP logposterior_R;
    PROTECT(logposterior_R = allocVector(REALSXP, n_burnin + n_draws));
    protect_count++;
    for (int s = 0; s < n_burnin + n_draws; s++)
        REAL(logposterior_R)
    [s] = log_posterior_sample[s];

    SEXP rng_state_R = PROTECT(allocVector(RAWSXP, gsl_rng_size(r)));
    protect_count++;
    memcpy(RAW(rng_state_R), gsl_rng_state(r), gsl_rng_size(r));
    int listSize = 7;
    SEXP list;
    SEXP listNames;
    PROTECT(list = allocVector(VECSXP, listSize));
    protect_count++;

    // Add common elements
    SET_VECTOR_ELT(list, 0, GamMean_R);
    SET_VECTOR_ELT(list, 1, thetaMean_R);
    SET_VECTOR_ELT(list, 2, YMean_R);
    SET_VECTOR_ELT(list, 3, logposterior_R);
    SET_VECTOR_ELT(list, 4, gamma_sample_R);
    SET_VECTOR_ELT(list, 5, thetaSampleMatrix_R);
    SET_VECTOR_ELT(list, 6, rng_state_R);

    PROTECT(listNames = allocVector(STRSXP, listSize));
    protect_count++;
    SET_STRING_ELT(listNames, 0, mkChar("gam_mean"));
    SET_STRING_ELT(listNames, 1, mkChar("theta_mean"));
    SET_STRING_ELT(listNames, 2, mkChar("estimate_latent_y"));
    SET_STRING_ELT(listNames, 3, mkChar("log_posterior"));
    SET_STRING_ELT(listNames, 4, mkChar("gam_sample"));
    SET_STRING_ELT(listNames, 5, mkChar("theta_sample"));
    SET_STRING_ELT(listNames, 6, mkChar("rng_state"));
    setAttrib(list, R_NamesSymbol, listNames);

    /// We free memories ...
    for (int s = 0; s < n_draws; s++)
    {
        for (int l = 0; l < n_platforms; l++)
        {
            free_bmatrix(gamma_sample[s][l], 0, n_platform_models_c[l] - 1, 0, G[l] - 1);
        }
        free(gamma_sample[s]);
    }
    free(gamma_sample);

    for (int m = 0; m < n_subgroups; m++)
    {
        free(censored_index[m]);
        free(ylatent[m]);
       // free(yobs[m]);
    }
    free(censored_index);
    free(ylatent);
    if (outcome_type == IMR_OUTCOME_BINARY)
    {
        for (int m = 0; m < n_subgroups; m++)
        {
            free(yobsb[m]);
        }
        free(yobsb);
    } else {
        for (int m = 0; m < n_subgroups; m++)
        {
            free(yobs[m]);
        }
        free(yobs);
    }
    
    if ((outcome_type == IMR_OUTCOME_SURVIVAL) ||
        (outcome_type == IMR_OUTCOME_BINARY))
    {
        for (int m = 0; m < n_subgroups; m++)
        {
            free(accept_y[m]);
        }
        free(accept_y);
    }

    if (!thetaFreed)
    {
        for (int l = 0; l < n_platforms; l++)
        {
            int n_theta_pairs = n_platform_models_c[l] *
                                (n_platform_models_c[l] - 1) / 2;
            if (n_theta_pairs > 0)
                free_dmatrix(theta_sample[l], 0, n_draws - 1,
                             0, n_theta_pairs - 1);
        }
        free(theta_sample);
    }
    gsl_rng_free(r);
    free(mrf);

    for (int l = 0; l < n_platforms; l++)
    {

        for (int m = 0; m < n_platform_models_c[l]; m++)
        {
            free(gamma[l][m]);
            free(gamma_mean[l][m]);
            free(accept_theta[l][m]);
            free(theta[l][m]);
            free(betaTh[l][m]);
        }
        free(gamma[l]);
        free(gamma_mean[l]);
        free(accept_theta[l]);
        free(accept_gamma[l]);
        free(theta[l]);
        free(betaTh[l]);
    }
    free(gamma);
    free(accept_gamma);
    free(gamma_mean);
    free(accept_theta);
    free(theta);
    free(betaTh);
    free_r_list_list_matrix_to_c(X0, n_subgroups, X1_filtered);
    X0 = NULL;

    for (int m = 0; m < n_subgroups; m++)
    {
        free(ymean[m]);
        free_dmatrix(newCC_ptrs[m], 0, sample_size_ptr[m] - 1, 0, K - 1);
        free_dmatrix(newYY_arr[m], 0, sample_size_ptr[m] - 1, 0, 2 - 1);
    }
    free(newCC_ptrs);
    free(newYY_arr);
    free(ymean);

    free(model_platforms_c);
    free(n_model_platforms_c);
    free(platform_models_c);
    free(n_platform_models_c);
    free(h);
    free(log_posterior_sample);

    PutRNGstate();
    Rprintf("\n");
    clock_t t1 = clock() - t;
    double time_taken = ((double)t1) / CLOCKS_PER_SEC; // in seconds
    Rprintf("\nTime taken in minutes before prediction is %f\n", time_taken / 60);
    Rprintf("\nTime taken in hours before prediction is  %f\n", time_taken / 3600);
    /* Keep the return value and all of its components protected across
     * PutRNGstate() and Rprintf(), both of which may allocate and trigger GC. */
    UNPROTECT(protect_count);
    return list;
}
