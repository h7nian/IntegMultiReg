/* Test-only adapter: conditional targets are enumerated outside the sampler. */
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include "my_header.h"
#include "utils.h"

/* Two groups, one fixed excluded feature: numerical quadrature in R provides
 * an independent reference for the conditional theta posterior. */
SEXP review_theta(SEXP start_R, SEXP draws_R) {
  int draws = asInteger(draws_R), burnin = 10000;
  double start = asReal(start_R), normalizer;
  double row0[2] = {0, start}, row1[2] = {start, 0};
  double *theta[2] = {row0, row1};
  double a0[2] = {0, 0}, a1[2] = {0, 0}, *accept[2] = {a0, a1};
  double r0[2] = {0, 10}, r1[2] = {10, 0}, *rate[2] = {r0, r1};
  _Bool g0[1] = {0}, g1[1] = {0}, *gamma[2] = {g0, g1};
  SEXP result = PROTECT(allocVector(REALSXP, draws));
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_rand48);
  gsl_rng_set(rng, 913);
  compute_mrf_log_normalizer(2, theta, -3, &normalizer);
  for (int i = 0; i < burnin + draws; ++i) {
    sample_mrf_theta(1, 2, theta, accept, &normalizer, gamma, -3, 40, rate, rng);
    if (i >= burnin) REAL(result)[i - burnin] = theta[0][1];
  }
  gsl_rng_free(rng);
  UNPROTECT(1);
  return result;
}

SEXP review_loglik(SEXP x, SEXP y, SEXP h0R, SEXP hR, SEXP alphaR, SEXP psiR) {
  int n = INTEGER(getAttrib(x, R_DimSymbol))[0];
  int p = INTEGER(getAttrib(x, R_DimSymbol))[1];
  int k = p + 1;
  double h0 = asReal(h0R), h = asReal(hR);
  double **design = malloc(n * sizeof(double *));
  for (int i = 0; i < n; i++) {
    design[i] = malloc(p * sizeof(double));
    for (int j = 0; j < p; j++) design[i][j] = REAL(x)[i + n*j];
  }
  imr_numerical_control numerical = {0, 25, 40, 25, 40, .001};
  double *precision = build_posterior_precision(k, 0, p, n, h, h0, h0, h, design, &numerical);
  double *copy = malloc(k*k*sizeof(double));
  memcpy(copy, precision, k*k*sizeof(double));
  gsl_set_error_handler_off();
  gsl_matrix_view m = gsl_matrix_view_array(precision, k, k);
  gsl_linalg_cholesky_decomp(&m.matrix);
  double *beta = malloc(k*sizeof(double));
  double ll = log_likelihood_nonlocal(k, 0, p, n, asReal(alphaR), asReal(psiR), REAL(y), design,
    copy, &m.matrix, beta, 1, h, h0, h0, h, 40, 1e-3, 0, &numerical, IMR_LAPLACE_PREDICTION);
  for (int i = 0; i < n; i++) free(design[i]);
  free(design); free(precision); free(copy); free(beta);
  return ScalarReal(ll);
}

SEXP review_prior(SEXP thetaR, SEXP stateR, SEXP methodR) {
  double t = asReal(thetaR), nu[1] = {-3}, mrf[1], ll[2] = {0,0};
  int G[1] = {1}, counts[1] = {2};
  double t0[2] = {0,t}, t1[2] = {t,0}, *tm[2] = {t0,t1}, **theta[1] = {tm};
  double r0[2] = {10,10}, r1[2] = {10,10}, *rm[2] = {r0,r1}, **rate[1] = {rm};
  _Bool g0[1] = {INTEGER(stateR)[0]}, g1[1] = {INTEGER(stateR)[1]}, *gm[2] = {g0,g1}, **gamma[1] = {gm};
  compute_mrf_log_normalizer(2, tm, nu[0], mrf);
  return ScalarReal(log_posterior(ll, gamma, nu, theta, mrf, 40, rate, 2, 1, G, counts, asInteger(methodR)));
}

SEXP review_precision(SEXP historical_R) {
  double row[3] = {0, 0, 0}, *design[1] = {row};
  imr_numerical_control numerical = {asInteger(historical_R), 25, 40, 25, 40, .001};
  double *matrix = build_posterior_precision(4, 1, 1, 1, 2, 3, 5, 7, design, &numerical);
  SEXP result = PROTECT(allocVector(REALSXP, 4));
  for (int i = 0; i < 4; ++i) REAL(result)[i] = matrix[4*i+i];
  free(matrix);
  UNPROTECT(1);
  return result;
}
