#ifndef INTEG_MULTI_REG_UTILS_H
#define INTEG_MULTI_REG_UTILS_H

#include <stdbool.h>

#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_rng.h>

/* R/C conversion helpers. R matrices are copied from column-major storage into
 * row-addressable C arrays because the sampler indexes observations first. */
SEXP array_to_r_list(double **array, int rows, int *cols);
SEXP c_array_to_r_matrix(double **array, int rows, int cols);
SEXP c_array_to_r_matrix_int(_Bool **array, int rows, int cols);

double **r_list_vector_double_to_c(int list_length, SEXP list_vector);
double ***r_list_matrix_to_c(int list_length, SEXP list_matrix);
double ****r_list_list_matrix_to_c(int list_length, SEXP list_list_matrix);
_Bool ****r_list_list_matrix_to_c_bool(int list_length, SEXP list_list_matrix);
void free_r_list_list_matrix_to_c(
    double ****array, int list_length, SEXP list_list_matrix);

/* Regression and distribution helpers. */
void ridge_predict_only(const double *x, const double *y, int n, int p,
                        double lambda, double *y_hat);
double r_lefttruncnorm(double lower, double mean, double sd);
double r_righttruncnorm(double upper, double mean, double sd);

/* Accuracy metrics and small vector summaries. */
double norm(int n, double *x);
double sum(int n, double *x);
double mean(int n, double *x);
double var(int n, double *x);

/* Matrix and array helpers. */
void mean_array_columns(int n, int n_cols, double **x, double *mean_out);
_Bool bool_vectors_equal(int n, _Bool *left, _Bool *right);

/* Numerical Recipes style allocators used by the original sampler code. */
double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
_Bool **bmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);

#endif
