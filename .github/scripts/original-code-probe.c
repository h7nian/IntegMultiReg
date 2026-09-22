/* Adapter for the unmodified 2017 functions, compiled only by the reference
 * verification script. G/M/D retain the dimensions in the original header;
 * unused molecular columns are zero and unselected. */
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include "myheader.h"
#include "utils.h"

double h0 = 10000, alpha = .37, psi = .21;
int K = 0;

SEXP original_cv(SEXP design_R, SEXP response_R, SEXP states_R,
                 SEXP train_R, SEXP test_R)
{
    int n = INTEGER(getAttrib(design_R, R_DimSymbol))[0];
    int p = INTEGER(getAttrib(design_R, R_DimSymbol))[1];
    int draws = INTEGER(getAttrib(states_R, R_DimSymbol))[0];
    int n_test = LENGTH(test_R);
    if (p > G || p != INTEGER(getAttrib(states_R, R_DimSymbol))[1])
        Rf_error("Invalid reference design");
    gsl_set_error_handler_off();
    double **x = dmatrix(0, n - 1, 0, G - 1);
    double **dummy = dmatrix(0, n - 1, 0, 0);
    for (int i = 0; i < n; ++i) {
        dummy[i][0] = 0;
        for (int j = 0; j < G; ++j) x[i][j] = j < p ? REAL(design_R)[i + n * j] : 0;
    }
    _Bool ***genes = malloc(draws * sizeof(*genes));
    _Bool ***mirna = malloc(draws * sizeof(*mirna));
    _Bool ***methylation = malloc(draws * sizeof(*methylation));
    int *indices = (int *)R_alloc(draws, sizeof(int));
    for (int draw = 0; draw < draws; ++draw) {
        genes[draw] = bmatrix(0, 3, 0, G - 1);
        mirna[draw] = bmatrix(0, 1, 0, M - 1);
        methylation[draw] = bmatrix(0, 1, 0, D - 1);
        indices[draw] = draw;
        for (int group = 0; group < 4; ++group)
            for (int j = 0; j < G; ++j)
                genes[draw][group][j] = j < p ? INTEGER(states_R)[draw + draws * j] : 0;
        for (int group = 0; group < 2; ++group) {
            for (int j = 0; j < M; ++j) mirna[draw][group][j] = 0;
            for (int j = 0; j < D; ++j) methylation[draw][group][j] = 0;
        }
    }
    double *prediction = pred_aftcv(4, n, n_test, INTEGER(test_R), INTEGER(train_R),
        REAL(response_R), dummy, x, dummy, dummy, genes, mirna, methylation,
        draws, indices, indices);
    SEXP result = PROTECT(allocVector(REALSXP, n_test));
    for (int i = 0; i < n_test; ++i) REAL(result)[i] = prediction[i];
    free(prediction);
    for (int draw = 0; draw < draws; ++draw) {
        for (int group = 0; group < 4; ++group) free(genes[draw][group]);
        for (int group = 0; group < 2; ++group) {
            free(mirna[draw][group]);
            free(methylation[draw][group]);
        }
        free(genes[draw]); free(mirna[draw]); free(methylation[draw]);
    }
    free(genes); free(mirna); free(methylation);
    free_dmatrix(x, 0, n - 1, 0, G - 1);
    free_dmatrix(dummy, 0, n - 1, 0, 0);
    UNPROTECT(1);
    return result;
}
