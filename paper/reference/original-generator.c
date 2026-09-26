/* Test/research adapter for archived ReadData.c; never linked into the package. */
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "myheader.h"
#include "utils.h"
#include <stdlib.h>
int K = 0;
static Matrix copy_matrix(SEXP value, int columns) {
    SEXP dim = getAttrib(value, R_DimSymbol);
    if (!isReal(value) || LENGTH(dim) != 2 || INTEGER(dim)[1] != columns)
        Rf_error("Original generator requires the archived feature dimensions");
    int n = INTEGER(dim)[0];
    Matrix matrix = {.nR = n, .nC = columns, .Mat = dmatrix(0, n-1, 0, columns-1)};
    for (int i = 0; i < n; ++i) for (int j = 0; j < columns; ++j)
        matrix.Mat[i][j] = REAL(value)[i+n*j];
    Normalize(n, columns, matrix.Mat);
    return matrix;
}
SEXP original_generate(SEXP groups, SEXP scenario_R, SEXP seed_R, SEXP rho_R, SEXP half_R) {
    Matrix X[4], Z[2], U[2]; SurvTime T[4];
    int selection_x[4][10], selection_z[2][6], selection_u[2][10];
    _Bool truth_x[4][G], truth_z[2][M], truth_u[2][D];
    for (int g = 0; g < 4; ++g) {
        X[g] = copy_matrix(VECTOR_ELT(VECTOR_ELT(groups, 0), g), G);
        T[g].n = X[g].nR;
        T[g].Time = (double *)R_alloc(T[g].n, sizeof(double));
        T[g].delta = (_Bool *)R_alloc(T[g].n, sizeof(_Bool));
    }
    for (int g = 0; g < 2; ++g) {
        Z[g] = copy_matrix(VECTOR_ELT(VECTOR_ELT(groups, 1), g), M);
        U[g] = copy_matrix(VECTOR_ELT(VECTOR_ELT(groups, 2), g), D);
    }
    SelectedBiomark(asInteger(scenario_R), 10, 6, selection_x, selection_z, selection_u);
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_rand48);
    gsl_rng_set(rng, (unsigned long)asReal(seed_R));
    generdata(asInteger(half_R), 10, 6, selection_x, selection_z, selection_u,
        T, X, Z, U, -1, 1, truth_x, truth_z, truth_u, asReal(rho_R), rng);
    gsl_rng_free(rng);
    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP response = PROTECT(allocVector(VECSXP, 4));
    SEXP features = PROTECT(allocVector(VECSXP, 3));
    SEXP truth = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result, 0, response); SET_VECTOR_ELT(result, 1, features); SET_VECTOR_ELT(result, 2, truth);
    for (int g = 0; g < 4; ++g) {
        SEXP y = PROTECT(allocMatrix(REALSXP, T[g].n, 2));
        for (int i = 0; i < T[g].n; ++i) {
            REAL(y)[i] = T[g].Time[i]; REAL(y)[i+T[g].n] = T[g].delta[i];
        }
        SET_VECTOR_ELT(response, g, y); UNPROTECT(1);
    }
    Matrix *platform[3] = {X, Z, U}; int counts[3] = {4, 2, 2}, columns[3] = {G, M, D};
    for (int p = 0; p < 3; ++p) {
        SEXP matrices = PROTECT(allocVector(VECSXP, counts[p]));
        SEXP selected = PROTECT(allocMatrix(INTSXP, counts[p], columns[p]));
        for (int g = 0; g < counts[p]; ++g) {
            Matrix matrix = platform[p][g];
            SEXP value = PROTECT(allocMatrix(REALSXP, matrix.nR, matrix.nC));
            for (int j = 0; j < matrix.nC; ++j) {
                INTEGER(selected)[g+counts[p]*j] = p == 0 ? truth_x[g][j] : p == 1 ? truth_z[g][j] : truth_u[g][j];
                for (int i = 0; i < matrix.nR; ++i) REAL(value)[i+matrix.nR*j] = matrix.Mat[i][j];
            }
            SET_VECTOR_ELT(matrices, g, value); UNPROTECT(1);
            free_dmatrix(matrix.Mat, 0, matrix.nR-1, 0, matrix.nC-1);
        }
        SET_VECTOR_ELT(features, p, matrices); SET_VECTOR_ELT(truth, p, selected); UNPROTECT(2);
    }
    UNPROTECT(4); return result;
}
