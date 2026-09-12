#include "../../src/cv_precision.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
static double *reference(int n_coefficients, int train_sample_size, const int *train_index, double **design) {
 int j,j1,i1,i2;
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


 return precision;
}
static double *candidate(int n_coefficients, int train_sample_size,
                         const int *train_index, double **design) {
  double *precision = malloc((size_t)n_coefficients * n_coefficients * sizeof(double));
  if (!precision) abort();
  postfit_build_precision(precision, design, n_coefficients, train_sample_size,
                          train_index, IMR_CV_WORKSPACE_BYTES);
  return precision;
}
static double *forced_fallback(int n_coefficients, int train_sample_size,
                               const int *train_index, double **design) {
  double *precision = malloc((size_t)n_coefficients * n_coefficients * sizeof(double));
  if (!precision) abort();
  postfit_build_precision(precision, design, n_coefficients, train_sample_size,
                          train_index, 0);
  return precision;
}
static void verify(int rows, int columns, int pattern) {
  int coefficients = columns + 1;
  int *indices = malloc((size_t)(rows + 1) * sizeof(int));
  double **design = malloc((size_t)(rows + 1) * sizeof(double *));
  double *values = malloc((size_t)(rows + 1) * (columns + 1) * sizeof(double));
  if (!indices || !design || !values) abort();
  for (int row=0; row<rows; ++row) {
    indices[row] = rows - 1 - row;
    design[row] = values + (size_t)row * (columns + 1);
    for (int column=0; column<columns; ++column)
      design[row][column] = pattern == 0 ? 1.0 :
        pattern == 1 ? (row % 7) * .1 + column * 1e-12 :
        ((row * 17 + column * 11) % 101 - 50) * .03;
  }
  double *expected = reference(coefficients,rows,indices,design);
  double *actual = candidate(coefficients,rows,indices,design);
  double *fallback = forced_fallback(coefficients,rows,indices,design);
  for (size_t i=0; i<(size_t)coefficients*coefficients; ++i) {
    if (!isfinite(actual[i]) || fabs(actual[i]-expected[i]) > 1e-10 + 1e-8*fabs(expected[i]) ||
        memcmp(fallback+i,expected+i,sizeof(double)) != 0) {
      fprintf(stderr,"Mismatch rows=%d columns=%d pattern=%d cell=%zu\n",rows,columns,pattern,i);
      exit(1);
    }
  }
  gsl_matrix_view expected_matrix = gsl_matrix_view_array(expected, coefficients, coefficients);
  gsl_matrix_view actual_matrix = gsl_matrix_view_array(actual, coefficients, coefficients);
  int expected_status = gsl_linalg_cholesky_decomp(&expected_matrix.matrix);
  int actual_status = gsl_linalg_cholesky_decomp(&actual_matrix.matrix);
  if (expected_status != actual_status) {
    fprintf(stderr, "Cholesky status changed: rows=%d columns=%d pattern=%d\n", rows, columns, pattern);
    exit(1);
  }
  free(expected);free(actual);free(fallback);free(indices);free(design);free(values);
}
int main(void) {
  gsl_set_error_handler_off();
  const int rows[] = {0,1,2,7,50,200};
  const int columns[] = {0,1,4,29};
  for (int r=0;r<6;++r)
    for (int c=0;c<4;++c)
      for (int p=0;p<3;++p) verify(rows[r],columns[c],p);
  verify(600000,2,2); /* Actual 8 MiB budget fallback, not a changed test budget. */
  puts("73 kernel cases pass, including forced and actual-budget fallbacks.");
  return 0;
}
