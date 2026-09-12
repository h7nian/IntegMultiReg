#include "cv_precision.h"
#include <stdlib.h>
#include <gsl/gsl_blas.h>

void postfit_build_precision(double *precision, double **design,
                             int n_coefficients, int train_sample_size,
                             const int *train_index, size_t workspace_bytes)
{
    int row_coefficient, column_coefficient, source_row, training_row;
    /* Training-only matrix product; no full-data subtraction.
     * Large designs or allocation failure retain the original indirect path. */
    double *train_columns = NULL;
    size_t column_count = (size_t)n_coefficients - 1;
    size_t workspace_limit = workspace_bytes / sizeof(double);
    if (train_sample_size > 0 && column_count > 0 &&
        column_count <= workspace_limit / (size_t)train_sample_size) {
      train_columns = malloc(column_count * (size_t)train_sample_size * sizeof(double));
      if (train_columns != NULL) {
        for (size_t column = 0; column < column_count; ++column)
          for (int row = 0; row < train_sample_size; ++row)
            train_columns[column * (size_t)train_sample_size + row] =
              design[train_index[row]][column];
      }
    }

    if (train_columns != NULL) {
      gsl_matrix_view columns = gsl_matrix_view_array(train_columns, column_count,
                                                       train_sample_size);
      gsl_matrix_view products = gsl_matrix_view_array_with_tda(
        precision + n_coefficients + 1, column_count, column_count, n_coefficients);
      gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &columns.matrix, 0.0,
                     &products.matrix);
    }
    for (row_coefficient = 0; row_coefficient < n_coefficients; row_coefficient++)
    {
      for (column_coefficient = 0; column_coefficient <= row_coefficient; column_coefficient++)
      {
        double crossproduct = 0;
        if ((row_coefficient == 0) && (column_coefficient == 0))
        {
          crossproduct = train_sample_size;
        }
        else if (column_coefficient == 0)
        {
          for (training_row = 0; training_row < train_sample_size; training_row++)
          {
            source_row = train_index[training_row];
            crossproduct += train_columns != NULL ?
              train_columns[(size_t)(row_coefficient - 1) * train_sample_size + training_row] : design[source_row][row_coefficient - 1];
          }
        }
        else if ((row_coefficient != 0) && (column_coefficient != 0))
        {
          if (train_columns != NULL) {
            crossproduct = precision[row_coefficient * n_coefficients + column_coefficient];
          } else {
            for (training_row = 0; training_row < train_sample_size; training_row++) {
              source_row = train_index[training_row];
              crossproduct += design[source_row][row_coefficient - 1] * design[source_row][column_coefficient - 1];
            }
          }
        }
        if (row_coefficient == column_coefficient)
          crossproduct += .001; // Preserve the existing diagonal ridge exactly.
        precision[row_coefficient * n_coefficients + column_coefficient] = precision[column_coefficient * n_coefficients + row_coefficient] = crossproduct;
      }
    }

    free(train_columns);
}
