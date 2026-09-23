#ifndef IMR_CV_PREDICTION_ROWS_H
#define IMR_CV_PREDICTION_ROWS_H

#include <stddef.h>

typedef struct {
  double **rows;
  size_t n_rows;
  size_t row_bytes;
  const int *representatives;
  void *(*allocate)(size_t);
  void (*release)(void *);
} postfit_prediction_rows;

int postfit_rows_init(postfit_prediction_rows *store, size_t n_rows,
                      size_t n_columns, const int *representatives,
                      void *(*allocate)(size_t), void (*release)(void *));
double *postfit_rows_allocate(postfit_prediction_rows *store, size_t row);
int postfit_rows_alias(postfit_prediction_rows *store, size_t row);
void postfit_rows_destroy(postfit_prediction_rows *store);

#endif
