#include "cv_prediction_rows.h"
#include <stdint.h>
#include <string.h>

int postfit_rows_init(postfit_prediction_rows *store, size_t n_rows,
                      size_t n_columns, const int *representatives,
                      void *(*allocate)(size_t), void (*release)(void *))
{
  memset(store, 0, sizeof(*store));
  if (!allocate || !release || n_rows == 0 ||
      n_rows > SIZE_MAX / sizeof(double *) || n_columns > SIZE_MAX / sizeof(double))
    return 0;
  store->n_rows = n_rows;
  store->row_bytes = (n_columns > 0 ? n_columns : 1) * sizeof(double);
  store->representatives = representatives;
  store->allocate = allocate;
  store->release = release;
  store->rows = allocate(n_rows * sizeof(double *));
  if (!store->rows) return 0;
  memset(store->rows, 0, n_rows * sizeof(double *));
  return 1;
}

double *postfit_rows_allocate(postfit_prediction_rows *store, size_t row)
{
  if (!store->rows || row >= store->n_rows || store->rows[row]) return NULL;
  store->rows[row] = store->allocate(store->row_bytes);
  return store->rows[row];
}

int postfit_rows_alias(postfit_prediction_rows *store, size_t row)
{
  if (!store->rows || row >= store->n_rows || store->rows[row] ||
      !store->representatives) return 0;
  int representative = store->representatives[row];
  if (representative < 0 || (size_t)representative >= row ||
      !store->rows[representative]) return 0;
  store->rows[row] = store->rows[representative];
  return 1;
}

void postfit_rows_destroy(postfit_prediction_rows *store)
{
  if (!store->rows) return;
  /* Representatives precede aliases, so reverse traversal never examines a
   * representative after its owned allocation has been released. */
  for (size_t count = store->n_rows; count > 0; --count) {
    size_t row = count - 1;
    if (!store->rows[row]) continue;
    int representative = store->representatives ? store->representatives[row] : -1;
    if (representative >= 0 && (size_t)representative < row &&
        store->rows[row] == store->rows[representative]) continue;
    store->release(store->rows[row]);
  }
  store->release(store->rows);
  store->rows = NULL;
}
