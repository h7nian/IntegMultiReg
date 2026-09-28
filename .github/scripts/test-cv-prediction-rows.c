#include "../../src/cv_prediction_rows.h"
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

/* Fault injection and allocation accounting exist only in this test binary. */
static void *live[32];
static size_t n_live, calls, fail_at;
static void *allocate_test(size_t bytes) {
  if (++calls == fail_at) return NULL;
  void *pointer = malloc(bytes);
  assert(pointer && n_live < 32);
  live[n_live++] = pointer;
  return pointer;
}
static void release_test(void *pointer) {
  size_t index = 0;
  while (index < n_live && live[index] != pointer) ++index;
  assert(index < n_live); /* Detect double frees and unowned pointers. */
  live[index] = live[--n_live];
  free(pointer);
}
static void verify(size_t columns, size_t failure, int use_aliases) {
  const int representatives[] = {0, 0, 2, 0, 4, 4, 6, 0};
  postfit_prediction_rows store;
  calls = n_live = 0;
  fail_at = failure;
  if (postfit_rows_init(&store, 8, columns, representatives, allocate_test, release_test)) {
    for (size_t row = 0; row < 8; ++row) {
      /* Row 5 deliberately owns a distinct allocation despite its representative,
       * exercising a cache miss after an unusable representative log-weight. */
      if (use_aliases && row != 5 && representatives[row] < (int)row) {
        assert(postfit_rows_alias(&store, row));
        assert(store.rows[row] == store.rows[representatives[row]]);
      } else {
        double *values = postfit_rows_allocate(&store, row);
        if (!values) break;
        for (size_t column = 0; column < columns; ++column) values[column] = row + column;
        assert(!postfit_rows_allocate(&store, row));
      }
    }
    assert(!postfit_rows_alias(&store, 8));
    assert(!postfit_rows_allocate(&store, 8));
  }
  postfit_rows_destroy(&store);
  postfit_rows_destroy(&store);
  assert(n_live == 0);
}
int main(void) {
  for (size_t columns = 0; columns <= 3; columns += 3)
    for (size_t failure = 0; failure <= 10; ++failure)
      for (int aliases = 0; aliases <= 1; ++aliases) verify(columns, failure, aliases);
  postfit_prediction_rows store;
  assert(!postfit_rows_init(&store, SIZE_MAX, 1, NULL, malloc, free));
  postfit_rows_destroy(&store);
  assert(!postfit_rows_init(&store, 1, SIZE_MAX, NULL, malloc, free));
  postfit_rows_destroy(&store);
  puts("44 ownership/failure scenarios and size-overflow guards pass.");
  return 0;
}
