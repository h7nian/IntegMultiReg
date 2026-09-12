#ifndef IMR_CV_PRECISION_H
#define IMR_CV_PRECISION_H

#include <stddef.h>

#define IMR_CV_WORKSPACE_BYTES ((size_t)8 * 1024 * 1024)

/* Fills caller-owned storage; optional scratch failure uses the scalar path. */
void postfit_build_precision(double *precision, double **design,
                             int n_coefficients, int train_sample_size,
                             const int *train_index, size_t workspace_bytes);

#endif
