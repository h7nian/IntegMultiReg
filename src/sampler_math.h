#ifndef IMR_SAMPLER_MATH_H
#define IMR_SAMPLER_MATH_H

#include <stdbool.h>

/* The symmetric MRF convention in Appendix D counts each edge twice. */
double imr_gamma_log_odds(int n_groups, int group, int feature,
                          const double *interaction, _Bool **selection,
                          double nu);

/* Reverse/forward probability for the existing boundary-aware flip/swap
 * proposal. Swaps leave the selected count unchanged and are symmetric. */
double imr_gamma_log_hastings(int n_features, int current_count,
                              int proposed_count, double flip_probability);

#endif
