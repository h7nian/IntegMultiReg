#include <math.h>
#include "sampler_math.h"

double imr_gamma_log_odds(int n_groups, int group, int feature,
                          const double *interaction, _Bool **selection,
                          double nu, int sampler_method)
{
    double neighbors = 0;
    for (int other = 0; other < n_groups; ++other)
        if (other != group)
            neighbors += interaction[other] * selection[other][feature];
    if (sampler_method == IMR_SAMPLER_PAPER) neighbors *= 2;
    return neighbors + nu;
}

double imr_gamma_log_hastings(int n_features, int current_count,
                              int proposed_count, double flip_probability)
{
    if (current_count == proposed_count) return 0;
    double forward = current_count == 0 || current_count == n_features ? 1 : flip_probability;
    double reverse = proposed_count == 0 || proposed_count == n_features ? 1 : flip_probability;
    return log(reverse / forward);
}
