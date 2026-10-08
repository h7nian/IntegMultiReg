#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include "pmom_proposal.h"

/* Accepted standardized draw for the scalar pMOM kernel. R retains input
 * validation, normalization and the final location/scale transformation. */
SEXP imr_pmom_standardized_draw(SEXP location_ratio_R, SEXP scale_ratio_R,
                                SEXP mixture_probability_R)
{
    if (!isReal(location_ratio_R) || XLENGTH(location_ratio_R) != 1 ||
        !isReal(scale_ratio_R) || XLENGTH(scale_ratio_R) != 1 ||
        !isReal(mixture_probability_R) || XLENGTH(mixture_probability_R) != 1)
        Rf_error("Invalid normalized pMOM inputs");
    double location_ratio = REAL(location_ratio_R)[0];
    double scale_ratio = REAL(scale_ratio_R)[0];
    double mixture_probability = REAL(mixture_probability_R)[0];
    if (!R_FINITE(location_ratio) || !R_FINITE(scale_ratio) ||
        !R_FINITE(mixture_probability) || fabs(location_ratio) > 1 ||
        scale_ratio < 0 || scale_ratio > 1 ||
        fmax(fabs(location_ratio), scale_ratio) != 1 ||
        mixture_probability < 0 || mixture_probability > 1)
        Rf_error("Invalid normalized pMOM inputs");

    GetRNGstate();
    double proposal;
    for (;;) {
        /* Use R's distribution functions, including runif's endpoint rule. */
        proposal = imr_pmom_proposal(mixture_probability);
        double acceptance_uniform = runif(0, 1);
        /* Match R's separate product/addition rounding, not a fused multiply-add.
         * Keep the original powers and arithmetic order in the acceptance ratio. */
        volatile double scaled_proposal = scale_ratio * proposal;
        double numerator = R_pow(location_ratio + scaled_proposal, 2);
        volatile double squared_scale_proposal =
            R_pow(scale_ratio, 2) * R_pow(proposal, 2);
        double denominator = 2 * (R_pow(location_ratio, 2) + squared_scale_proposal);
        double acceptance_probability = numerator / denominator;
        if (ISNAN(acceptance_probability)) {
            PutRNGstate();
            Rf_error("missing value where TRUE/FALSE needed");
        }
        if (acceptance_uniform < acceptance_probability) break;
    }
    PutRNGstate();
    return ScalarReal(proposal);
}
