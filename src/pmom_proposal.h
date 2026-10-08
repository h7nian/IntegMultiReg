#ifndef IMR_PMOM_PROPOSAL_H
#define IMR_PMOM_PROPOSAL_H

#include <Rmath.h>
#include <math.h>

/* Common mixture proposal. Each caller retains its original acceptance
   arithmetic and location/scale rounding, and owns the R RNG state. */
static double imr_pmom_proposal(double mixture_probability)
{
    if (runif(0, 1) < mixture_probability) return rnorm(0, 1);
    double sign = runif(0, 1) < .5 ? -1 : 1;
    return sign * sqrt(rchisq(3));
}

#endif
