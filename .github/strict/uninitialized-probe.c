/* Deliberately invalid R-heap read for an audit control; never packaged. */
#include <R.h>
#include <Rinternals.h>

SEXP imr_uninitialized_probe(void)
{
    SEXP value = PROTECT(allocVector(REALSXP, 1));
    volatile double observed = REAL(value)[0];
    if (observed > 0.0)
        Rprintf("Uninitialized value affected control flow.\n");
    UNPROTECT(1);
    return R_NilValue;
}
