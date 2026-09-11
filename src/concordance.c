#include <stdbool.h>

#include <R.h>
#include <Rinternals.h>

static double concordance_index(int n, const double *prediction,
                                const double *observed_time,
                                const _Bool *event)
{
    double comparable = 0.0;
    double concordant = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            int early = -1;
            int late = -1;
            if (event[i] && (observed_time[i] < observed_time[j] ||
                (observed_time[i] == observed_time[j] && !event[j])))
            {
                early = i;
                late = j;
            }
            else if (event[j] && (observed_time[j] < observed_time[i] ||
                     (observed_time[j] == observed_time[i] && !event[i])))
            {
                early = j;
                late = i;
            }
            if (early < 0)
                continue;
            comparable += 1.0;
            concordant += prediction[early] < prediction[late] ? 1.0 :
                          prediction[early] == prediction[late] ? 0.5 : 0.0;
        }
    }
    return comparable > 0.0 ? concordant / comparable : NA_REAL;
}

SEXP imr_concordance(SEXP prediction, SEXP time, SEXP status)
{
    const int n = LENGTH(prediction);
    if (!isReal(prediction) || !isReal(time) || !isInteger(status) ||
        LENGTH(time) != n || LENGTH(status) != n)
        Rf_error("Invalid concordance inputs");

    _Bool *event = (_Bool *)R_alloc(n, sizeof(_Bool));
    for (int i = 0; i < n; ++i)
    {
        if (!R_FINITE(REAL(prediction)[i]) || !R_FINITE(REAL(time)[i]) ||
            (INTEGER(status)[i] != 0 && INTEGER(status)[i] != 1))
            Rf_error("Non-finite or invalid concordance inputs");
        event[i] = (_Bool)INTEGER(status)[i];
    }
    return ScalarReal(concordance_index(
        n, REAL(prediction), REAL(time), event));
}
