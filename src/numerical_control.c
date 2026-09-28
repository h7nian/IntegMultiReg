#include <R.h>
#include <limits.h>
#include <math.h>
#include "numerical_control.h"

imr_numerical_control imr_read_numerical_control(SEXP control)
{
    if (!isReal(control) || XLENGTH(control) != 6)
        Rf_error("Invalid numerical controls");
    const double *value = REAL(control);
    for (int i = 0; i < 6; ++i)
        if (!R_FINITE(value[i])) Rf_error("Numerical controls must be finite");
    if (value[0] != 0 && value[0] != 1) Rf_error("Invalid prior indexing");
    for (int i = 1; i <= 4; ++i)
        if (value[i] < 1 || value[i] > INT_MAX || floor(value[i]) != value[i])
            Rf_error("Laplace iteration limits must be positive integers");
    if (value[5] <= 0) Rf_error("Laplace tolerance must be positive");
    return (imr_numerical_control){(int)value[0], (int)value[1],
        (int)value[2], (int)value[3], (int)value[4], value[5], NULL};
}

void imr_record_laplace(const imr_numerical_control *control, int stage, int measure)
{
    if (!control || !control->diagnostics || stage == IMR_LAPLACE_PREDICTION) return;
    imr_laplace_diagnostics *d = control->diagnostics;
    d->values[3 * d->subgroup + stage + (size_t)measure * 3 * d->n_subgroups] += 1;
}
