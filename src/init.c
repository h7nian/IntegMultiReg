#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_errno.h>

extern SEXP imr_joint_sample(SEXP, SEXP, SEXP, SEXP);
extern SEXP imr_marginal_sample(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imr_collapsed_sample(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imr_concordance(SEXP, SEXP, SEXP);
extern SEXP imr_pmom_standardized_draw(SEXP, SEXP, SEXP);
extern SEXP imr_predict_selection(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imr_joint_sample", (DL_FUNC) &imr_joint_sample, 4},
    {"imr_marginal_sample", (DL_FUNC) &imr_marginal_sample, 5},
    {"imr_collapsed_sample", (DL_FUNC) &imr_collapsed_sample, 25},
    {"imr_concordance", (DL_FUNC) &imr_concordance, 3},
    {"imr_pmom_standardized_draw", (DL_FUNC) &imr_pmom_standardized_draw, 3},
    {"imr_predict_selection", (DL_FUNC) &imr_predict_selection, 28},
    {NULL, NULL, 0}
};

void R_init_IntegMultiReg(DllInfo *dll)
{
    /* Numerical failures in the Laplace backend must not abort the R session. */
    gsl_set_error_handler_off();
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
