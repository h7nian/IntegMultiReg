#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP imr_joint_sample(SEXP, SEXP, SEXP, SEXP);
extern SEXP imr_concordance(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imr_joint_sample", (DL_FUNC) &imr_joint_sample, 4},
    {"imr_concordance", (DL_FUNC) &imr_concordance, 3},
    {NULL, NULL, 0}
};

void R_init_IntegMultiReg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
