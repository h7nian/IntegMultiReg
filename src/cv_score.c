/* Historical metrics shared by both CV engines. */
#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include "utils.h"
#include "my_header.h"
#include "cv_score.h"
/* Historical 0.1.0 concordance, including its original tie convention. */
double imr_legacy_concordance(int n, double *prediction, double *observed_time, _Bool *event)
{
  int i, j;
  double concordance_denominator = 0;
  double concordance_numerator = 0;
  double time1, time2, prediction1, prediction2;
  for (i = 0; i < n; i++)
  {
    time1 = observed_time[i];
    prediction1 = prediction[i];
    for (j = 0; j < n; j++)
    {
      if (i != j)
      {
        time2 = observed_time[j];
        prediction2 = prediction[j];
        concordance_numerator +=
            (prediction2 > prediction1) * (time2 > time1) * (event[i] == 1) +
            (prediction2 < prediction1) * (time2 < time1) * (event[j] == 1) +
            0.5 * ((prediction2 == prediction1) || (time2 == time1)) * (event[i] == 1) * (event[j] == 0) +
            0.5 * ((prediction2 == prediction1) || (time2 == time1)) * (event[j] == 1) * (event[i] == 0);
        concordance_denominator +=
            (time2 > time1) * (event[i] == 1) +
            (time2 < time1) * (event[j] == 1) +
            (time2 == time1) * (event[i] == 1) * (event[j] == 0) +
            (time2 == time1) * (event[i] == 0) * (event[j] == 1);
      }
    }
  }
  return concordance_denominator > 0 ?
      concordance_numerator / concordance_denominator : NA_REAL;
}


double imr_legacy_auc(int n, double *esti, _Bool * class)
{
    if (n == 0) return NA_REAL;
    int positives = 0;
    for (int i = 0; i < n; ++i) positives += class[i];
    if (positives == 0 || positives == n) return NA_REAL;
    double fpr[n + 2], tpr[n + 2];
    double auc1 = 0;
    int i, j;
    double esti1[n];
    for (i = 0; i < n; i++)
    {
        esti1[i] = esti[i];
    }
    int idx[n];
    sort_descending_index(n, esti1, idx);

    fpr[n + 1] = 1;
    tpr[n + 1] = 1;
    fpr[0] = 0;
    tpr[0] = 0;
    for (i = n; i >= 1; --i)
    {
        double af = 0;
        double at = 0;
        for (j = 0; j < n; j++)
        {
            if (esti[j] > esti1[i - 1])
            {
                if (class[j] == 0)
                {
                    af += 1;
                }
                else
                {
                    at += 1;
                }
            }
        }
        tpr[i] = at / positives;
        fpr[i] = af / (n - positives);
        auc1 += (fpr[i + 1] - fpr[i]) * (tpr[i + 1] + tpr[i]);
    }
    auc1 += (fpr[1] - fpr[0]) * (tpr[1] + tpr[0]);
    auc1 = 0.5 * (auc1);
    return auc1;
}



SEXP imr_cv_legacy_score(SEXP prediction_R, SEXP response_R, SEXP event_R, SEXP type_R)
{
    if (!isReal(prediction_R) || !isReal(response_R) ||
        XLENGTH(prediction_R) != XLENGTH(response_R) || XLENGTH(prediction_R) > INT_MAX ||
        !isInteger(type_R) || XLENGTH(type_R) != 1 ||
        INTEGER(type_R)[0] < 0 || INTEGER(type_R)[0] > 1)
        Rf_error("Invalid historical scoring inputs");
    int n = (int)XLENGTH(prediction_R);
    int survival = INTEGER(type_R)[0] == 0;
    if (survival && (!isInteger(event_R) || XLENGTH(event_R) != n))
        Rf_error("Invalid historical scoring events");
    _Bool *indicator = (_Bool *)R_alloc(n > 0 ? n : 1, sizeof(_Bool));
    for (int i = 0; i < n; ++i) {
        if (!R_FINITE(REAL(prediction_R)[i]) || !R_FINITE(REAL(response_R)[i]))
            Rf_error("Historical scoring inputs must be finite");
        if (survival) {
            int event = INTEGER(event_R)[i];
            if (event != 0 && event != 1) Rf_error("Historical scoring events must be 0 or 1");
            indicator[i] = event;
        } else {
            double response = REAL(response_R)[i];
            if (response != 0 && response != 1) Rf_error("Historical scoring outcomes must be 0 or 1");
            indicator[i] = response;
        }
    }
    return ScalarReal(survival ? imr_legacy_concordance(n, REAL(prediction_R), REAL(response_R), indicator)
                              : imr_legacy_auc(n, REAL(prediction_R), indicator));
}
