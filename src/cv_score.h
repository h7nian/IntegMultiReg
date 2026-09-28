#ifndef IMR_CV_SCORE_H
#define IMR_CV_SCORE_H
#include <stdbool.h>
double imr_legacy_concordance(int n, double *prediction, double *observed_time, _Bool *event);
double imr_legacy_auc(int n, double *prediction, _Bool *outcome);
#endif
