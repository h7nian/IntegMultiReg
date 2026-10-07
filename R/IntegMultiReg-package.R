#' Joint Bayesian Regression for Multi-Platform Biomarkers
#'
#' IntegMultiReg fits availability-subgroup regressions with shared molecular
#' selection through a Markov random field and non-local pMOM coefficient priors.
#' Continuous, binary and right-censored outcomes use one joint posterior
#' sampling engine. Clinical covariates are included in every subgroup model.
#'
#' Start with [imr_data()] and [imr()]. [imr_priors()] defines the statistical
#' model and [imr_mcmc()] defines the sampling budget. Standard [coef()],
#' [confint()], [summary()] and [predict()] methods use the stored joint draws.
#' [inclusion_probabilities()] reports feature selection; [posterior_draws()]
#' extracts samples; [mcmc_diagnostics()] reports their exploration and precision.
#' [cv_imr()] refits the complete workflow or performs diagnosed posterior
#' reweighting. [compare_fit_summaries()] gives descriptive comparisons only.
#'
#' @useDynLib IntegMultiReg, .registration = TRUE
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis box image par
#' @importFrom stats pnorm sd
#' @importFrom utils capture.output
#' @keywords internal
"_PACKAGE"
