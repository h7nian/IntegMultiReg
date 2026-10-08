#' Joint Bayesian Regression for Multi-Platform Biomarkers
#'
#' IntegMultiReg fits availability-subgroup regressions with shared molecular
#' selection through a Markov random field and non-local pMOM coefficient priors.
#' Continuous, binary and right-censored outcomes support four choices of
#' marginalized regression parameters. Clinical covariates are always included.
#'
#' Start with [imr_data()] and [imr()]. [imr_priors()] defines the statistical
#' model, [imr_mcmc()] defines the sampling budget and [imr_control()] defines
#' numerical settings. Standard [coef()],
#' [confint()], [summary()] and [predict()] methods use the stored joint draws.
#' [inclusion_probabilities()] reports feature selection; [posterior_draws()]
#' extracts samples; [mcmc_diagnostics()] reports their exploration and precision.
#' Laplace selection fits use [predict.imr_selection()] and optional
#' [sample_regression_posterior()] for regression uncertainty.
#' [cv_imr()] refits the complete workflow or performs diagnosed posterior
#' reweighting. [compare_fit_summaries()] gives descriptive comparisons only.
#'
#' @useDynLib IntegMultiReg, .registration = TRUE
#' @importFrom stats sd
#' @keywords internal
"_PACKAGE"
