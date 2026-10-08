# Run from the package root. Writes fresh independent references outside fixtures.
# For observed-data quadrature allow several minutes and several GiB of memory.
args <- commandArgs(TRUE)
out <- normalizePath(if (length(args)) args[1L] else tempdir(), mustWork = TRUE)
setwd('.github/reference')
source('exact-reference.R')
source('quadrature-reference.R')
source('availability-reference.R')
for (case in c('continuous', 'correlated-swap')) {
  groups <- readRDS(file.path('fixtures', paste0(case, '-data.rds')))
  reference <- exact_two_group_posterior(groups)
  saveRDS(reference, file.path(out, paste0(case, '-reference.rds')))
}
for (case in c('binary', 'right.censored')) {
  groups <- readRDS(file.path('fixtures', paste0(case, '-data.rds')))
  coarse <- observed_two_group_reference(groups, normal_order = 35, gamma_order = 60)
  fine <- observed_two_group_reference(groups, normal_order = 55, gamma_order = 100)
  vectorize <- function(r)c(r$probability, unlist(r$beta_mean), unlist(r$beta_second),
    r$variance_mean, r$theta_mean, unlist(r$latent_mean))
  max_difference <- max(abs(vectorize(coarse) - vectorize(fine)))
  stopifnot(max_difference < .001)
  saveRDS(list(coarse = coarse, fine = fine, max_difference = max_difference),
    file.path(out, paste0(case, '-quadrature.rds')))
}
groups <- readRDS('fixtures/availability-data.rds')
for (variant in c('imr', 'bms')) {
  reference <- gaussian_availability_reference(groups, c(1L, 2L), c(-1, -1.5), model_variant = variant)
  saveRDS(reference, file.path(out, paste0('availability-', variant, '-reference.rds')))
}
for (name in list.files(out, pattern = '-(reference|quadrature)[.]rds$')) {
  if (file.exists(file.path('fixtures', name)))
    stopifnot(isTRUE(all.equal(readRDS(file.path(out, name)), readRDS(file.path('fixtures', name)), tolerance = 1e-9)))
}
cat('Independent regeneration matches stored references within 1e-9.\n')
