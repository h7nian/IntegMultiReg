## Generates the `simIMR` example data set shipped with IntegMultiReg.
##
## Three platforms with deliberately overlapping but partially missing subjects
## induce several non-empty Venn subgroups; a handful of truly associated
## features per platform drive a latent linear predictor from which continuous,
## binary and right-censored outcomes are derived.  Run from the package root:
##   Rscript data-raw/make_simIMR.R

set.seed(2024)

N <- 300
id <- seq_len(N)

## ---- platform membership (induces the availability subgroups) -------------
## X (genomic)      : ids   1..240
## Z (proteomic)    : ids   1..180
## U (metabolomic)  : ids 121..300
idX <- 1:240
idZ <- 1:180
idU <- 121:300

p1 <- 20L  # genomic features
p2 <- 10L  # proteomic features
p3 <- 8L   # metabolomic features

mk <- function(ids, p, prefix) {
  m <- matrix(rnorm(length(ids) * p), nrow = length(ids), ncol = p)
  colnames(m) <- sprintf("%s%02d", prefix, seq_len(p))
  data.frame(id = ids, m, check.names = FALSE)
}

X <- mk(idX, p1, "G")  # genomic
Z <- mk(idZ, p2, "P")  # proteomic
U <- mk(idU, p3, "M")  # metabolomic

## ---- clinical covariates (available for everyone) -------------------------
CC <- data.frame(id = id, age = rnorm(N), sex = rnorm(N), stage = rnorm(N),
                 check.names = FALSE)

## ---- true signal ----------------------------------------------------------
truth <- list(
  genomic     = c(1L, 2L, 3L),  # effects  1.5, -1.5,  1.0
  proteomic   = c(1L, 2L),      # effects  1.2, -1.2
  metabolomic = c(1L)           # effect   1.0
)
bG <- c(1.5, -1.5, 1.0)
bP <- c(1.2, -1.2)
bM <- c(1.0)
bCov <- c(0.3, -0.2, 0.25)

## latent linear predictor; a missing platform contributes nothing
eta <- numeric(N)
eta[match(idX, id)] <- eta[match(idX, id)] +
  as.matrix(X[, 1 + truth$genomic]) %*% bG
eta[match(idZ, id)] <- eta[match(idZ, id)] +
  as.matrix(Z[, 1 + truth$proteomic]) %*% bP
eta[match(idU, id)] <- eta[match(idU, id)] +
  as.matrix(U[, 1 + truth$metabolomic]) %*% bM
eta <- eta + as.matrix(CC[, -1]) %*% bCov
eta <- as.numeric(scale(eta))  # centre/scale the signal

## ---- outcomes -------------------------------------------------------------
## continuous (Gaussian)
ycont <- eta + rnorm(N, sd = 0.5)
outcome.continuous <- data.frame(id = id, y = ycont)

## binary (probit-style latent threshold)
ybin <- as.integer((eta + rnorm(N)) > 0)
outcome.binary <- data.frame(id = id, y = ybin)

## right-censored survival (accelerated failure time)
logT <- eta + rnorm(N, sd = 0.5)
time <- exp(logT)
cens <- exp(rnorm(N, mean = mean(logT) + 0.5, sd = 0.7))
status <- as.integer(time <= cens)             # 1 = event, 0 = censored
obstime <- pmin(time, cens)
outcome.survival <- data.frame(id = id, time = obstime, status = status)

simIMR <- list(
  platforms = list(genomic = X, proteomic = Z, metabolomic = U),
  covariates = CC,
  outcome = outcome.binary,            # canonical demo outcome (binary)
  outcome.binary = outcome.binary,
  outcome.continuous = outcome.continuous,
  outcome.survival = outcome.survival,
  truth = truth
)

dir.create("data", showWarnings = FALSE)
save(simIMR, file = "data/simIMR.rda", compress = "xz")
cat("Wrote data/simIMR.rda\n")
cat("Event rate (survival):", round(mean(status), 2),
    "  Positive rate (binary):", round(mean(ybin), 2), "\n")
