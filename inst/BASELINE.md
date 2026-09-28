# 0.1.4 behavior baseline

The 0.2.0 cleanup started from Git commit `9eaa321`.

The source package built from that commit with `R CMD build
--no-build-vignettes --no-manual` has SHA256:

`7ef7466647a973fed96a3abdd4088ad663d7fd56b5226f5817f96ff0c5b3eae4`

Fixed-seed differential checks cover all three outcome types with IMR and BMS.
They compare inclusion probabilities, interaction means, latent responses,
log-posterior traces, and predictions. The MRF log-sum-exp rewrite changes only
round-off in log-posterior values (observed maximum: `1.14e-13`).
