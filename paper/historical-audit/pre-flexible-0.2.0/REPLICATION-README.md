# Version 0.1.4 revision — September 8, 2026

The current `paper/IntegMultiReg-replication.R` requires version 0.1.4.
It independently refits every CV training fold, corrects binary probability
averaging and fixes concordance. Its CV values cannot be compared against the
0.1.3 references. `paper/IntegMultiReg-replication/` remains the frozen 0.1.3
archive, including historical checksums and validation records.

For the current workflow:

```sh
R CMD INSTALL IntegMultiReg_0.1.4.tar.gz
Rscript paper/IntegMultiReg-replication.R --quick
Rscript paper/IntegMultiReg-replication.R
```

The normal run compares against `paper/expected-results-0.1.4/`. To intentionally
create a new reference, run with `--update-reference`, then repeat independently
without that flag. Keep the source tarball and session information with the
reference. Short runs never create numerical references. The 400000-draw KIRC
appendix remains a historical 0.1.3 audit record and has not been regenerated
with training-fold refits. It must not be described as current validation.

The material below documents the frozen 0.1.3 archive.

---

# Replication materials for the IntegMultiReg manuscript

These files reproduce the empirical results in:

> Zhang S, Wang J, Chekouo T. *IntegMultiReg: An R Package for Integrative
> Bayesian Multi-Regression of Multi-Platform Biomarkers*.

## Contents

- `IntegMultiReg-replication.R`: the single, commented standalone replication
  script.
- `IntegMultiReg_0.1.3.tar.gz`: the exact R package source used by the script.
- `data/kirc_table1_full.rda`: the converted Biometrics supplementary KIRC
  object used by the long-chain appendix reanalysis.
- `expected-results/`: compact reference tables from the manuscript run.

The KIRC data are replication data derived from the supplementary matrices for
Chekouo, Stingo, Doecke and Do (2017), DOI `10.1111/biom.12587`. They are kept
in this replication archive rather than in the CRAN package because the source
data have a separate distribution and provenance boundary.

## Requirements

- R 4.0.0 or later.
- GNU Scientific Library (GSL) 2.0 or later.
- The R package `survival` for the optional historical KIRC workflow.

Install the included package source before running the replication:

```sh
R CMD INSTALL IntegMultiReg_0.1.3.tar.gz
```

## Reproduce the manuscript simulation figures and tables

Run from the extracted archive directory:

```sh
Rscript IntegMultiReg-replication.R
```

This creates `replication-output/`, including the three manuscript figures,
CSV copies of the numerical tables, an RDS object containing all fitted results
and `sessionInfo.txt`.

For a short end-to-end code-path check that does not reproduce the manuscript
numbers:

```sh
Rscript IntegMultiReg-replication.R --quick
```

## Run the historical TCGA-KIRC package-based reanalysis

The appendix configuration is computationally expensive (400000 retained
MCMC draws, 50000 burn-in draws and 10 rounds of 10-fold cross-validation):

```sh
Rscript IntegMultiReg-replication.R --full-kirc
```

To verify the complete KIRC code path using short chains:

```sh
Rscript IntegMultiReg-replication.R --quick --smoke-kirc
```

The smoke run is only a software check; its numerical output must not be
compared with the manuscript table.

The parenthetic values in the KIRC output are sample standard deviations across
the ten cross-validation rounds. The original Biometrics article labelled its
parenthetic values as standard errors, but the original study author has
confirmed that they were standard deviations and that the SE label was a typo.

This is a package-based reanalysis rather than an exact numerical reproduction
of the historical C program. In particular, the package constructs
cross-validation partitions separately within each availability subgroup, and
its parameter-sampling order differs from the original C implementation. These
differences can change individual C-index values while preserving the broad
comparative trends.

## Reproducibility controls

All stochastic analyses use fixed seeds. The script requires the exact archived
package version and checks the normal-run simulation results against the
manuscript values; it stops if any reported value differs beyond the stated
tolerance. The long KIRC reanalysis results are
also supplied in `expected-results/table_kirc_supplement.csv` so that a
completed long run can be compared without rerunning it during a quick audit.

## Response scale and historical results

The 0.1.3 default simulation and reduced KIRC runs use log survival time,
matching the original Biometrics supplementary code. Their reference tables
were regenerated and checked with the exact included source package.

The optional `--full-kirc` / `--smoke-kirc` path is explicitly a **legacy
raw-time** reproduction (`survival_scale = "identity"`). Its previously
published long-chain table is retained as an audit record, not as evidence for
the corrected log-time AFT model. A new long-chain log-time analysis is pending.
The original 0.1.2 replication archive is preserved separately.

The coefficient/predictive extension is approximate: model weights retain
Laplace approximation error. See the package documentation and accompanying
validation report for scope and diagnostics.

## Optional coefficient and predictive interval run

After the default replication, run `Rscript IntegMultiReg-uncertainty.R
replication-output` (as one command). It uses the validated 5000 burn-in / 4000
retained iterations per conditional chain and seed 609. The reference run had
155 distinct subgroup models and maximum classical split R-hat 1.015019.
Inspect diagnostics on your own run; this does not diagnose the selection chain.
