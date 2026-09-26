# Current manuscript reference

Promoted after two complete runs of the same validated 0.2.0 candidate.
A established an isolated reference; B ran without `--update-reference`.
The six prediction configurations retain 4000 draws after 1000 burn-in and
use the paper sampler with 5-fold x 10-round refit CV. The 23 saved CSV/RDS
artifacts agree at tolerance 1e-12 (21 are R-identical); all three figures
render pixel-identically at 150 dpi. Independent score/fold audits passed.

Environment: R 4.4.2, macOS 26.3.1 ARM64, GSL 2.8. Checksums and evidence
locations are in provenance.json. This directory is a numerical reference,
not a release approval or proof of convergence. It does not replace or explain
the historical 0.1.4 values, or establish original-scale study reproduction.
