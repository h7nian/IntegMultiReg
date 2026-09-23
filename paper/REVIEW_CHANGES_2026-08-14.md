# IntegMultiReg reviewer-change record — 2026-08-14

This record maps the 19 unresolved Overleaf Review entries on
`IntegMultiReg.tex` (17 comments and 2 tracked edits) to manuscript changes.
The project Overview also lists one tracked addition in
`Archives/Recovery2026August13.tex`, for 20 project-level entries in total.
All Review entries were deliberately left unresolved so the authors can inspect
and close them individually.

## Recovery points

- Pre-review manuscript commit: `a0e5900`.
- Merge containing the latest web-side Overleaf edits: `4770d3b`.
- Exact pre-revision source: `Archives/IntegMultiReg-before-review-2026-08-14.tex`.
  SHA-256: `cfdc2ccb6bfe61389b3eaee3bd560a379bf1f78a79d1cf08c9f64ec3db1c0548`.
- Exact pre-revision bibliography:
  `Archives/IntegMultiReg-before-review-2026-08-14.bib`.
  SHA-256: `38ddc37ddcf787a51d901efc6f15516530f6551a608547b49db745858d193eda`.
- `Archives/Recovery2026August13.tex`, already created in Overleaf, is also
  retained unchanged.

## Comment-by-comment disposition

1. **Cite `spikeSlabGAM`.** Added Scheipl (2011) and described its
   single-design structured additive spike-and-slab target in the Introduction.
2. **Discuss MOFA, BIPnet, SIDA and RandMVLearn.** Added all four, cited their
   primary papers, and stated the precise distinction: none implements IMR's
   availability-subgroup outcome regressions with MRF-linked feature indicators.
   The text explicitly notes that MOFA can accommodate incomplete entries but
   is unsupervised.
3. **CRAN availability.** Added the canonical CRAN package URL in the
   Introduction.
4. **Define beta and beta-prime in Equation 1.** Replaced the ambiguous symbols
   with `beta_C` (always-included clinical covariates) and `beta_X` (candidate
   molecular features), and defined both immediately after the equation. The
   appendix now uses `beta_C` consistently.
5. **Introduce the outcome list.** Added “Depending on the outcome type,
   y-star is defined as follows.”
6. **Define t.** Defined `t_n^(s)` as the observed event or censoring time for
   subject `n` in subgroup `s`.
7. **Clarify borrowing-strength interaction matrices.** Expanded the
   `theta_mean` object-table entry to identify `theta_ss'^(k)`, the per-platform
   matrix layout, and the pair of linked subgroups.
8. **Clarify responses y-star.** Expanded `estimate_latent_y` to distinguish
   working responses, augmented log-times, and latent Gaussian utilities.
9. **Configuration/model terminology.** Standardized posterior gamma objects
   to “selection model”; retained “configuration” only where it denotes a
   general run/setup rather than a gamma state.
10. **`cv_imr()` matches the paper.** Corrected the text: it starts from one
    full-data fit, re-estimates coefficients in training folds for sampled
    selection models, and combines held-out predictions with importance
    weights; the full MCMC sampler is not rerun within each fold.
11. **Remove the claim that the paper used complete/nested CV.** Deleted that
    claim and the recommendation that nested CV is required to match the paper.
12. **Compare refitting and post-fitting procedures.** Added that full refitting
    avoids importance-sampling approximation error when computation permits,
    whereas `cv_imr()` is the faster paper-aligned procedure.
13. **Table 6 zero SEs and numerical differences.** Initially replaced the
    two-decimal values with three-decimal values from
    `figures/table1_supplement_long.rds`, making the underlying variation
    visible. Added the exact long-chain setting and clarified that this package
    reanalysis is not a literal duplication of the historical study-specific
    C program. Item 19 records the subsequent correction from SE to SD.
14. **Re-check reproducibility statements.** Stated explicitly that each
    IMR/BMS specification uses one 400,000-retained/50,000-burn-in fit followed
    by 10 rounds of 10-fold post-fitting CV; the MCMC fit is not repeated for
    all 100 fold assignments.
15. **Platform-specific hyperparameters already exist.** Removed that proposed
    extension and explicitly acknowledged the existing platform-specific
    sparsity vector `nu_k`.
16. **Add Dai et al. (2023).** Added the full citation and identified its
    feature-ordering/population-structure construction as a possible extension.
17. **Binary latent-variable precedent.** Retained the existing
    Chekouo et al. AOAS citation and clarified Equation 10's exponential
    parameterization as a one-sided proposal with mean `1/abs(z)`.
18. **Do not overstate the historical KIRC comparison.** Reframed the long-chain
    Table 6 analysis as a package-based reanalysis rather than an exact
    reproduction. The manuscript now reports the mean and maximum absolute
    C-index differences from the published table.
19. **Move the complete KIRC experiment and report SD.** Moved the entire
    TCGA-KIRC real-data experiment to a dedicated appendix: the complete-case
    sample-retention comparison, data provenance, the reduced runnable package
    example, subgroup construction, survival-model code, long-chain setup,
    comparison table and interpretation. The main text now retains only a
    concise pointer to the appendix and develops the simulation experiment.
    Replaced the displayed
    SEs by sample SDs across the ten CV rounds, and documented the two
    implementation differences identified by the original study author:
    subgroup-specific CV partitions and a different parameter-sampling order
    from the historical C program.
20. **Use a single Appendix A numbering scheme.** Combined the supplementary
    empirical analyses and outcome extensions under one Appendix A. Appendix
    headings now use A.1--A.5, appendix tables use Table A.1--A.2, and
    appendix equations use (A.1) onward; all cross-references remain automatic.

## Web-side and tracked edits retained and reconciled

- The web-side source change “via importance sampling,” already present in the
  latest Overleaf source, is incorporated into the corrected `cv_imr()`
  description.
- The web-side computational-scale edit is retained and expanded to distinguish
  post-fitting coefficient/model-averaging cost from fully nested MCMC cost.
- The web-side “cross-validation” wording is retained in the revised section.

## Historical KIRC table audit evidence

- Source artifact: `figures/table1_supplement_long.rds`.
- MCMC: 400,000 retained draws plus 50,000 burn-in draws.
- Validation: 10 rounds of 10-fold post-fitting CV, seed 1.
- Final parenthetic values are sample SDs across the ten CV rounds, not SEs.
  They are obtained from the same saved results as `SE * sqrt(10)` and range
  from approximately 0.003 to 0.030.
- The underlying SDs are nonzero in both E1 (approximately 0.008--0.015) and
  E2 (approximately 0.017--0.030).
- Compared with the published table, the mean absolute C-index difference is
  0.023 and the maximum is 0.062 (BMS + C + M in E2). The table is therefore
  described as a package-based reanalysis, not an exact numerical reproduction.
- The Biometrics table described its parenthetic values as SEs, but the original
  study author confirmed that they were SDs and that the label was a typo. The
  manuscript and replication code now use SD consistently.

## Validation checklist

- [x] BibTeX resolves all six newly added citations.
- [x] LaTeX compiles without errors.
- [x] No undefined citations or references remain.
- [x] The long-chain KIRC table fits within the page width.
- [x] Revised PDF is visually inspected (Introduction, model definitions,
  fitted-object table, CV section, KIRC table, Summary, bibliography, and
  binary-outcome appendix).
- [x] Revised source and this record are committed and pushed to Overleaf.
- [ ] Authors review and resolve Overleaf comments individually.
