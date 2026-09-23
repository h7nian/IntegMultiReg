# Independent JSS/code review — 2026-09-09

## Fix follow-up

Both findings below have now been fixed in the local working tree:

- `imr_data()` and `validate_imr_data()` reject reserved platform names `id` and
  `subgroup` through a shared check. Roxygen/Rd documentation and NEWS explain
  the restriction. New regressions cover construction, saved objects, CV,
  prediction and a successful ordinary-name CV control.
- The binary appendix now specifies the weighted sum of model-specific probit
  probabilities, and distinguishes conditional posterior prediction. Both local
  manuscript sources and the manuscript PDFs in paper/, overleaf/ and output/pdf/
  have been synchronized. Pages 28–29 were rendered and visually checked.
- The full package-context test suite passes 479 expectations (12 new), with no
  failures, warnings or skips. LaTeX compiles with no undefined references or
  citations; pre-existing overfull-box messages elsewhere remain. Logs are in
  `output/jss-review-fixes-2026-09-09/`.

The review below records the pre-fix findings and evidence. Existing release
archives and previous CI results remain historical; this follow-up does not
claim a new CI run, remote push or online Overleaf comment resolution.

## Conclusion

The seven concerns in `Review/JSS Comments.rtf` have substantive responses in
the current 0.1.4 implementation and submission materials. However, complete
closure is premature: the two reproducible findings below remain. This review
does not modify the implementation, manuscript, remote branches or review status.

## Findings

### P2 — Platform names can collide with availability metadata

Location: `IntegMultiReg/R/imr-data.R:44–53, 97–102`.

The constructor only checks uniqueness between platform names. A platform named
`subgroup` produces two availability columns with that name: the logical presence
column and the subgroup bitstring. Validation reconstructs the same ambiguous
table and accepts it. `cv_imr()` reads the first column with `$subgroup`, loses
the eligible subjects and incorrectly reports that k exceeds subgroup size.
A platform named `id` similarly creates duplicate metadata column names.

Reproduced with the freshly compiled workspace package:

```r
library(IntegMultiReg)
x <- data.frame(id = 1:20, marker = sin(1:20))
y <- data.frame(id = 1:20, y = cos(1:20))
d <- imr_data(list(subgroup = x), y, type_outcome = "continuous")
names(d$availability) # "id" "subgroup" "subgroup"
validate_imr_data(d)  # accepted
f <- imr(d, ssize = 1, h0 = 1, sample_mcmc = c(10, 5), seed = 1)
cv_imr(f, k = 2, rounds = 1) # erroneous subgroup-size error
```

The otherwise identical example named `assay` succeeds. Reject reserved names
in construction and validation, or keep presence columns in a namespace that
cannot collide with metadata. Add a regression for constructor/validator/CV.

### P2 — Binary appendix still describes the corrected-away algorithm

Location: `paper/IntegMultiReg.tex:1554–1556`; the local Overleaf copy is identical.

The appendix says to transform the averaged latent predictor through the probit
link. The actual implementation, `src/prediction_cv.c:118`, and manuscript main
text correctly average model-specific probabilities. In general,
`sum(w * pnorm(eta)) != pnorm(sum(w * eta))`. This is directly relevant to the
editor's numerical-reproducibility concern and contradicts the response letter's
claim that the methodological description has been synchronized.

Replace the appendix description with the weighted model-specific probability
formula, distinguishing the point-prediction approximation from conditional
posterior prediction, then regenerate the manuscript PDF.

## Editor-comment disposition

| Concern | Review result |
| --- | --- |
| Validated data class and reusable prediction inputs | Implemented, but reserved-name validation defect remains above. |
| Fit inspection/comparison, traces, uncertainty, covariate choice | Implemented: diagnostics, posterior methods and paired/nested formula-comparison example. Conditional-chain diagnostics and Laplace approximation limits are stated. Clinical variables remain forced within a candidate formula. |
| Replication contains manuscript commands and prints results | Nine CodeInput blocks are represented in the current script/companion; both fit-printing calls and displayed prediction output are present in full-run logs. |
| Prediction-table discrepancy | Current versioned results match the two saved complete runs; CV independently refits each training fold. Binary appendix still needs the correction above. |
| Related software outside R | A cross-environment table covers R/Python MOFA2, Python mvlearn, MATLAB/R SNF and Julia CCA with cited sources. Source currency was not independently re-researched in this review. |
| Full prediction precision | No rounding in returned point predictions; probability-averaging regression exists. |
| Formula/data interface | Implemented with model.frame/model.matrix, stored training terms, levels and contrasts, and training-fold refits. |

## Independently checked evidence

- Built and installed the current working-tree source into a fresh temporary R library.
- Independently reran `testthat::test_local("IntegMultiReg")`: all 467 expectations
  passed. Log: `/private/tmp/imr-review-sep09-tests-local.log`. An initial bare
  `test_dir()` invocation failed to resolve one package-internal helper because
  it lacked the package test environment; the standard package-context run passed.
- Compared all 14 CSV outputs directly: both saved complete runs and all final
  reference CSVs are byte-identical. This checks saved evidence; a new complete
  manuscript-scale simulation was not run in this review.
- Current R, C and test sources match the final distribution tarball. Differences
  found under src/inst are a Makevars.win final newline and generated inst/doc files.
- Current manuscript TeX matches its local Overleaf copy.
- Public GitHub API confirms revision branch SHA
  `e59ed766d8916dc05dd6fd040d1b490d084ff385` and successful completed runs
  [Linux release](https://github.com/h7nian/IntegMultiReg/actions/runs/34309998775)
  and [cross-platform](https://github.com/h7nian/IntegMultiReg/actions/runs/34310000246)
  at that SHA. No authentication token was needed for these public reads.
- Existing complete R CMD check log reports Status OK; existing package-context
  test log records 467 passing expectations with no failures, warnings or skips.

## Scope and version cautions

The local package checkout has uncommitted 0.1.4 changes on top of an older
commit; remote revision evidence must be tied to source contents, not local HEAD.
The August review-change record is historical and still describes the old
post-fitting CV approach. It also says Overleaf review entries were deliberately
left unresolved for author inspection. This review does not establish the live
Overleaf comment-resolution state or editorial acceptance.
