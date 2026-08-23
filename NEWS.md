# PACE 0.99.1

Correctness fixes from an internal audit of the fitting engine. All of these
produced wrong numbers with no error raised. The analyses in the accompanying
manuscript were produced with the pre-fix engine and were checked to be
unaffected.

Reported output is unchanged. On the shipped breast cancer subset the shrunken
slopes move by at most `1.6e-06` and the local false sign rates by `3.2e-05`
across all 17,792 (gene, focal, neighbour) rows, with the same 60 calls at
`lfsr < 0.05`, no flipped calls, and no sign changes. The raw solver output
moves at the `1e-04` level, and moves *towards* the full-double reference
(fixed effects `8.4e-04 -> 5.3e-05` from it), because the ridge no longer
drops into single precision. Every other fix is exactly reproducing: with the
precision guard disabled the refit is identical to 0.99.0 (`max|diff| = 0`).

* An `NA` in any random-effect variable no longer corrupts the design matrix.
  `model.matrix()` drops such rows, which left the sparse assembly recycling
  covariates so that cells received other cells' values; this now stops with a
  message naming the column and the number of rows dropped.
* The condition interaction term is derived from `condition_col` instead of a
  hardcoded `"Responder"` prefix. Cohorts whose condition column had any other
  name silently produced an all-zero condition indicator, and therefore an
  identically zero responder spatial block, while the output looked complete.
  A `resp_term` that names no level of `condition_col` is now refused.
* Standard errors are computed on every iteration, so a fit that exits early at
  convergence no longer returns all-`NA` `se_B`/`se_U` (which the shrinkage
  layer then discarded). Convergence on a float interior takes one further pass
  in double so the retained standard errors are not float-quantised.
* Non-finite genes are counted and reported instead of being masked out of the
  convergence metric, and a fit that has lost genes is no longer reported as
  converged. In the streaming engine the NaN repair is re-checked after it runs,
  and genes that remain non-finite are reported on every iteration, including
  the last, where the interior is already double and the repair cannot help.
* The ridge is kept out of single precision. It is accumulated inside the float
  kernel, where the smallest representable increment on a diagonal of order
  `sum(w)` is `eps_float * sum(w)`; below that the regularisation is quantised
  away, at a threshold that falls as `1/n`, so the same model and design lose
  their ridge purely by growing the cohort. The interior is now forced to double
  whenever the ridge is small relative to the diagonal scale.
* The variance components are bounded from above by a new `tau_max` argument
  (default `100`). Every previous tau write was a `pmax` against a floor with no
  ceiling, so a column identified only by the ridge climbed by `mean(u^2)` per
  iteration while its inflated posterior variance fed the next update. The
  default sits roughly 6-30x above the largest component observed on healthy
  cohorts and just below `var(z)`, the variance of the working response that the
  components decompose; it does not bind on any fit in the test suite. A binding
  cap is reported, since it means the term is not identified by the data.
* The empirical-Bayes floor in the adaptive tau update is scaled by the median
  across genes rather than the mean, so a single runaway gene can no longer lift
  the floor for every gene in its row. The shrinkage target itself is unchanged.
* The per-iteration largest variance component is recorded in the fit history and
  printed alongside the block medians, since a runaway confined to one term does
  not move a median taken over the whole block.
* mash is now given the correlation between conditions under the null. The
  focal cell types are conditions estimated in the same per-gene solve, so
  their errors are correlated, and the previous default of `V = I` treated them
  as independent and mis-calibrated the local false sign rates. The estimate is
  taken over the conditions that vary: sparse focal-neighbour pairs are zeroed
  by the kernel drop without their column being removed, so those conditions
  have a constant z column, and `cov2cor()` on the full matrix returns `NaN`
  for every entry. Degenerate conditions are left uncorrelated, which is the
  correct limit, and the estimate falls back to independence when fewer than
  two conditions vary. Set `null_correlation = FALSE` in `paceShrink()` to
  restore the old behaviour.

  On the shipped breast cancer subset the estimated correlations reach 0.20 in
  absolute value, the shrunken slopes correlate at 0.9954 with the independent
  fit, and one call changes net at `lfsr < 0.05` (60 to 61: CD163 in macrophages
  near tumour and PECAM1 in B cells near T cells gained, ERBB2 in tumour near
  T cells lost, all three within 0.02 of the threshold). The 1,294 nominal sign
  changes are numerical noise: none is significant in either pass and their
  median absolute effect is 4e-21.
* Leave-one-patient-out variance decomposition now subsets the contamination
  offset matrix alongside `mu`. Leaving it at full length offset every cell
  after the dropped patient's first row, so each leave-one-out row was wrong
  while the full-data row was correct.
* The condition-by-proximity pair contributions are matched to the neighbours
  whose interaction terms are actually present, rather than to the first
  `n` entries of the neighbour list, and pair rows are labelled from the same
  matching vector. The previous indexing misaligned neighbour columns whenever
  an interaction term was missing.

# PACE 0.99.0

* Initial package: `paceFit()` for SpatialExperiment input, the `PACEFit` S4
  class, and the `neighbourSlopes()`, `varianceDecomposition()`, and
  `topDrivers()` accessors, wrapping the locked PACE statistical engine.
