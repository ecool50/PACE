# Changelog

## PACE 0.99.18

- Condition signals build their messages without
  [`paste()`](https://rdrr.io/r/base/paste.html), and the ridge check
  handles empty input directly rather than suppressing the warning it
  would raise.

## PACE 0.99.17

- New
  [`compactFit()`](https://ecool50.github.io/PACE/reference/compactFit.md)
  drops the per-cell matrices before saving a fit, which the vignette
  previously did by assigning into a slot.
- The vignettes read fits through accessors rather than `@`, which
  clears BiocCheck’s note about slot access in examples and vignettes.

## PACE 0.99.16

- New
  [`stashPACE()`](https://ecool50.github.io/PACE/reference/stashPACE.md)
  copies a fit’s readouts into `metadata(spe)`, so results can travel
  with the data they came from. The solver state is not copied.

## PACE 0.99.15

- Records the Australian Government Research Training Program as funder
  in `Authors@R`.

## PACE 0.99.14

- `PACEFit` gains a validity method. Replacing a slot wholesale used to
  pass silently and break every readout downstream; `fit@fit <- list()`
  is now refused.
- New accessors
  [`cellTypes()`](https://ecool50.github.io/PACE/reference/cellTypes.md),
  [`paceParams()`](https://ecool50.github.io/PACE/reference/paceParams.md)
  and
  [`fitState()`](https://ecool50.github.io/PACE/reference/fitState.md),
  and validated replacement methods for
  [`neighbourSlopes()`](https://ecool50.github.io/PACE/reference/neighbourSlopes.md),
  [`varianceDecomposition()`](https://ecool50.github.io/PACE/reference/varianceDecomposition.md)
  and
  [`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md).

## PACE 0.99.13

- `verbose = FALSE` is now silent. The shrinkage step reported its
  progress regardless, because it never received the flag, and `mashr`
  printed its own on top of that.
- Installation instructions use
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).
- The vignette evaluates every chunk.
- Minimum R version is 4.6.0.

## PACE 0.99.12

- [`plotDecomposition()`](https://ecool50.github.io/PACE/reference/plotDecomposition.md)
  draws the zoom panel’s percentage labels upright, so that neighbouring
  bars of similar height no longer overlap each other’s text.

## PACE 0.99.11

- Added a canonical multivariate adaptive shrinkage implementation to
  the compiled core, alongside the existing `mashr` path. It reproduces
  `mashr`’s calls exactly on both reference cohorts.

## PACE 0.99.10

- The driver score columns are named for what the method calls them:
  `driver_score` and `driver_score_spec4`, previously `MCSD` and
  `MCSD4`.
  [`plotDrivers()`](https://ecool50.github.io/PACE/reference/plotDrivers.md)
  takes `panels = "score"` in place of `"mcsd"`.

## PACE 0.99.9

- [`topDrivers()`](https://ecool50.github.io/PACE/reference/topDrivers.md)
  documents which model term each column holds. On a fit with a
  condition the table carries two different slopes, and they are not a
  raw and shrunken pair of the same quantity.
- The condition vignette shows the responder decomposition and driver
  figures, and prints the driver table.

## PACE 0.99.6

- The per-pair gene score is called the driver score throughout the
  documentation and plot labels, matching the method description.

## PACE 0.99.4

- Fixed portability defects found by the Bioconductor build matrix: bit
  identity on arm64, a crash on Windows from `static thread_local`
  destructors, an unsupported `system2(env=)` on Windows, and missing
  `-pthread` on compile flags for Wasm.

## PACE 0.99.0

- Initial Bioconductor submission.
