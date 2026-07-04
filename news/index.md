# Changelog

## funconstrain 0.1.0

### Bug fixes and minor improvements

- Big spring (summer) clean to bring repo up to some modern standards.
- Fixed narrow Hessian edge cases in
  [`bard()`](https://jlmelville.github.io/funconstrain/reference/bard.md),
  [`meyer()`](https://jlmelville.github.io/funconstrain/reference/meyer.md),
  and
  [`kow_osb()`](https://jlmelville.github.io/funconstrain/reference/kow_osb.md).
- Removed a duplicate Hessian assignment in
  [`rosen()`](https://jlmelville.github.io/funconstrain/reference/rosen.md).
- Improved documentation and tests.

## funconstrain 0.0.0.9004 and earlier:

### 2024-04-09

- Added `inst/doc/Examples/` with a program to apply `optimx` solvers to
  functions from funconstrain.
- Added `optimx` to `Suggests` in `DESCRIPTION`.
- Compacted the example `RunFunconstrainTests.Rmd` output with
  `tools::compactPDF("RunFunconstrainTests.pdf", gs_quality = "ebook")`.

### 2022-11-24

- Added `fmin` and `xmin` to each result to give an example solution for
  each problem. Some problems have variable dimensions and inputs, or
  may have multiple minima.

### 2022-05-19

- Added a `he()` function to each case to provide the Hessian at given
  parameters.
