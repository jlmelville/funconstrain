# funconstrain 0.1.1

## Bug fixes and minor improvements

- Added `funconstrain_catalog()` and `funconstrain_problem()` for ordered
  problem discovery and strict, solver-neutral resolution of concrete test
  configurations.
- Fixed gradient edge cases in `brown_al()` and `helical()`.
- Corrected the standard starting points for `disc_bv()` and `disc_ie()`.
- Improved dimension and input validation, including consistent errors where
  Hessians are undefined.
- Extended `chebyquad()` to support the full `m >= n` problem family while
  retaining `m = n` when `m` is omitted.
- Updated the mize integration harness and article for stable mize 0.3.0;
  release installs no longer need GitHub commit metadata.
- Removed the duplicate `"L-BFGS-B"` entry from the optimizer methods returned
  by `fufn()`.
- `fufn()` now rejects malformed, non-whole, and out-of-range problem numbers
  with clear validation errors.
- Documentation has been clarified for problems with configurable dimensions and the reported minima.

# funconstrain 0.1.0

## Bug fixes and minor improvements

- Big spring (summer) clean to bring repo up to some modern standards.
- Fixed narrow Hessian edge cases in `bard()`, `meyer()`, and `kow_osb()`.
- Removed a duplicate Hessian assignment in `rosen()`.
- Improved documentation and tests.

# funconstrain 0.0.0.9004 and earlier:

## 2024-04-09

- Added `inst/doc/Examples/` with a program to apply `optimx` solvers to functions from 
funconstrain.
- Added `optimx` to `Suggests` in `DESCRIPTION`.
- Compacted the example `RunFunconstrainTests.Rmd` output with
  `tools::compactPDF("RunFunconstrainTests.pdf", gs_quality = "ebook")`.

## 2022-11-24

- Added `fmin` and `xmin` to each result to give an example solution for each problem. Some
  problems have variable dimensions and inputs, or may have multiple minima.

## 2022-05-19

- Added a `he()` function to each case to provide the Hessian at given parameters. Thank you to
  [John C. Nash](https://github.com/nashjc) for the contribution.
