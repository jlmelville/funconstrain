# Linear Function - Rank 1 with Zero Columns and Rows

Test function 34 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
linfun_r1z(m = 100)
```

## Arguments

- m:

  Number of summand functions in the objective function. Should be equal
  to or greater than `n`.

## Value

A list containing the core problem contract:

- `fn`: Objective function. It takes a numeric parameter vector and
  returns the scalar objective value.

- `gr`: Gradient function. It takes a numeric parameter vector and
  returns the gradient vector.

- `he`: Hessian function. It takes a numeric parameter vector and
  returns the Hessian matrix of second derivatives.

- `fg`: Combined objective and gradient function. It takes a numeric
  parameter vector and returns a list with members `fn` and `gr`.

- `x0`: Suggested starting point. For fixed-dimension problems this is a
  numeric vector; for variable-dimension problems this is a function
  that returns a numeric vector for a requested `n`.

- `fmin`: Reported minimum objective value.

- `xmin`: Numeric vector at a reported minimum.

Some factories also include `m`, a metadata field for the number of
summand functions. It is absent for most factories, `NA` for several
legacy fixed problems, and numeric for
[`linfun_fr()`](https://jlmelville.github.io/funconstrain/reference/linfun_fr.md),
[`linfun_r1()`](https://jlmelville.github.io/funconstrain/reference/linfun_r1.md),
`linfun_r1z()`, and
[`trigon()`](https://jlmelville.github.io/funconstrain/reference/trigon.md).
It is not part of the core return contract needed by optimizers.

## Details

The objective function is the sum of `m` functions, each of `n`
parameters.

- Dimensions: Number of parameters `n` variable, number of summand
  functions `m >= n`.

- Minima: `f = (m * m + 3 * m - 6) / (2 * (2 * m - 3))` at any set of
  points `x[j]` with `j = 2, ..., n - 1` where the sum of
  `j * x[j] = 3 / (2 * m - 3)`.

The number of parameters, `n`, in the objective function is not
specified when invoking this function. It is implicitly set by the
length of the parameter vector passed to the objective and gradient
functions that this function creates. See the 'Examples' section.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41. <https://doi.org/10.1145/355934.355936>

## Examples

``` r
linr1z <- linfun_r1z(m = 10)
# 6 variable problem using the standard starting point
x0_6 <- linr1z$x0(n = 6)
res_6 <- stats::optim(x0_6, linr1z$fn, linr1z$gr, method = "L-BFGS-B")
# Standard starting point with 8 variables
res_8 <- stats::optim(linr1z$x0(8), linr1z$fn, linr1z$gr, method =
"L-BFGS-B")
# Create your own 4 variable starting point
res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), linr1z$fn, linr1z$gr,
                      method = "L-BFGS-B")
# Use 20 summand functions
linr1z_m20 <- linfun_r1z(m = 20)
# Repeat 4 parameter optimization with new test function
res_n4_m20 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), linr1z_m20$fn,
linr1z_m20$gr, method = "L-BFGS-B")
```
