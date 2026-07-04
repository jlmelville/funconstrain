# Variably Dimensioned Function

Test function 25 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
var_dim()
```

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
[`linfun_r1z()`](https://jlmelville.github.io/funconstrain/reference/linfun_r1z.md),
and
[`trigon()`](https://jlmelville.github.io/funconstrain/reference/trigon.md).
It is not part of the core return contract needed by optimizers.

## Details

The objective function is the sum of `m` functions, each of `n`
parameters.

- Dimensions: Number of parameters `n` variable, number of summand
  functions `m = n + 2`.

- Minima: `f = 0` at `rep(1, n)`.

The number of parameters, `n`, in the objective function is not
specified when invoking this function. It is implicitly set by the
length of the parameter vector passed to the objective and gradient
functions that this function creates. See the 'Examples' section.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41.
[doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936)

## Examples

``` r
vdim <- var_dim()
# 6 variable problem using the standard starting point
x0_6 <- vdim$x0(6)
res_6 <- stats::optim(x0_6, vdim$fn, vdim$gr, method = "L-BFGS-B")
# Standard starting point with 8 variables
res_8 <- stats::optim(vdim$x0(8), vdim$fn, vdim$gr, method = "L-BFGS-B")
# Create your own 4 variable starting point
res_4 <- stats::optim(c(0.1, 0.2, 0.3, 0.4), vdim$fn, vdim$gr,
                      method = "L-BFGS-B")
```
