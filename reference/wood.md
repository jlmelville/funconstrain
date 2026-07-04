# Wood Function

Test function 14 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
wood()
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

- Dimensions: Number of parameters `n = 4`, number of summand functions
  `m = 6`.

- Minima: `f = 0` at `rep(1, 4)`.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41.
[doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936)

Colville, A. R. (1968) *A comparative study on nonlinear programming
codes* (Report No. 320-2949). New York, NY: IBM.

## Examples

``` r
fun <- wood()
# Optimize using the standard starting point
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B")
# Use your own starting point
res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
"L-BFGS-B")
```
