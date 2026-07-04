# Jennrich and Sampson Function

Test function 6 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
jenn_samp(m = 10)
```

## Arguments

- m:

  Number of summand functions in the objective function. Should be equal
  to or greater than 2.

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

- Dimensions: Number of parameters `n = 2`, number of summand functions
  `m >= n`.

- Minima: `f = 124.362...` at `(x1 = x2 = 0.2578)` for `m = 10`,

## Note

This test problem isn't really unconstrained. `x1` must take a value
between `(-1, 1)`. Included for the sake of completeness.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41.
[doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936)

Jennrich, R. I., & Sampson, P. F. (1968). Application of stepwise
regression to non-linear estimation. *Technometrics*, *10*(1), 63-72.

## Examples

``` r
# Use m = 10 summand functions
fun <- jenn_samp(m = 10)
# Optimize using the standard starting point
# Set 'lower' and 'upper' parameter to constrain par[1]. Only works with
# L-BFGS-B.
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B", lower = -1, upper = 1)
# Use your own starting point
res <- stats::optim(c(0.1, 0.2), fun$fn, fun$gr, method = "L-BFGS-B",
lower = -1, upper = 1)

# Use 20 summand functions
fun20 <- jenn_samp(m = 20)
res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B",
lower = -1, upper = 1)
```
