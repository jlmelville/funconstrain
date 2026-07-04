# Bard Function

Test function 8 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
bard()
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

- Dimensions: Number of parameters `n = 3`, number of summand functions
  `m = 15`.

- Minima: `f = 8.214877e-3` at c(0.08241056, 1.133036, 2.343695) Solvers
  terminate with `f` near 17 for parameter 1 in 0.84 to 0.89
  approximately and large negative values of the other two parameters.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41.
[doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936)

Bard, Y. (1970). Comparison of gradient methods for the solution of
nonlinear parameter estimation problems. *SIAM Journal on Numerical
Analysis*, *7*(1), 157-186.
[doi:10.1137/0707011](https://doi.org/10.1137/0707011)

## Examples

``` r
fun <- bard()
# Optimize using the standard starting point
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B")
# Use your own starting point
res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
```
