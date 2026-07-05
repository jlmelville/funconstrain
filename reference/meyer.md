# Meyer Function

Test function 10 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
meyer()
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
  `m = 16`.

- Minima: the MGH (1981) only provides the optimal value, with
  `f = 87.9458...`. Meyer and Roth (1972) give the optimal parameter
  values as `(0.0056, 6181.4, 345.2)`, with `f = 88`.

## Note

The gradient is large even at the optimal value, and enormous if using
the level of precision given by Meyer and Roth (smallest gradient
component is at least 1e4). It is not recommended to rely on the typical
gradient norm termination conditions if using this test function.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41. <https://doi.org/10.1145/355934.355936>

Meyer, R. R. (1970). Theoretical and computational aspects of nonlinear
regression. In J. B. Rosen, O. L. Mangasarian, and K. Ritter (Eds.)
*Nonlinear programming* (pp465-496). New York: Academic Press.

Meyer, R. R., & Roth, P. M. (1972). Modified damped least squares: an
algorithm for non-linear estimation. *IMA Journal of Applied
Mathematics*, *9*(2), 218-233. <https://doi.org/10.1093/imamat/9.2.218>

## Examples

``` r
fun <- meyer()
# Optimize using the standard starting point
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B")
# Use your own starting point
res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")
```
