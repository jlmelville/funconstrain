# Biggs EXP6 Function

Test function 18 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
biggs_exp6(m = 13)
```

## Arguments

- m:

  Number of summand functions in the objective function. Should be equal
  to or greater than 6.

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

- Dimensions: Number of parameters `n = 6`, number of summand functions
  `m >= n`.

- Minima: `f = 5.65565...e-3` if `m = 13`; not reported in the
  MGH (1981) paper is `(f = 0)` at `c(1, 10, 1, 5, 4, 3)` and
  `c(4, 10, 3, 5, 1, 1)` (and probably others) for all `m` (probably: I
  stopped testing after `m = 1000`).

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41. <https://doi.org/10.1145/355934.355936>

Biggs, M. C. (1971). Minimization algorithms making use of non-quadratic
properties of the objective function. *IMA Journal of Applied
Mathematics*, *8*(3), 315-327.

## Examples

``` r
fun <- biggs_exp6()
# Optimize using the standard starting point
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B")
# Use your own starting point
res <- stats::optim(1:6, fun$fn, fun$gr, method = "L-BFGS-B")

# Use 20 summand functions
fun20 <- biggs_exp6(m = 20)
res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
```
