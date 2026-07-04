# Watson Function

Test function 20 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
watson()
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

- Dimensions: Number of parameters `2 <= n <= 31`, number of summand
  functions `m = 31`.

- Minima: `f = 2.28767...e-3` if `n = 6`; `f = 1.39976...e-6` if
  `n = 9`; `f = 4.72238...e-10` if `n = 12`

The number of parameters, `n`, in the objective function is not
specified when invoking this function. It is implicitly set by the
length of the parameter vector passed to the objective and gradient
functions that this function creates. See the 'Examples' section.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41.
[doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936)

Kowalik, J. S., & Osborne, M. R. (1968). *Methods for unconstrained
optimization problems.* New York, NY: Elsevier North-Holland.

## Examples

``` r
wat <- watson()
# 6 variable problem using the standard starting point
x0_6 <- wat$x0(6)
res_6 <- stats::optim(x0_6, wat$fn, wat$gr, method = "L-BFGS-B")
# Standard starting point with 9 variables
x0_9 <- wat$x0(9)
res_9 <- stats::optim(x0_9, wat$fn, wat$gr, method = "L-BFGS-B")
# Create your own 3 variable starting point
res_3 <- stats::optim(c(0.1, 0.2, 0.3), wat$fn, wat$gr, method = "L-BFGS-B")
```
