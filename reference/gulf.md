# Gulf Research and Development Function

Test function 11 from the Moré, Garbow and Hillstrom paper.

## Usage

``` r
gulf(m = 99)
```

## Arguments

- m:

  Number of summand functions in the objective function. Should be
  between 3 and 100, according to the MGH paper. Default value is 99,
  which Jamil and Xang (2013) list as the only valid value.

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
  `3 <= m <= 100`.

- Minima: `f = 0` at `(50, 25, 1.5)`.

Note that the equation as published by Moré, Garbow and Hillstrom (1981)
contains an error, where the symbol 'mi' should be interpreted as a
minus sign. The corrected version can be found in Jamil and Xang (2013),
and no doubt several other publications. The Jamil and Xang equation
unfortunately contains its own minor errors, but you can piece together
the correct equation from these two sources without too much trouble.

## References

Moré, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing
unconstrained optimization software. *ACM Transactions on Mathematical
Software (TOMS)*, *7*(1), 17-41. <https://doi.org/10.1145/355934.355936>

Jamil, M., & Yang, X. S. (2013). A literature survey of benchmark
functions for global optimisation problems. *International Journal of
Mathematical Modelling and Numerical Optimisation*, *4*(2), 150-194.
<https://doi.org/10.1504/IJMMNO.2013.055204>
<https://arxiv.org/abs/1308.4008>

## Examples

``` r
# Use 10 summand functions
fun <- gulf(m = 10)
# Optimize using the standard starting point
x0 <- fun$x0
res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
"L-BFGS-B")
# Use your own starting point
res <- stats::optim(c(0.1, 0.2, 0.3), fun$fn, fun$gr, method = "L-BFGS-B")

# Use 20 summand functions
fun20 <- gulf(m = 20)
res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
```
