# funconstrain

An R Package of Functions for Testing Unconstrained Numerical Optimization.

[![R-CMD-check](https://github.com/jlmelville/funconstrain/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jlmelville/funconstrain/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/jlmelville/funconstrain/graph/badge.svg?token=eR44zzwo9V)](https://app.codecov.io/gh/jlmelville/funconstrain)

`funconstrain` is a pure R implementation of the 35 test functions in the
paper by
[Moré, Garbow, and Hillstrom](https://doi.org/10.1145/355934.355936), useful
(to varying degrees) for testing unconstrained numerical optimization methods
such as steepest descent, Newton, BFGS, L-BFGS, and conjugate gradient.

## Install

```R
# install.packages("pak")
pak::pak("jlmelville/funconstrain")
```

## Documentation

```R
package?funconstrain
```

The [getting-started
article](https://jlmelville.github.io/funconstrain/articles/getting-started.html)
shows how to discover problems, resolve dimensions, run an optimizer, and
interpret stored references.

## Quick start

Use `funconstrain_problem()` to resolve a named problem to one numeric starting
point, dimension-checked callbacks, and metadata for the effective
configuration:

```R
library(funconstrain)

problem <- funconstrain_problem("rosen")
result <- stats::optim(
  par = problem$x0,
  fn = problem$fn,
  gr = problem$gr,
  method = "BFGS"
)

result$convergence
problem$fn(result$par)
max(abs(problem$gr(result$par)))
```

Direct factories remain available when you want the lower-level interface:

```R
raw_problem <- rosen()
names(raw_problem)
```

See `?funconstrain_catalog` for the resolver contract and
`package?funconstrain` for the raw factory contract and full function list.

## Why do this?

For testing numerical optimization routines, the go-to set of test problems is
[CUTEst](https://github.com/ralna/CUTEst). However, if you aren't compiling
and linking to it directly, you'd have to write a parser for the SIF file format it uses. Neither
of those possibilities appealed.

Instead, I re-implemented the functions as provided in the paper and also calculated analytical
gradients for them. I did look to see if all 35 problems were implemented in one place in R, but
failed to find such a package.

[John Nash](https://github.com/nashjc) contributed the Hessians.

## Are the functions correct?

There are unit tests for each test problem which ensure that:

- The analytical gradients match finite difference estimates at the suggested starting point.
- The Hessians have the expected shape and symmetry at the suggested starting point.
- If the location of a minima was given in the paper, that the analytical gradient is close to zero
at that location.
- If the location of a minima was given in the paper, that the objective function has the correct
value at that location.
- Running the L-BFGS or BFGS method as implemented in the `stats::optim` function gets to the
specified minimum.
- *July 7 2026*: I verified the gradient and Hessian values at the test locations via an independent
implementation in PyTorch, comparing the analytic results with the autodiff output.

## Is the implementation efficient?

Not really. My goal was correctness, and to make the code clear. Also perhaps, to be useful if
anyone ever wants to translate these into other languages without having to know a lot of idiomatic
R.

I have made use of vectorized arithmetic operation rather than explicit `for` loops where possible
as well as using functions like `sum`. Also, I am pretty profligate in storing pre-computed
vectors, trading off memory consumption for clarity and potentially fast vectorized computations (I
have not done any profiling). But I consciously eschewed the use of  `apply` `sweep` or other
cleverness.

I think I have elided the most gratuitous inefficiencies, such as unnecessary recomputation of
values inside loops.

## See also

- The aforementioned [CUTEst](https://github.com/ralna/CUTEst). I believe all
or nearly all of the test problems in this package are implemented in CUTEst, but I make no
representation that you will get the same results (if there are any differences, assume it's a bug
in `funconstrain`).
- I made ample use of the excellent [Derivative Calculator](https://www.derivative-calculator.net/)
to calculate the analytical gradients.
- Shameless plug: I wrote this package to test [mize](https://github.com/jlmelville/mize), an R
package for doing numerical optimization.

## License

[MIT License](https://opensource.org/licenses/MIT).
