# Getting started with funconstrain

`funconstrain` collects the 35 nonlinear least-squares problems used by
Moré, Garbow, and Hillstrom to test optimization software. We will start
by fitting Rosenbrock with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), then see how to
choose dimensions and use the stored reference values.

## Solve your first problem

[`funconstrain_problem()`](https://jlmelville.github.io/funconstrain/reference/funconstrain_catalog.md)
returns a starting point and the functions needed by an optimizer:

``` r

library(funconstrain)

problem <- funconstrain_problem("rosen")
fit <- stats::optim(
  par = problem$x0,
  fn = problem$fn,
  gr = problem$gr,
  method = "BFGS"
)

final_objective <- problem$fn(fit$par)
gradient_inf_norm <- max(abs(problem$gr(fit$par)))

data.frame(
  convergence = fit$convergence,
  final_objective = signif(final_objective, 3),
  gradient_inf_norm = signif(gradient_inf_norm, 3)
)
#>   convergence final_objective gradient_inf_norm
#> 1           0        9.59e-18          2.48e-09
```

A convergence code of zero means that
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) reports
successful completion. The small objective and gradient norm provide a
second check on the returned point.

The three fields used by [`optim()`](https://rdrr.io/r/stats/optim.html)
have direct roles: `x0` is the starting vector, `fn` is the objective,
and `gr` is its gradient.

## Browse the catalogue

[`funconstrain_catalog()`](https://jlmelville.github.io/funconstrain/reference/funconstrain_catalog.md)
lists all 35 problems in paper order. Here `n` is the number of
parameters and `m` is the number of residual functions:

``` r

catalog <- funconstrain_catalog()
catalog_view <- catalog[
  catalog$name %in% c("rosen", "jenn_samp", "ex_rosen", "chebyquad"),
]

knitr::kable(data.frame(
  number = catalog_view$number,
  name = catalog_view$name,
  title = catalog_view$title,
  n = ifelse(
    catalog_view$n_kind == "fixed",
    as.character(catalog_view$n_default),
    paste0("variable; default ", catalog_view$n_default)
  ),
  m = ifelse(
    catalog_view$m_configurable,
    paste0("choose; default ", catalog_view$m_default),
    paste0("fixed or follows n; default ", catalog_view$m_default)
  ),
  row.names = NULL
))
```

| number | name | title | n | m |
|---:|:---|:---|:---|:---|
| 1 | rosen | Rosenbrock Function | 2 | fixed or follows n; default 2 |
| 6 | jenn_samp | Jennrich and Sampson Function | 2 | choose; default 10 |
| 21 | ex_rosen | Extended Rosenbrock Function | variable; default 8 | fixed or follows n; default 8 |
| 35 | chebyquad | Chebyquad Function | variable; default 50 | fixed or follows n; default 50 |

The [mathematical
definitions](https://jlmelville.github.io/funconstrain/articles/problem-definitions.md)
give the supported dimensions, standard start, and equations for every
catalogue entry.

## Choose `n` or `m`

For a variable-dimension problem, supply `n`. For a problem whose
residual count can be chosen independently, supply `m`:

``` r

extended <- funconstrain_problem("ex_rosen", n = 4)
sampled <- funconstrain_problem("jenn_samp", m = 20)

data.frame(
  problem = c(extended$problem$name, sampled$problem$name),
  supplied = c("n = 4", "m = 20"),
  n = c(
    extended$configuration$effective$n,
    sampled$configuration$effective$n
  ),
  m = c(
    extended$configuration$effective$m,
    sampled$configuration$effective$m
  )
)
#>     problem supplied n  m
#> 1  ex_rosen    n = 4 4  4
#> 2 jenn_samp   m = 20 2 20
```

Extended Rosenbrock derives `m = 4` from the chosen `n = 4`. Jennrich
and Sampson keeps its fixed `n = 2` and uses the requested `m = 20`.

## Check a reference before using it

Before comparing a result with a stored minimum, check that the
reference applies to the problem you resolved:

``` r

fmin_reference <- problem$reference$fmin
objective_gap <- if (identical(fmin_reference$status, "applicable")) {
  final_objective - fmin_reference$value
} else {
  NA_real_
}

data.frame(
  reference_status = fmin_reference$status,
  objective_gap = signif(objective_gap, 3)
)
#>   reference_status objective_gap
#> 1       applicable      9.59e-18
```

Minimum values and minimizer locations are checked separately. Extended
Rosenbrock at `n = 4` shows why:

``` r

data.frame(
  reference = c("Minimum value", "Minimizer location"),
  status = c(
    extended$reference$fmin$status,
    extended$reference$xmin$status
  ),
  value_stored = c(
    !is.null(extended$reference$fmin$value),
    !is.null(extended$reference$xmin$value)
  )
)
#>            reference         status value_stored
#> 1      Minimum value     applicable         TRUE
#> 2 Minimizer location not_applicable         TRUE
```

The minimum value zero applies at `n = 4`, while the stored minimizer
belongs to the package’s `n = 8` example. Check the status of the
particular reference you intend to use. The [`funconstrain_problem()`
reference](https://jlmelville.github.io/funconstrain/reference/funconstrain_problem.md)
documents the complete status and diagnostic vocabulary.

## Next steps

For equations and supported dimensions, continue with [Mathematical
definitions of the test
problems](https://jlmelville.github.io/funconstrain/articles/problem-definitions.md).
For optimizer examples, see [Using funconstrain with
mize](https://jlmelville.github.io/funconstrain/articles/mize-integration.md)
or [Using the optimx
bridge](https://jlmelville.github.io/funconstrain/articles/optimx-bridge.md).
