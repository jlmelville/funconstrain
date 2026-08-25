# Getting started with resolved problems

`funconstrain` provides the 35 nonlinear least-squares test problems
from Moré, Garbow, and Hillstrom. The resolver is the simplest entry
point for solver-neutral code: it turns a problem name and optional
dimensions into one concrete starting point, a set of dimension-checked
callbacks, and metadata describing the exact configuration.

## Find a problem

[`funconstrain_catalog()`](https://jlmelville.github.io/funconstrain/reference/funconstrain_catalog.md)
lists the available problems in paper order. Its columns distinguish
fixed and variable parameter dimensions (`n`) and show whether the
number of summand functions (`m`) can be configured.

``` r

library(funconstrain)

catalog <- funconstrain_catalog()
catalog[
  catalog$name %in% c("rosen", "jenn_samp", "ex_rosen", "chebyquad"),
  c(
    "number", "name", "n_kind", "n_default", "m_configurable", "m_default"
  )
]
#>    number      name   n_kind n_default m_configurable m_default
#> 1       1     rosen    fixed         2          FALSE         2
#> 6       6 jenn_samp    fixed         2           TRUE        10
#> 21     21  ex_rosen variable         8          FALSE         8
#> 35     35 chebyquad variable        50          FALSE        50
```

Names are matched exactly. This makes a selected problem unambiguous in
scripts and saved experiment configurations.

## Resolve one exact problem

The Rosenbrock problem has fixed `n = 2` and `m = 2`, so no dimensions
need to be supplied:

``` r

rosen_problem <- funconstrain_problem("rosen")
names(rosen_problem)
#> [1] "fn"            "gr"            "he"            "fg"           
#> [5] "x0"            "problem"       "configuration" "reference"
rosen_problem$x0
#> [1] -1.2  1.0
rosen_problem$configuration
#> $requested
#> $requested$n
#> NULL
#> 
#> $requested$m
#> NULL
#> 
#> 
#> $effective
#> $effective$n
#> [1] 2
#> 
#> $effective$m
#> [1] 2
```

The resolved callbacks are:

- `fn(x)`, the scalar objective;
- `gr(x)`, its gradient;
- `he(x)`, its Hessian; and
- `fg(x)`, a combined objective-and-gradient calculation.

Unlike the raw factory interface, a resolved `x0` is always a numeric
vector. Each callback is pinned to the effective value of `n`, so an
accidental vector of another length fails at the problem boundary.

## Run an optimizer

The callbacks can be passed directly to an optimizer. Here,
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) uses its BFGS
method and the analytic gradient:

``` r

fit <- stats::optim(
  par = rosen_problem$x0,
  fn = rosen_problem$fn,
  gr = rosen_problem$gr,
  method = "BFGS"
)

final_objective <- rosen_problem$fn(fit$par)
gradient_inf_norm <- max(abs(rosen_problem$gr(fit$par)))

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
successful completion. The final objective and gradient infinity norm
are evaluated again through the resolved problem, independently of the
optimizer’s stored objective value. Their small magnitudes provide a
useful quality check for this run; exact values and iteration counts can
vary across R versions and platforms.

## Request dimensions explicitly

For a variable-dimension problem, pass `n`. For a problem with
configurable summand count, pass `m`:

``` r

extended <- funconstrain_problem("ex_rosen", n = 4)
sampled <- funconstrain_problem("jenn_samp", m = 20)

data.frame(
  problem = c(extended$problem$name, sampled$problem$name),
  requested = c("n = 4", "m = 20"),
  effective_n = c(
    extended$configuration$effective$n,
    sampled$configuration$effective$n
  ),
  effective_m = c(
    extended$configuration$effective$m,
    sampled$configuration$effective$m
  )
)
#>     problem requested effective_n effective_m
#> 1  ex_rosen     n = 4           4           4
#> 2 jenn_samp    m = 20           2          20
```

`configuration$requested` records normalized values supplied by the
caller; an omitted value remains `NULL`. `configuration$effective`
records the concrete integer values used to create the starting point
and callbacks. Here, Extended Rosenbrock derives `m = 4` from the
requested `n = 4`, while Jennrich and Sampson combines the requested
`m = 20` with its default `n = 2`.

## Decide whether a reference applies

Each resolved problem has separate `fmin` and `xmin` reference records.
A record contains a `value`, a `status`, a machine-readable `reason`,
and `source_configuration`, the encoded configuration rule associated
with the record. The status says whether a value is available and
whether that rule establishes applicability. The table below is built
from four current problem configurations.

| Problem configuration | Reference | Status | Value stored | Source configuration |
|:---|:---|:---|:---|:---|
| rosen (defaults) | fmin | applicable | yes | all supported configurations |
| rosen (defaults) | xmin | applicable | yes | all supported configurations |
| ex_rosen (n = 4) | fmin | applicable | yes | all supported configurations |
| ex_rosen (n = 4) | xmin | not_applicable | yes | n = 8 |
| jenn_samp (m = 20) | fmin | not_applicable | yes | m = 10 |
| jenn_samp (m = 20) | xmin | not_applicable | yes | m = 10 |
| box_3d (defaults) | fmin | applicable | yes | all supported configurations |
| box_3d (defaults) | xmin | unavailable | no | all supported configurations |

Interpret the four possible statuses before using a stored value as a
target:

- `applicable`: the value applies to the effective configuration;
- `not_applicable`: the stored value is retained for inspection, but its
  source configuration does not match this problem instance;
- `unavailable`: the package has no stored value, so `value` is `NULL`;
  and
- `unknown`: a value is stored, but the package has not encoded enough
  information to decide applicability.

No current catalogue entry returns `unknown`, but it remains part of the
public record vocabulary so consumers do not have to treat missing
knowledge as either applicable or inapplicable. Notice also that `fmin`
and `xmin` can have different statuses: at `ex_rosen(n = 4)`, the
minimum value applies but the stored minimizer is for `n = 8`. Always
inspect the status of the specific reference you intend to use.

## When to use a raw factory

Direct factories such as
[`rosen()`](https://jlmelville.github.io/funconstrain/reference/rosen.md)
remain useful when you want to inspect the lower-level problem
definition or manage dimensions yourself. Their starting point is a
numeric vector for fixed-dimension problems and a function for
variable-dimension problems. The resolver is usually more convenient for
solver-neutral code because it returns one numeric start, pinned
callbacks, and explicit configuration and reference metadata.

See the [problem resolver
reference](https://jlmelville.github.io/funconstrain/reference/funconstrain_catalog.md)
for the complete return contract. To run the problems through the
optional `optimx` bridge, continue with [Using the optimx
bridge](https://jlmelville.github.io/funconstrain/articles/optimx-bridge.md).
